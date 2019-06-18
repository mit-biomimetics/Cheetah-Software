/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "clang_compiler.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/casadi_meta.hpp"
#include <fstream>

// To be able to get the plugin path
#ifdef _WIN32 // also for 64-bit
#define NOMINMAX
#include <windows.h>
#include <shlwapi.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_IMPORTER_CLANG_EXPORT
  casadi_register_importer_clang(ImporterInternal::Plugin* plugin) {
    plugin->creator = ClangCompiler::creator;
    plugin->name = "clang";
    plugin->doc = ClangCompiler::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ClangCompiler::options_;
    return 0;
  }

  extern "C"
  void CASADI_IMPORTER_CLANG_EXPORT casadi_load_importer_clang() {
    ImporterInternal::registerPlugin(casadi_register_importer_clang);
  }

  ClangCompiler::ClangCompiler(const std::string& name) :
    ImporterInternal(name) {

    myerr_ = nullptr;
    executionEngine_ = nullptr;
    context_ = nullptr;
    act_ = nullptr;
  }

  ClangCompiler::~ClangCompiler() {
    if (act_) delete act_;
    if (myerr_) delete myerr_;
    if (executionEngine_) delete executionEngine_;
    if (context_) delete context_;
  }

  Options ClangCompiler::options_
  = {{&ImporterInternal::options_},
     {{"include_path",
       {OT_STRING,
        "Include paths for the JIT compiler. "
        "The include directory shipped with CasADi will be automatically appended."}},
      {"flags",
       {OT_STRINGVECTOR,
        "Compile flags for the JIT compiler. Default: None"}}
     }
  };

  void ClangCompiler::init(const Dict& opts) {
    // Base class
    ImporterInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="include_path") {
        include_path_ = op.second.to_string();
      } else if (op.first=="flags") {
        flags_ = op.second;
      }
    }

    // Arguments to pass to the clang frontend
    vector<const char *> args(1, name_.c_str());
    for (auto&& f : flags_) {
      args.push_back(f.c_str());
    }

    // Create the compiler instance
    clang::CompilerInstance compInst;

    // A symbol in the DLL
    void *addr = reinterpret_cast<void*>(&casadi_register_importer_clang);

    // Get runtime include path
    std::string jit_include, filesep;
#ifdef _WIN32
    char buffer[MAX_PATH];
    HMODULE hm = NULL;
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR)addr, &hm)) {
      casadi_error("GetModuleHandle failed");
    }
    GetModuleFileNameA(hm, buffer, sizeof(buffer));
    PathRemoveFileSpecA(buffer);
    jit_include = buffer;
    filesep = "\\";
#else // _WIN32
    Dl_info dl_info;
    if (!dladdr(addr, &dl_info)) {
      casadi_error("dladdr failed");
    }
    jit_include = dl_info.dli_fname;
    jit_include = jit_include.substr(0, jit_include.find_last_of('/'));
    filesep = "/";
#endif // _WIN32
    jit_include += filesep + "casadi" + filesep + "jit";

#if 0
    // Initialize target info with the default triple for our platform.
    auto targetoptions = std::make_shared<clang::Taroptions>();
    targetoptions->Triple = llvm::sys::getDefaultTargetTriple();
    clang::TargetInfo *targetInfo =
      clang::TargetInfo::CreateTargetInfo(compInst.get_diagnostics(), targetoptions);
    compInst.setTarget(targetInfo);
#endif

    // The compiler invocation needs a DiagnosticsEngine so it can report problems
    clang::DiagnosticOptions* diagOpts = new clang::DiagnosticOptions();
    myerr_ = new llvm::raw_os_ostream(uerr());
    clang::TextDiagnosticPrinter *diagClient = new clang::TextDiagnosticPrinter(*myerr_, diagOpts);

    clang::DiagnosticIDs* diagID = new clang::DiagnosticIDs();
    // This object takes ownerships of all three passed-in pointers
    clang::DiagnosticsEngine diags(diagID, diagOpts, diagClient);

    // Create the compiler invocation
    #if LLVM_VERSION_MAJOR>=4
    std::shared_ptr<clang::CompilerInvocation> compInv(new clang::CompilerInvocation());
    #else
    clang::CompilerInvocation* compInv = new clang::CompilerInvocation();
    #endif
    clang::CompilerInvocation::CreateFromArgs(*compInv, &args[0],
                                              &args[0] + args.size(), diags);
    compInst.setInvocation(compInv);

    // Get ready to report problems
    compInst.createDiagnostics();
    if (!compInst.hasDiagnostics())
      casadi_error("Cannot create diagnostics");

    // Set resource directory
    std::string resourcedir = jit_include + filesep + "clang" + filesep + CLANG_VERSION_STRING;
    compInst.getHeaderSearchOpts().ResourceDir = resourcedir;

    // Read the system includes (C or C++)
    vector<pair<string, bool> > system_include = getIncludes("system_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::System, i->second, false);
    }

    // Read the system includes (C only)
    system_include = getIncludes("csystem_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::CSystem, i->second, false);
    }

    // Read the system includes (C++ only)
    system_include = getIncludes("cxxsystem_includes.txt", jit_include);
    for (auto i=system_include.begin(); i!=system_include.end(); ++i) {
      compInst.getHeaderSearchOpts().AddPath(i->first,
                                             clang::frontend::CXXSystem, i->second, false);
    }

#ifdef _WIN32
    char pathsep = ';';
#else
    char pathsep = ':';
#endif

    // Search path
    std::stringstream paths;
    paths << include_path_ << pathsep;
    std::string path;
    while (std::getline(paths, path, pathsep)) {
      compInst.getHeaderSearchOpts().AddPath(path.c_str(), clang::frontend::System, false, false);
    }

    // Create an LLVM context (NOTE: should use a static context instead?)
    context_ = new llvm::LLVMContext();

    // Create an action and make the compiler instance carry it out
    act_ = new clang::EmitLLVMOnlyAction(context_);
    if (!compInst.ExecuteAction(*act_))
      casadi_error("Cannot execute action");

    // Grab the module built by the EmitLLVMOnlyAction
    #if LLVM_VERSION_MAJOR>=4 || (LLVM_VERSION_MAJOR==3 && LLVM_VERSION_MINOR>=5)
    std::unique_ptr<llvm::Module> module = act_->takeModule();
    module_ = module.get();
    #else
    llvm::Module* module = act_->takeModule();
    module_ = module;
    #endif

    llvm::InitializeNativeTarget();
    llvm::InitializeNativeTargetAsmPrinter();

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    executionEngine_ =
      llvm::EngineBuilder(std::move(module)).setEngineKind(llvm::EngineKind::JIT)
      .setErrorStr(&ErrStr).create();
    if (!executionEngine_) {
      casadi_error("Could not create ExecutionEngine: " + ErrStr);
    }

    executionEngine_->finalizeObject();
  }

  signal_t ClangCompiler::get_function(const std::string& symname) {
    llvm::Function* f = module_->getFunction(symname);
    if (f) {
      return reinterpret_cast<signal_t>(executionEngine_->getPointerToFunction(f));
    } else {
      return nullptr;
    }
  }

  std::vector<std::pair<std::string, bool> > ClangCompiler::
  getIncludes(const std::string& file, const std::string& path) {
    // File separator
#ifdef _WIN32
    const char sep = '\\';
#else // _WIN32
    const char sep = '/';
#endif // _WIN32

    // Return value
    vector<pair<string, bool> > ret;

    // Read line-by-line
    std::ifstream setup_file(path + sep + file);
    std::string line;
    while (std::getline(setup_file, line)) {
      // Skip empty lines
      if (line.empty()) continue;

      // Check if framework
      size_t loc = line.find(" (framework directory)");
      bool isframework = loc != string::npos;
      if (isframework) {
        // Truncate path
        line = line.substr(0, loc);
      }

      // Check if the path is absolute or relative
#ifdef _WIN32
      bool relative = PathIsRelative(TEXT(line.c_str()));
#else // _WIN32
      bool relative = line.at(0)!=sep;
#endif // _WIN32

      if (relative) {
        // Relative path, make absolute
        ret.push_back(make_pair(path + sep + line, isframework));
      } else {
        // Absolute path
        ret.push_back(make_pair(line, isframework));
      }
    }

    return ret;
  }

} // namespace casadi
