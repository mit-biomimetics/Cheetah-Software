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


#ifndef CASADI_CLANG_COMPILER_HPP
#define CASADI_CLANG_COMPILER_HPP

#include "casadi/core/importer_internal.hpp"
#include <casadi/interfaces/clang/casadi_importer_clang_export.h>

#include <clang/CodeGen/CodeGenAction.h>
#include <clang/Basic/DiagnosticOptions.h>
#include <clang/Basic/TargetInfo.h>
#include <clang/Basic/Version.h>
#include <clang/Driver/Compilation.h>
#include <clang/Driver/Driver.h>
#include <clang/Driver/Tool.h>
#include <clang/Frontend/CompilerInstance.h>
#include <clang/Frontend/CompilerInvocation.h>
#include <clang/Frontend/FrontendDiagnostic.h>
#include <clang/Frontend/TextDiagnosticPrinter.h>

#include <llvm/ADT/SmallString.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#if (LLVM_VERSION_MAJOR>=4) || (LLVM_VERSION_MAJOR==3 && LLVM_VERSION_MINOR>=5)
#include <llvm/ExecutionEngine/MCJIT.h>
#else
#include "llvm/ExecutionEngine/JIT.h"
#endif
#include <llvm/IR/Module.h>
#include "llvm/IR/LLVMContext.h"
//#include "llvm/IR/Verifier.h"
#include <llvm/Support/FileSystem.h>
#include <llvm/Support/Host.h>
#include <llvm/Support/ManagedStatic.h>
#include <llvm/Support/Path.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/raw_os_ostream.h>
//#include <llvm/ExecutionEngine/ExecutionEngine.h>

/** \defgroup plugin_Importer_clang
      Interface to the JIT compiler CLANG
*/

/** \pluginsection{Importer,clang} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{Importer,clang}


   \author Joris Gillis
   \date 2015
   *
   @copydoc Importer_doc
   @copydoc plugin_Importer_clang
   * */
  class CASADI_IMPORTER_CLANG_EXPORT ClangCompiler : public ImporterInternal {
  public:

    /** \brief Constructor */
    explicit ClangCompiler(const std::string& name);

    /** \brief  Create a new JIT function */
    static ImporterInternal* creator(const std::string& name) {
      return new ClangCompiler(name);
    }

    /** \brief Destructor */
    ~ClangCompiler() override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /// A documentation string
    static const std::string meta_doc;

    /// Get name of plugin
    const char* plugin_name() const override { return "clang";}

    // Get name of the class
    std::string class_name() const override { return "ClangCompiler";}

    /// Get a function pointer for numerical evaluation
    signal_t get_function(const std::string& symname) override;

    // Helper function for reading includes
    static std::vector<std::pair<std::string, bool> >
      getIncludes(const std::string& file, const std::string& path);

    // Options
    std::string include_path_;
    std::vector<std::string> flags_;

  protected:
    clang::EmitLLVMOnlyAction* act_;
    llvm::ExecutionEngine* executionEngine_;
    llvm::LLVMContext* context_;
    llvm::raw_ostream* myerr_;
    llvm::Module* module_; // owned by executionEngine_
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CLANG_COMPILER_HPP
