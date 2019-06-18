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


#include "importer_internal.hpp"

using namespace std;
namespace casadi {

  ImporterInternal::ImporterInternal(const std::string& name) : name_(name) {
    verbose_ = false;
  }

  ImporterInternal::~ImporterInternal() {
  }

  void ImporterInternal::disp(ostream &stream, bool more) const {
  }

  std::map<std::string, ImporterInternal::Plugin> ImporterInternal::solvers_;

  const std::string ImporterInternal::infix_ = "importer";

  Options ImporterInternal::options_
  = {{},
     {{"verbose",
       {OT_BOOL,
        "Verbose evaluation -- for debugging"}}
      }
    };

  void ImporterInternal::construct(const Dict& opts) {
    // Sanitize dictionary is needed
    if (!Options::is_sane(opts)) {
      // Call recursively
      return construct(Options::sanitize(opts));
    }

    // Make sure all options exist
    get_options().check(opts);

    // Initialize object
    init(opts);
  }

  void ImporterInternal::init(const Dict& opts) {
    // Read meta information from file
    if (can_have_meta()) {
      casadi_int offset = 0;
      ifstream file(name_);
      std::string line;
      while (getline(file, line)) {
        // Update offset
        offset++;

        // Try to find a /*CASADIMETA delimiter
        string cmd = "/*CASADIMETA";
        size_t pos = line.find(cmd);
        if (pos != string::npos) {
          read_meta(file, offset);
          continue;
        }

        // Try to find a /*CASADIEXTERNAL delimiter
        cmd = "/*CASADIEXTERNAL";
        pos = line.find(cmd);
        if (pos != string::npos) {
          istringstream ss(line.substr(pos+cmd.size()));
          // Read name
          string sym;
          ss >> sym;
          casadi_assert_dev(ss.good());
          // Default attributes
          bool inlined = false;

          // Read attributes: FIXME(@jaeandersson): Hacky
          size_t eqpos = line.find('=', pos+cmd.size());
          if (eqpos != string::npos) {
            string attr = "inline";
            if (line.compare(eqpos-attr.size(), attr.size(), attr)==0) {
              casadi_assert_dev(line.size()>eqpos+1);
              if (line.at(eqpos+1)=='1') {
                inlined=true;
              } else {
                casadi_assert_dev(line.at(eqpos+1)=='0');
              }
            }
          }

          read_external(sym, inlined, file, offset);
          continue;
        }
      }
    }

    // Read options
    for (auto&& op : opts) {
      if (op.first=="verbose") {
        verbose_ = op.second;
      }
    }
  }

  void ImporterInternal::read_meta(istream& file, casadi_int& offset) {
    // Loop over the lines
    std::string line;
    while (getline(file, line)) {
      offset++;

      // End of meta found?
      if (line.find("*/") != string::npos) return;

      // If comment or empty line, skip
      if (line.empty() || line.at(0)=='#') continue;

      // Get command string
      casadi_assert(line.at(0)==':',
                            "Syntax error: " + line + " is not a command string");
      string cmd = line.substr(1, line.find(' ')-1);

      // New entry
      stringstream ss;

      // Collect the meta data
      line = line.substr(cmd.size()+2);
      while (true) {
        // Find the backslash, if any
        size_t stop = line.find('\\');

        // Add to entry
        ss << line.substr(0, stop);

        // Break if not multiline
        if (stop == string::npos) break;

        // Read another line
        ss << std::endl;
        if (!getline(file, line)) {
          casadi_error("Failed to read \"" + cmd + "\"");
        }
        offset++;
      }

      // Insert new element in map
      auto new_el = meta_.insert(make_pair(cmd, make_pair(offset, ss.str())));
      casadi_assert(new_el.second, "Duplicate entry: \"" + cmd + "\"");
    }
    casadi_error("End-of-file reached while searching for \"*/\"");
  }

  void ImporterInternal::
  read_external(const string& sym, bool inlined, istream& file, casadi_int& offset) {
    // New entry
    stringstream ss;

    // Are we still in the function declaration
    bool in_declaration = true;

    // Loop over the lines
    std::string line;
    while (getline(file, line)) {
      offset++;

      // Skip line if still in declaration
      if (in_declaration) {
        size_t stop = line.find('{');
        if (stop != string::npos) in_declaration = false;
        continue;
      }

      // End of declaration found?
      if (line.find("/*CASADIEXTERNAL") != string::npos) {
        auto new_el = external_.insert(make_pair(sym, make_pair(inlined, ss.str())));
        casadi_assert(new_el.second, "Duplicate symbol: \"" + sym + "\"");
        return;
      }

      // Add to entry
      if (inlined) {
        ss << line << endl;
      }
    }
    casadi_error("End-of-file reached while searching for \"/*CASADIEXTERNAL\"");
  }

  bool ImporterInternal::has_function(const std::string& symname) const {
    // Check if in meta information
    if (external_.find(symname)!=external_.end()) return true;

    // Convert to a dummy function pointer
    return const_cast<ImporterInternal*>(this)->get_function(symname)!=nullptr;
  }

  DllLibrary::DllLibrary(const std::string& bin_name)
    : ImporterInternal(bin_name), handle_(nullptr) {
#ifdef WITH_DL
#ifdef _WIN32
    handle_ = LoadLibrary(TEXT(name_.c_str()));
    casadi_assert(handle_!=0,
      "CommonExternal: Cannot open \"" + name_ + "\". "
      "Error code (WIN32): " + str(GetLastError()));
#else // _WIN32
    handle_ = dlopen(name_.c_str(), RTLD_LAZY);
    casadi_assert(handle_!=nullptr,
      "CommonExternal: Cannot open \"" + name_ + "\". "
      "Error code: " + str(dlerror()));
    // reset error
    dlerror();
#endif // _WIN32
#else // WITH_DL
    casadi_error("CommonExternal: WITH_DL  not activated");
#endif // WITH_DL
  }

  DllLibrary::~DllLibrary() {
#ifdef WITH_DL
    // close the dll
#ifdef _WIN32
    if (handle_) FreeLibrary(handle_);
#else // _WIN32
    if (handle_) dlclose(handle_);
#endif // _WIN32
#endif // WITH_DL
  }

  signal_t DllLibrary::get_function(const std::string& sym) {
#ifdef WITH_DL
#ifdef _WIN32
    return (signal_t)GetProcAddress(handle_, TEXT(sym.c_str()));
#else // _WIN32
    signal_t fcnPtr = (signal_t)dlsym(handle_, sym.c_str());
    if (dlerror()) {
      fcnPtr=nullptr;
      dlerror(); // Reset error flags
    }
    return fcnPtr;
#endif // _WIN32
#endif // WITH_DL
  }

  std::string ImporterInternal::get_meta(const std::string& cmd, casadi_int ind) const {
    if (ind>=0) return get_meta(indexed(cmd, ind));
    casadi_assert(has_meta(cmd), "No such command: " + cmd);
    return meta_.at(cmd).second;
  }

  bool ImporterInternal::has_meta(const std::string& cmd, casadi_int ind) const {
    if (ind>=0) return has_meta(indexed(cmd, ind));
    return meta_.find(cmd) != meta_.end();
  }

  bool ImporterInternal::inlined(const std::string& symname) const {
    auto it = external_.find(symname);
    return it!=external_.end() && it->second.first;
  }

  std::string ImporterInternal::body(const std::string& symname) const {
    auto it = external_.find(symname);
    casadi_assert_dev(it!=external_.end() && it->second.first);
    return it->second.second;
  }

} // namespace casadi
