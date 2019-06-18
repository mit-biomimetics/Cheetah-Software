/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef CASADI_PLUGIN_INTERFACE_HPP
#define CASADI_PLUGIN_INTERFACE_HPP

#include "function_internal.hpp"
#include "global_options.hpp"

#include <stdlib.h>

/// \cond INTERNAL

// For dynamic loading
#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32

// Set default shared library prefix
#ifndef SHARED_LIBRARY_PREFIX
#define SHARED_LIBRARY_PREFIX "lib"
#endif // SHARED_LIBRARY_PREFIX

// Set default shared library suffix
#ifndef SHARED_LIBRARY_SUFFIX
#define SHARED_LIBRARY_SUFFIX ".so"
#endif // SHARED_LIBRARY_SUFFIX

#ifdef _WIN32
    #define DL_HANDLE_TYPE HINSTANCE
#else // _WIN32
    #define DL_HANDLE_TYPE void *
#endif

// http://stackoverflow.com/questions/303562/c-format-macro-inline-ostringstream
#define STRING(ITEMS) \
  ((dynamic_cast<std::ostringstream &>(std::ostringstream() \
    . seekp(0, std::ios_base::cur) << ITEMS)) . str())

#endif // WITH_DL

namespace casadi {
  // Avoid segmentation faults when exposed function not implemented
  template<typename T>
  T check_exposed(T t) {
    casadi_assert(t!=0, "Static function not implemented for plugin");
    return t;
  }

  /** \brief Interface for accessing input and output data structures
      \author Joel Andersson
      \date 2013
  */
  template<class Derived>
  class PluginInterface {
    public:

    /// Fields
    struct Plugin{
      typename Derived::Creator creator;
      const char* name;
      const char* doc;
      int version;
      typename Derived::Exposed exposed;
      const Options* options;
      // Constructor
      Plugin() : creator(nullptr), name(nullptr), doc(nullptr), version(0), options(nullptr) {}
    };

    // Plugin registration function
    typedef int (*RegFcn)(Plugin* plugin);

    /// Check if a plugin is available or can be loaded
    static bool has_plugin(const std::string& pname, bool verbose=false);

    /// Get the plugin options
    static const Options& plugin_options(const std::string& pname);

    /// Instantiate a Plugin struct from a factory function
    static Plugin pluginFromRegFcn(RegFcn regfcn);

    /// Load a plugin dynamically
    static Plugin load_plugin(const std::string& pname, bool register_plugin=true);

    /// Load a library dynamically
    static DL_HANDLE_TYPE load_library(const std::string& libname, std::string& resultpath,
      bool global);

    /// Register an integrator in the factory
    static void registerPlugin(const Plugin& plugin);

    /// Register an integrator in the factory
    static void registerPlugin(RegFcn regfcn);

    /// Load and get the creator function
    static Plugin& getPlugin(const std::string& pname);

    // Create solver instance
    template<class Problem>
      static Derived* instantiate(const std::string& fname,
                                        const std::string& pname, Problem problem);
    // Get name of the plugin
    virtual const char* plugin_name() const = 0;
  };

  template<class Derived>
  bool PluginInterface<Derived>::has_plugin(const std::string& pname, bool verbose) {

    // Quick return if available
    if (Derived::solvers_.find(pname) != Derived::solvers_.end()) {
      return true;
    }

    // Try loading the plugin
    try {
      (void)load_plugin(pname, false);
      return true;
    } catch (CasadiException& ex) {
      if (verbose) {
        casadi_warning(ex.what());
      }
      return false;
    }
  }

  template<class Derived>
  const Options& PluginInterface<Derived>::plugin_options(const std::string& pname) {
    const Options *op = getPlugin(pname).options;
    casadi_assert(op!=nullptr, "Plugin \"" + pname + "\" does not support options");
    return *op;
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin
      PluginInterface<Derived>::pluginFromRegFcn(RegFcn regfcn) {
    // Create a temporary struct
    Plugin plugin;

    // Set the fields
    int flag = regfcn(&plugin);
    casadi_assert_dev(flag==0);

    return plugin;
  }


  template<class Derived>
  DL_HANDLE_TYPE PluginInterface<Derived>::load_library(const std::string& libname,
    std::string& resultpath, bool global) {

#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL

    // Get the name of the shared library
    std::string lib = SHARED_LIBRARY_PREFIX + libname + SHARED_LIBRARY_SUFFIX;

    // Build up search paths;
    std::vector<std::string> search_paths;

    #ifdef _WIN32
    char pathsep = ';';
    const std::string filesep("\\");
    #else
    char pathsep = ':';
    const std::string filesep("/");
    #endif

    // Search path: global casadipath option
    std::stringstream casadipaths(GlobalOptions::getCasadiPath());
    std::string casadipath;
    while (std::getline(casadipaths, casadipath, pathsep)) {
      search_paths.push_back(casadipath);
    }

    // Search path: CASADIPATH env variable
    char* pLIBDIR;
    pLIBDIR = getenv("CASADIPATH");

    if (pLIBDIR!=nullptr) {
      std::stringstream casadipaths(pLIBDIR);
      std::string casadipath;
      while (std::getline(casadipaths, casadipath, pathsep)) {
        search_paths.push_back(casadipath);
      }
    }

    // Search path: bare
    search_paths.push_back("");

    // Search path : PLUGIN_EXTRA_SEARCH_PATH
    #ifdef PLUGIN_EXTRA_SEARCH_PATH
    search_paths.push_back(
      std::string("") + PLUGIN_EXTRA_SEARCH_PATH);
    #endif // PLUGIN_EXTRA_SEARCH_PATH

    // Search path : current directory
    search_paths.push_back(".");

    // Prepare error string
    std::stringstream errors;
    errors << "PluginInterface::load_plugin: Cannot load shared library '"
           << lib << "': " << std::endl;
    errors << "   (\n"
           << "    Searched directories: 1. casadipath from GlobalOptions\n"
           << "                          2. CASADIPATH env var\n"
           << "                          3. PATH env var (Windows)\n"
           << "                          4. LD_LIBRARY_PATH env var (Linux)\n"
           << "                          5. DYLD_LIBRARY_PATH env var (osx)\n"
           << "    A library may be 'not found' even if the file exists:\n"
           << "          * library is not compatible (different compiler/bitness)\n"
           << "          * the dependencies are not found\n"
           << "   )";

    // Alocate a handle pointer
#ifdef _WIN32
    HINSTANCE handle = 0;
#else // _WIN32
    void * handle = nullptr;
#endif

    // Alocate a handle pointer
#ifndef _WIN32
    int flag;
    if (global) {
      flag = RTLD_NOW | RTLD_GLOBAL;
    } else {
      flag = RTLD_LAZY | RTLD_LOCAL;
    }
#ifdef WITH_DEEPBIND
#ifndef __APPLE__
    flag |= RTLD_DEEPBIND;
#endif
#endif
#endif

    std::string searchpath;

    // Try getting a handle
    for (casadi_int i=0;i<search_paths.size();++i) {
      searchpath = search_paths[i];
#ifdef _WIN32
      SetDllDirectory(TEXT(searchpath.c_str()));
      handle = LoadLibrary(TEXT(lib.c_str()));
      SetDllDirectory(NULL);
#else // _WIN32
      std::string libname = searchpath.size()==0 ? lib : searchpath + filesep + lib;
      handle = dlopen(libname.c_str(), flag);
#endif // _WIN32
      if (handle) {
        break;
      } else {
        errors << std::endl << "  Tried '" << searchpath << "' :";
#ifdef _WIN32
        errors << std::endl << "    Error code (WIN32): " << STRING(GetLastError());
#else // _WIN32
        errors << std::endl << "    Error code: " << dlerror();
#endif // _WIN32
      }
    }

    resultpath = searchpath;

    casadi_assert(handle!=nullptr, errors.str());

    return handle;

#endif // WITH_DL
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin
      PluginInterface<Derived>::load_plugin(const std::string& pname, bool register_plugin) {
    // Issue warning and quick return if already loaded
    if (Derived::solvers_.find(pname) != Derived::solvers_.end()) {
      casadi_warning("PluginInterface: Solver " + pname + " is already in use. Ignored.");
      return Plugin();
    }


#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL
    // Retrieve the registration function
    RegFcn reg;

    // Load the dll
    std::string regName = "casadi_register_" + Derived::infix_ + "_" + pname;

    std::string searchpath;
    DL_HANDLE_TYPE handle = load_library("casadi_" + Derived::infix_ + "_" + pname, searchpath,
      false);

#ifdef _WIN32
    reg = (RegFcn)GetProcAddress(handle, TEXT(regName.c_str()));
#else // _WIN32
    // Reset error
    dlerror();

    // Load creator
    reg = (RegFcn)dlsym(handle, regName.c_str());
#endif // _WIN32
    casadi_assert(reg!=nullptr,
      "PluginInterface::load_plugin: no \"" + regName + "\" found in " + searchpath + ".");

    // Create a temporary struct
    Plugin plugin = pluginFromRegFcn(reg);
    // Register the plugin
    if (register_plugin) {
      registerPlugin(plugin);
    }

    return plugin;

#endif // WITH_DL
  }

  template<class Derived>
  void PluginInterface<Derived>::registerPlugin(RegFcn regfcn) {
    registerPlugin(pluginFromRegFcn(regfcn));
  }

  template<class Derived>
  void PluginInterface<Derived>::registerPlugin(const Plugin& plugin) {

    // Check if the solver name is in use
    typename std::map<std::string, Plugin>::iterator it=Derived::solvers_.find(plugin.name);
    casadi_assert(it==Derived::solvers_.end(),
      "Solver " + str(plugin.name) + " is already in use");

    // Add to list of solvers
    Derived::solvers_[plugin.name] = plugin;
  }

  template<class Derived>
  typename PluginInterface<Derived>::Plugin&
  PluginInterface<Derived>::getPlugin(const std::string& pname) {

    // Check if the solver has been loaded
    auto it=Derived::solvers_.find(pname);

    // Load the solver if needed
    if (it==Derived::solvers_.end()) {
      load_plugin(pname);
      it=Derived::solvers_.find(pname);
    }
    casadi_assert_dev(it!=Derived::solvers_.end());
    return it->second;
  }

  template<class Derived>
  template<class Problem>
  Derived* PluginInterface<Derived>::
  instantiate(const std::string& fname,
              const std::string& pname, Problem problem) {

    // Assert the plugin exists (needed for adaptors)
    if (!has_plugin(pname, true)) {
      casadi_error("Plugin '" + pname + "' is not found.");
    }
    return getPlugin(pname).creator(fname, problem);
  }

} // namespace casadi

/// \endcond

#endif // CASADI_PLUGIN_INTERFACE_HPP
