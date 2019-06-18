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


#ifndef CASADI_IMPORTER_HPP
#define CASADI_IMPORTER_HPP

#include "function.hpp"
#include "printable.hpp"

namespace casadi {

#ifndef SWIG
    /** \brief Convert to a type */
    template<typename T>
    T text2type(const std::string& text) {
      std::istringstream ss(text);
      T ret;
      ss >> ret;
      return ret;
    }

    /** \brief Get entry as a vector */
    template<typename T>
    std::vector<T> text2vector(const std::string& text) {
      std::istringstream ss(text);
      T tmp;
      std::vector<T> ret;
      while (ss >> tmp) ret.push_back(tmp);
      return ret;
    }

    /** \brief Get entry as a set */
    template<typename T>
    std::set<T> text2set(const std::string& text) {
      std::vector<T> v = text2vector<T>(text);
      std::set<T> ret;
      for (auto&& i : v) ret.insert(i);
      return ret;
    }
#endif // SWIG

  // Forward declaration of internal class
  class ImporterInternal;

  /** \brief Importer

      Just-in-time compilation of code

      \generalsection{Importer}
      \pluginssection{Importer}

      \author Joris Gillis
      \date 2015
  */
  class CASADI_EXPORT Importer
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Importer>) {
  public:
    /** \brief Get type name */
    static std::string type_name() {return "Importer";}

    /// Default constructor
    Importer();

    /// Importer factory
    explicit Importer(const std::string& name,
                      const std::string& compiler,
                      const Dict& opts=Dict());

    /// Access functions of the node
    ImporterInternal* operator->();
    const ImporterInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /// Check if a plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    // Check if symbol exists
    bool has_function(const std::string& symname) const;

#ifndef SWIG
    /// Get a function pointer for numerical evaluation
    signal_t get_function(const std::string& symname);
#endif // SWIG

    /** \brief Does a meta entry exist? */
    bool has_meta(const std::string& cmd, casadi_int ind=-1) const;

    /** \brief Get entry as a text */
    std::string get_meta(const std::string& cmd, casadi_int ind=-1) const;

    /// Check if a function is inlined
    bool inlined(const std::string& symname) const;

    /// Get the function body, if inlined
    std::string body(const std::string& symname) const;

#ifndef SWIG
    /** Convert indexed command */
    static inline std::string indexed(const std::string& cmd, casadi_int ind) {
      std::stringstream ss;
      ss << cmd << "[" << ind << "]";
      return ss.str();
    }

    /** \brief Convert to a type */
    template<typename T>
    T to(const std::string& cmd, casadi_int ind=-1) const {
      return text2type<T>(get_meta(cmd, ind));
    }

    /** \brief Get entry as a string */
    std::string meta_string(const std::string& cmd, casadi_int ind=-1) const {
      return to<std::string>(cmd, ind);
    }

    /** \brief Get entry as a vector */
    template<typename T>
    std::vector<T> meta_vector(const std::string& cmd, casadi_int ind=-1) const {
      return text2vector<T>(get_meta(cmd, ind));
    }

    /** \brief Get entry as a set */
    template<typename T>
    std::set<T> meta_set(const std::string& cmd, casadi_int ind=-1) const {
      return text2set<T>(get_meta(cmd, ind));
    }

    /** \brief Get entry as an integer */
    casadi_int meta_int(const std::string& cmd, casadi_int ind=-1) const {
      return to<casadi_int>(cmd, ind);
    }
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_IMPORTER_HPP
