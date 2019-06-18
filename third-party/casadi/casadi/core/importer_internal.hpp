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


#ifndef CASADI_IMPORTER_INTERNAL_HPP
#define CASADI_IMPORTER_INTERNAL_HPP

#include "importer.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief Importer internal class

  @copydoc Importer_doc
  \author Joel Andersson
  \date 2010-2013
*/
  class CASADI_EXPORT
  ImporterInternal : public SharedObjectInternal,
                     public PluginInterface<ImporterInternal> {

  public:
    /// Constructor
    explicit ImporterInternal(const std::string& name);

    /// Destructor
    ~ImporterInternal() override;

    /** \brief Get type name */
    std::string class_name() const override { return "ImporterInternal";}

    /** \brief Print */
    void disp(std::ostream& stream, bool more) const override;

    // Creator function for internal class
    typedef ImporterInternal* (*Creator)(const std::string& name);

    /** \brief Construct
        Prepares the function for evaluation
     */
    void construct(const Dict& opts);

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief Initialize */
    virtual void init(const Dict& opts);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "importer";}

    /// Queery plugin name
    const char* plugin_name() const override { return "none";}

    /// Get a function pointer for numerical evaluation
    virtual signal_t get_function(const std::string& symname) { return nullptr;}

    /// Get a function pointer for numerical evaluation
    bool has_function(const std::string& symname) const;

    /** \brief Does an entry exist? */
    bool has_meta(const std::string& cmd, casadi_int ind=-1) const;

    /** \brief Get entry as a text */
    std::string get_meta(const std::string& cmd, casadi_int ind=-1) const;

    /// Get meta information
    void read_meta(std::istream& file, casadi_int& offset);

    /// Get an external function declaration
    void read_external(const std::string& sym, bool inlined,
                       std::istream& file, casadi_int& offset);

    // Check if a function is inlined
    bool inlined(const std::string& symname) const;

    /// Get the function body, if inlined
    std::string body(const std::string& symname) const;

    /// Can meta information be read?
    virtual bool can_have_meta() const { return true;}

    /** \brief Get entry as a text */
    std::string to_text(const std::string& cmd, casadi_int ind=-1) const;

    /** Convert indexed command */
    static inline std::string indexed(const std::string& cmd, casadi_int ind) {
      std::stringstream ss;
      ss << cmd << "[" << ind << "]";
      return ss.str();
    }

    /// C filename
    std::string name_;

    /// Meta data
    std::map<std::string, std::pair<casadi_int, std::string> > meta_;

    /// External functions
    std::map<std::string, std::pair<bool, std::string> > external_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;
  };

  /** \brief Dynamically linked library
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT
  DllLibrary : public ImporterInternal {
  private:
#if defined(WITH_DL) && defined(_WIN32) // also for 64-bit
    typedef HINSTANCE handle_t;
#else
    typedef void* handle_t;
#endif
    handle_t handle_;
  public:

    // Constructor
    explicit DllLibrary(const std::string& bin_name);

    // Destructor
    ~DllLibrary() override;

    /** \brief Get type name */
    std::string class_name() const override { return "DllLibrary";}

    // Dummy type
    signal_t get_function(const std::string& symname) override;

    /// Can meta information be read?
    bool can_have_meta() const override { return false;}
  };

} // namespace casadi
/// \endcond
#endif // CASADI_IMPORTER_INTERNAL_HPP
