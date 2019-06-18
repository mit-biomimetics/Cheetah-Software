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


#ifndef CASADI_SHELL_INTERFACE_HPP
#define CASADI_SHELL_INTERFACE_HPP

#include "casadi/core/importer_internal.hpp"
#include <casadi/solvers/casadi_importer_shell_export.h>
#include "casadi/core/plugin_interface.hpp"

/** \defgroup plugin_Importer_shell
      Interface to the JIT compiler SHELL
*/

/** \pluginsection{Importer,shell} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{Importer,shell}


   \author Joel Andersson
   \date 2015
   *
   @copydoc Importer_doc
   @copydoc plugin_Importer_shell
   * */
  class CASADI_IMPORTER_SHELL_EXPORT ShellCompiler : public ImporterInternal {
  public:

    /** \brief Constructor */
    explicit ShellCompiler(const std::string& name);

    /** \brief  Create a new JIT function */
    static ImporterInternal* creator(const std::string& name) {
      return new ShellCompiler(name);
    }

    /** \brief Destructor */
    ~ShellCompiler() override;

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
    const char* plugin_name() const override { return "shell";}

    // Get name of the class
    std::string class_name() const override { return "ShellCompiler";}

    /// Get a function pointer for numerical evaluation
    signal_t get_function(const std::string& symname) override;
  protected:
    std::string base_name_;

    /// Temporary file
    std::string bin_name_;

    /// Temporary file
    std::string obj_name_;

    /// Extra files
    std::vector<std::string> extra_suffixes_;

    /// Cleanup temporary files when unloading
    bool cleanup_;

    // Shared library handle
    typedef DL_HANDLE_TYPE handle_t;
    handle_t handle_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SHELL_INTERFACE_HPP
