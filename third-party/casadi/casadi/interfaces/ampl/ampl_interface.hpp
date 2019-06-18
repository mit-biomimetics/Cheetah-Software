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


#ifndef CASADI_AMPL_INTERFACE_HPP
#define CASADI_AMPL_INTERFACE_HPP

#include <casadi/interfaces/ampl/casadi_nlpsol_ampl_export.h>
#include "casadi/core/nlpsol_impl.hpp"

/** \defgroup plugin_Nlpsol_ampl
  * Interface to AMPL
  *
  * \author Joel Andersson
  * \date 2017
*/

/** \pluginsection{Nlpsol,AmplInterface} */

/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class AmplInterface;

  struct CASADI_NLPSOL_AMPL_EXPORT AmplInterfaceMemory : public NlpsolMemory {
  };

  /** \brief \pluginbrief{Nlpsol,AmplInterface}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_AmplInterface
  */
  class CASADI_NLPSOL_AMPL_EXPORT AmplInterface : public Nlpsol {
  public:
    explicit AmplInterface(const std::string& name, const Function& nlp);
    ~AmplInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "ampl";}

    // Get name of the class
    std::string class_name() const override { return "AmplInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new AmplInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new AmplInterfaceMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<AmplInterfaceMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Solver binary
    std::string solver_;

    // Construction of the NL problem
    std::stringstream nl_init_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_AMPL_INTERFACE_HPP
