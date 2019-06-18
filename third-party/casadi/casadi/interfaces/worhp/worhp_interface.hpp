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


#ifndef CASADI_WORHP_INTERFACE_HPP
#define CASADI_WORHP_INTERFACE_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/interfaces/worhp/casadi_nlpsol_worhp_export.h>

// Workaround for Clang, but should not be a problem for other compilers, #771
#define _Bool bool

#include <worhp.h>

// MACROs that pollute our code
#undef Q
/**\defgroup plugin_Nlpsol_worhp
 WORHP interface

 Designed for Worhp 1.12

*/
/** \pluginsection{Nlpsol,worhp} **/

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_WORHP_EXPORT WorhpMemory : public NlpsolMemory {
    OptVar    worhp_o;
    Workspace worhp_w;
    Params    worhp_p;
    Control   worhp_c;

    // Stats
    casadi_int iter;
    casadi_int iter_sqp;
    double inf_pr;
    double inf_du;
    double alpha_pr;
    casadi_int return_code;
    const char* return_status;

    bool init_;

    /// Constructor
    WorhpMemory();

    /// Destructor
    ~WorhpMemory();
  };

  /** \brief \pluginbrief{Nlpsol,worhp}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_worhp
  */
  class CASADI_NLPSOL_WORHP_EXPORT WorhpInterface : public Nlpsol {
  public:
    // NLP functions
    Function f_fcn_;
    Function g_fcn_;
    Function grad_f_fcn_;
    Function jac_g_fcn_;
    Function hess_l_fcn_;
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;

    // Constructor
    explicit WorhpInterface(const std::string& name, const Function& nlp);

    // Destructor
    ~WorhpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "worhp";}

    // Get name of the class
    std::string class_name() const override { return "WorhpInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new WorhpInterface(name, nlp);
    }

    // Reset solver
    void reset();

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new WorhpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<WorhpMemory*>(mem);}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    // Options
    std::map<std::string, bool> bool_opts_;
    std::map<std::string, casadi_int> int_opts_;
    std::map<std::string, double> double_opts_;
    Dict qp_opts_;

    // WORHP return codes
    static const char* return_codes(casadi_int flag);

    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_WORHP_INTERFACE_HPP
