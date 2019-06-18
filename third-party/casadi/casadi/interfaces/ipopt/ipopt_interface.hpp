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


#ifndef CASADI_IPOPT_INTERFACE_HPP
#define CASADI_IPOPT_INTERFACE_HPP

#include <casadi/interfaces/ipopt/casadi_nlpsol_ipopt_export.h>
#include "casadi/core/nlpsol_impl.hpp"
#include "casadi/core/timing.hpp"

/** \defgroup plugin_Nlpsol_ipopt
 *
 * When in warmstart mode, output NLPSOL_LAM_X may be used as input
 *
 * NOTE: Even when max_iter == 0, it is not guaranteed that
 * input(NLPSOL_X0) == output(NLPSOL_X).
 * Indeed if bounds on X or constraints are unmet, they will differ.
 *
 * For a good tutorial on IPOPT, see
 * http://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf
 *
 * A good resource about the algorithms in IPOPT is: Wachter and L. T. Biegler,
 * On the Implementation of an Interior-Point Filter Line-Search Algorithm for
 * Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57,
 * 2006 (As Research Report RC 23149, IBM T. J. Watson Research Center, Yorktown, USA
 *
 * Caveats:
 * * with default options, multipliers for the decision variables are wrong for equality
 * constraints.
 * Change the 'fixed_variable_treatment' to 'make_constraint' or 'relax_bounds' to obtain
 * correct results.
 *
 */

/** \pluginsection{Nlpsol,ipopt} **/

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_IPOPT_EXPORT IpoptMemory : public NlpsolMemory {
    /** NOTE:
     * To allow this header file to be free of IPOPT types
     * (that are sometimes declared outside their scope!) and after
     * experiencing problems with working with IPOPT classes without
     * IPOPT smart pointers, we work with dynamically allocated IPOPT
     * smart pointers in this interface, that are stored as void
     * pointers in the interface.
     *
     */
    void* userclass;
    void* app;

    // Current calculated quantities
    double *gk, *grad_fk, *jac_gk, *hess_lk, *grad_lk;

    // Stats
    std::vector<double> inf_pr, inf_du, mu, d_norm, regularization_size,
      obj, alpha_pr, alpha_du;
    std::vector<int> ls_trials;
    const char* return_status;
    int iter_count;

    // Meta-data
    std::map<std::string, std::vector<std::string> > var_string_md;
    std::map<std::string, std::vector<int> > var_integer_md;
    std::map<std::string, std::vector<double> > var_numeric_md;
    std::map<std::string, std::vector<std::string> > con_string_md;
    std::map<std::string, std::vector<int> > con_integer_md;
    std::map<std::string, std::vector<double> > con_numeric_md;

    /// Constructor
    IpoptMemory();

    /// Destructor
    ~IpoptMemory();
  };

  /** \brief \pluginbrief{Nlpsol,ipopt}

      @copydoc Nlpsol_doc
      @copydoc plugin_Nlpsol_ipopt
  */
  class CASADI_NLPSOL_IPOPT_EXPORT IpoptInterface : public Nlpsol {
    friend class IpoptUserClass;
  public:
    Sparsity jacg_sp_;
    Sparsity hesslag_sp_;

    explicit IpoptInterface(const std::string& name, const Function& nlp);
    ~IpoptInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "ipopt";}

    // Get name of the class
    std::string class_name() const override { return "IpoptInterface";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new IpoptInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new IpoptMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<IpoptMemory*>(mem);}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// Exact Hessian?
    bool exact_hessian_;

    /// All IPOPT options
    Dict opts_;

    // Ipopt callback functions
    void finalize_solution(IpoptMemory* m, const double* x, const double* z_L, const double* z_U,
                           const double* g, const double* lambda, double obj_value,
                           int iter_count) const;
    bool get_bounds_info(IpoptMemory* m, double* x_l, double* x_u,
                         double* g_l, double* g_u) const;
    bool get_starting_point(IpoptMemory* m, bool init_x, double* x,
                            bool init_z, double* z_L, double* z_U,
                            bool init_lambda, double* lambda) const;
    void get_nlp_info(IpoptMemory* m, int& nx, int& ng,
                      int& nnz_jac_g, int& nnz_h_lag) const;
    int get_number_of_nonlinear_variables() const;
    bool get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const;
    bool intermediate_callback(IpoptMemory* m, const double* x, const double* z_L,
                               const double* z_U, const double* g,
                               const double* lambda, double obj_value, int iter,
                               double inf_pr, double inf_du, double mu, double d_norm,
                               double regularization_size, double alpha_du, double alpha_pr,
                               int ls_trials, bool full_callback) const;
    bool get_var_con_metadata(std::map<std::string, std::vector<std::string> >& var_string_md,
                              std::map<std::string, std::vector<int> >& var_integer_md,
                              std::map<std::string, std::vector<double> >& var_numeric_md,
                              std::map<std::string, std::vector<std::string> >& con_string_md,
                              std::map<std::string, std::vector<int> >& con_integer_md,
                              std::map<std::string, std::vector<double> >& con_numeric_md) const;

    /// A documentation string
    static const std::string meta_doc;

    // Options
    bool pass_nonlinear_variables_;
    std::vector<bool> nl_ex_;
    Dict var_string_md_, var_integer_md_, var_numeric_md_,
      con_string_md_, con_integer_md_, con_numeric_md_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_IPOPT_INTERFACE_HPP
