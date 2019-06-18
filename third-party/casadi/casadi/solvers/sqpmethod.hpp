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


#ifndef CASADI_SQPMETHOD_HPP
#define CASADI_SQPMETHOD_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/solvers/casadi_nlpsol_sqpmethod_export.h>

/** \defgroup plugin_Nlpsol_sqpmethod
 A textbook SQPMethod
*/

/** \pluginsection{Nlpsol,sqpmethod} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_SQPMETHOD_EXPORT SqpmethodMemory : public NlpsolMemory {
    /// Current and previous linearization point and candidate
    double *x_cand;

    /// Lagrange gradient in the next iterate
    double *gLag, *gLag_old;

    /// Constraint function value
    double *g_cand;

    /// Gradient of the objective function
    double *gf;

    // Bounds of the QP
    double *qp_LBA, *qp_UBA, *qp_LBX, *qp_UBX;

    // QP solution
    double *dx, *qp_DUAL_X, *qp_DUAL_A;

    // Current Jacobian
    double *Jk;

    /// Current Hessian approximation
    double *Bk;

    /// Hessian regularization
    double reg;

    /// Linesearch parameters
    double sigma;

    // Storage for merit function
    double* merit_mem;
    size_t merit_ind;

    /// Last return status
    const char* return_status;

    /// Iteration count
    int iter_count;
  };

  /** \brief  \pluginbrief{Nlpsol,sqpmethod}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_sqpmethod
  */
  class CASADI_NLPSOL_SQPMETHOD_EXPORT Sqpmethod : public Nlpsol {
  public:
    explicit Sqpmethod(const std::string& name, const Function& nlp);
    ~Sqpmethod() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "sqpmethod";}

    // Name of the class
    std::string class_name() const override { return "Sqpmethod";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Sqpmethod(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SqpmethodMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SqpmethodMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// QP solver for the subproblems
    Function qpsol_;

    /// Exact Hessian?
    bool exact_hessian_;

    /// Maximum, minimum number of SQP iterations
    casadi_int max_iter_, min_iter_;

    /// Memory size of L-BFGS method
    casadi_int lbfgs_memory_;

    /// Tolerance of primal and dual infeasibility
    double tol_pr_, tol_du_;

    /// Minimum step size allowed
    double min_step_size_;

    /// Linesearch parameters
    ///@{
    double c1_;
    double beta_;
    casadi_int max_iter_ls_;
    casadi_int merit_memsize_;
    ///@}

    // Print options
    bool print_header_, print_iteration_, print_status_;

    // Hessian sparsity
    Sparsity Hsp_;

    // Jacobian sparsity
    Sparsity Asp_;

    /// Regularization
    bool regularize_;

    /// Access Conic
    const Function getConic() const { return qpsol_;}

    /// Print iteration header
    void print_iteration() const;

    /// Print iteration
    void print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf,
                         double dx_norm, double reg, casadi_int ls_trials, bool ls_success) const;

    // Solve the QP subproblem
    virtual void solve_QP(SqpmethodMemory* m, const double* H, const double* g,
                          const double* lbx, const double* ubx,
                          const double* A, const double* lbA, const double* ubA,
                          double* x_opt, double* lambda_x_opt, double* lambda_A_opt) const;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SQPMETHOD_HPP
