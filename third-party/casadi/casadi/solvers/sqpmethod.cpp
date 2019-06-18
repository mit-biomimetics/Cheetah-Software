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


#include "sqpmethod.hpp"

#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/calculus.hpp"
#include "casadi/core/conic.hpp"

#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_SQPMETHOD_EXPORT
      casadi_register_nlpsol_sqpmethod(Nlpsol::Plugin* plugin) {
    plugin->creator = Sqpmethod::creator;
    plugin->name = "sqpmethod";
    plugin->doc = Sqpmethod::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Sqpmethod::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SQPMETHOD_EXPORT casadi_load_nlpsol_sqpmethod() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_sqpmethod);
  }

  Sqpmethod::Sqpmethod(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  Sqpmethod::~Sqpmethod() {
    clear_mem();
  }

  Options Sqpmethod::options_
  = {{&Nlpsol::options_},
     {{"qpsol",
       {OT_STRING,
        "The QP solver to be used by the SQP method [qpoases]"}},
      {"qpsol_options",
       {OT_DICT,
        "Options to be passed to the QP solver"}},
      {"hessian_approximation",
       {OT_STRING,
        "limited-memory|exact"}},
      {"max_iter",
       {OT_INT,
        "Maximum number of SQP iterations"}},
      {"min_iter",
       {OT_INT,
        "Minimum number of SQP iterations"}},
      {"max_iter_ls",
       {OT_INT,
        "Maximum number of linesearch iterations"}},
      {"tol_pr",
       {OT_DOUBLE,
        "Stopping criterion for primal infeasibility"}},
      {"tol_du",
       {OT_DOUBLE,
        "Stopping criterion for dual infeasability"}},
      {"c1",
       {OT_DOUBLE,
        "Armijo condition, coefficient of decrease in merit"}},
      {"beta",
       {OT_DOUBLE,
        "Line-search parameter, restoration factor of stepsize"}},
      {"merit_memory",
       {OT_INT,
        "Size of memory to store history of merit function values"}},
      {"lbfgs_memory",
       {OT_INT,
        "Size of L-BFGS memory."}},
      {"regularize",
       {OT_BOOL,
        "Automatic regularization of Lagrange Hessian."}},
      {"print_header",
       {OT_BOOL,
        "Print the header with problem statistics"}},
      {"print_iteration",
       {OT_BOOL,
        "Print the iterations"}},
      {"print_status",
       {OT_BOOL,
        "Print a status message after solving"}},
      {"min_step_size",
       {OT_DOUBLE,
        "The size (inf-norm) of the step size should not become smaller than this."}},
     }
  };

  void Sqpmethod::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Default options
    min_iter_ = 0;
    max_iter_ = 50;
    max_iter_ls_ = 3;
    c1_ = 1e-4;
    beta_ = 0.8;
    merit_memsize_ = 4;
    lbfgs_memory_ = 10;
    tol_pr_ = 1e-6;
    tol_du_ = 1e-6;
    regularize_ = false;
    string hessian_approximation = "exact";
    min_step_size_ = 1e-10;
    string qpsol_plugin = "qpoases";
    Dict qpsol_options;
    print_header_ = true;
    print_iteration_ = true;
    print_status_ = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="min_iter") {
        min_iter_ = op.second;
      } else if (op.first=="max_iter_ls") {
        max_iter_ls_ = op.second;
      } else if (op.first=="c1") {
        c1_ = op.second;
      } else if (op.first=="beta") {
        beta_ = op.second;
      } else if (op.first=="merit_memory") {
        merit_memsize_ = op.second;
      } else if (op.first=="lbfgs_memory") {
        lbfgs_memory_ = op.second;
      } else if (op.first=="tol_pr") {
        tol_pr_ = op.second;
      } else if (op.first=="tol_du") {
        tol_du_ = op.second;
      } else if (op.first=="hessian_approximation") {
        hessian_approximation = op.second.to_string();
      } else if (op.first=="min_step_size") {
        min_step_size_ = op.second;
      } else if (op.first=="qpsol") {
        qpsol_plugin = op.second.to_string();
      } else if (op.first=="qpsol_options") {
        qpsol_options = op.second;
      } else if (op.first=="regularize") {
        regularize_ = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="print_iteration") {
        print_iteration_ = op.second;
      } else if (op.first=="print_status") {
        print_status_ = op.second;
      }
    }

    // Use exact Hessian?
    exact_hessian_ = hessian_approximation =="exact";

    // Get/generate required functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    // First order derivative information
    Function jac_g_fcn = create_function("nlp_jac_fg", {"x", "p"},
                                        {"f", "grad:f:x", "g", "jac:g:x"});
    Asp_ = jac_g_fcn.sparsity_out(3);

    if (exact_hessian_) {
      Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                           {"sym:hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
      Hsp_ = hess_l_fcn.sparsity_out(0);
    } else {
      Hsp_ = Sparsity::dense(nx_, nx_);
    }


    // Allocate a QP solver
    casadi_assert(!qpsol_plugin.empty(), "'qpsol' option has not been set");
    qpsol_ = conic("qpsol", qpsol_plugin, {{"h", Hsp_}, {"a", Asp_}},
                   qpsol_options);
    alloc(qpsol_);

    // BFGS?
    if (!exact_hessian_) {
      alloc_w(2*nx_); // casadi_bfgs
    }

    // Header
    if (print_header_) {
      print("-------------------------------------------\n");
      print("This is casadi::Sqpmethod.\n");
      if (exact_hessian_) {
        print("Using exact Hessian\n");
      } else {
        print("Using limited memory BFGS Hessian approximation\n");
      }
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", ng_);
      print("Number of nonzeros in constraint Jacobian: %9d\n", Asp_.nnz());
      print("Number of nonzeros in Lagrangian Hessian:  %9d\n", Hsp_.nnz());
      print("\n");
    }

    // Current linearization point
    alloc_w(nx_, true); // x_cand_

    // Lagrange gradient in the next iterate
    alloc_w(nx_, true); // gLag_
    alloc_w(nx_, true); // gLag_old_

    // Constraint function value
    alloc_w(ng_, true); // g_cand_

    // Gradient of the objective
    alloc_w(nx_, true); // gf_

    // Bounds of the QP
    alloc_w(ng_, true); // qp_LBA
    alloc_w(ng_, true); // qp_UBA
    alloc_w(nx_, true); // qp_LBX
    alloc_w(nx_, true); // qp_UBX

    // QP solution
    alloc_w(nx_, true); // dx_
    alloc_w(nx_, true); // qp_DUAL_X_
    alloc_w(ng_, true); // qp_DUAL_A_

    // Hessian approximation
    alloc_w(Hsp_.nnz(), true); // Bk_

    // Jacobian
    alloc_w(Asp_.nnz(), true); // Jk_

    // Line-search memory
    alloc_w(merit_memsize_, true);
  }

  void Sqpmethod::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<SqpmethodMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Current linearization point
    m->x_cand = w; w += nx_;

    // Lagrange gradient in the next iterate
    m->gLag = w; w += nx_;
    m->gLag_old = w; w += nx_;

    // Constraint function value
    m->g_cand = w; w += ng_;

    // Gradient of the objective
    m->gf = w; w += nx_;

    // Bounds of the QP
    m->qp_LBA = w; w += ng_;
    m->qp_UBA = w; w += ng_;
    m->qp_LBX = w; w += nx_;
    m->qp_UBX = w; w += nx_;

    // QP solution
    m->dx = w; w += nx_;
    m->qp_DUAL_X = w; w += nx_;
    m->qp_DUAL_A = w; w += ng_;

    // Hessian approximation
    m->Bk = w; w += Hsp_.nnz();

    // Jacobian
    m->Jk = w; w += Asp_.nnz();

    // merit_mem
    m->merit_mem = w; w += merit_memsize_;

    m->iter_count = -1;
  }

  int Sqpmethod::solve(void* mem) const {
    auto m = static_cast<SqpmethodMemory*>(mem);

    // Number of SQP iterations
    m->iter_count = 0;

    // Number of line-search iterations
    casadi_int ls_iter = 0;

    // Last linesearch successfull
    bool ls_success = true;

    // Reset
    m->merit_ind = 0;
    m->sigma = 0.;    // NOTE: Move this into the main optimization loop
    m->reg = 0;

    // Default stepsize
    double t = 0;

    // For seeds
    const double one = 1.;

    // MAIN OPTIMIZATION LOOP
    while (true) {
      // Evaluate f, g and first order derivative information
      m->arg[0] = m->x;
      m->arg[1] = m->p;
      m->res[0] = &m->f;
      m->res[1] = m->gf;
      m->res[2] = m->g;
      m->res[3] = m->Jk;
      if (calc_function(m, "nlp_jac_fg")) return 1;

      // Evaluate the gradient of the Lagrangian
      casadi_copy(m->gf, nx_, m->gLag);
      casadi_mv(m->Jk, Asp_, m->lam_g, m->gLag, true);
      casadi_axpy(nx_, 1., m->lam_x, m->gLag);

      // Primal infeasability
      double pr_inf = std::fmax(casadi_max_viol(nx_, m->x, m->lbx, m->ubx),
                                casadi_max_viol(ng_, m->g, m->lbg, m->ubg));

      // inf-norm of lagrange gradient
      double gLag_norminf = casadi_norm_inf(nx_, m->gLag);

      // inf-norm of step
      double dx_norminf = casadi_norm_inf(nx_, m->dx);

      // Printing information about the actual iterate
      if (print_iteration_) {
        if (m->iter_count % 10 == 0) print_iteration();
        print_iteration(m->iter_count, m->f, pr_inf, gLag_norminf, dx_norminf,
                        m->reg, ls_iter, ls_success);
      }

      // Callback function
      if (callback(m, m->x, &m->f, m->g, m->lam_x, m->lam_g, nullptr)) {
        if (print_status_) print("WARNING(sqpmethod): Aborted by callback...\n");
        m->return_status = "User_Requested_Stop";
        break;
      }

      // Checking convergence criteria
      if (m->iter_count >= min_iter_ && pr_inf < tol_pr_ && gLag_norminf < tol_du_) {
        if (print_status_)
          print("MESSAGE(sqpmethod): Convergence achieved after %d iterations\n", m->iter_count);
        m->return_status = "Solve_Succeeded";
        m->success = true;
        break;
      }

      if (m->iter_count >= max_iter_) {
        if (print_status_) print("MESSAGE(sqpmethod): Maximum number of iterations reached.\n");
        m->return_status = "Maximum_Iterations_Exceeded";
        break;
      }

      if (m->iter_count >= 1 && m->iter_count >= min_iter_ && dx_norminf <= min_step_size_) {
        if (print_status_) print("MESSAGE(sqpmethod): Search direction becomes too small without "
              "convergence criteria being met.\n");
        m->return_status = "Search_Direction_Becomes_Too_Small";
        break;
      }

      if (exact_hessian_) {
        // Update/reset exact Hessian
        m->arg[0] = m->x;
        m->arg[1] = m->p;
        m->arg[2] = &one;
        m->arg[3] = m->lam_g;
        m->res[0] = m->Bk;
        if (calc_function(m, "nlp_hess_l")) return 1;

        // Determing regularization parameter with Gershgorin theorem
        if (regularize_) {
          m->reg = std::fmin(0, -casadi_lb_eig(Hsp_, m->Bk));
          if (m->reg > 0) casadi_regularize(Hsp_, m->Bk, m->reg);
        }
      } else if (m->iter_count==0) {
        // Initialize BFGS
        casadi_fill(m->Bk, Hsp_.nnz(), 1.);
        casadi_bfgs_reset(Hsp_, m->Bk);
      } else {
        // Update BFGS
        if (m->iter_count % lbfgs_memory_ == 0) casadi_bfgs_reset(Hsp_, m->Bk);
        // Update the Hessian approximation
        casadi_bfgs(Hsp_, m->Bk, m->dx, m->gLag, m->gLag_old, m->w);
      }

      // Formulate the QP
      casadi_copy(m->lbx, nx_, m->qp_LBX);
      casadi_axpy(nx_, -1., m->x, m->qp_LBX);
      casadi_copy(m->ubx, nx_, m->qp_UBX);
      casadi_axpy(nx_, -1., m->x, m->qp_UBX);
      casadi_copy(m->lbg, ng_, m->qp_LBA);
      casadi_axpy(ng_, -1., m->g, m->qp_LBA);
      casadi_copy(m->ubg, ng_, m->qp_UBA);
      casadi_axpy(ng_, -1., m->g, m->qp_UBA);

      // Intitial guess
      casadi_copy(m->lam_x, nx_, m->qp_DUAL_X);
      casadi_copy(m->lam_g, ng_, m->qp_DUAL_A);
      casadi_fill(m->dx, nx_, 0.);

      // Increase counter
      m->iter_count++;

      // Solve the QP
      solve_QP(m, m->Bk, m->gf, m->qp_LBX, m->qp_UBX, m->Jk, m->qp_LBA,
               m->qp_UBA, m->dx, m->qp_DUAL_X, m->qp_DUAL_A);

      // Detecting indefiniteness
      double gain = casadi_bilin(m->Bk, Hsp_, m->dx, m->dx);
      if (gain < 0) {
        if (print_status_) print("WARNING(sqpmethod): Indefinite Hessian detected\n");
      }

      // Calculate penalty parameter of merit function
      m->sigma = std::fmax(m->sigma, 1.01*casadi_norm_inf(nx_, m->qp_DUAL_X));
      m->sigma = std::fmax(m->sigma, 1.01*casadi_norm_inf(ng_, m->qp_DUAL_A));

      // Calculate L1-merit function in the actual iterate
      double l1_infeas = std::fmax(casadi_max_viol(nx_, m->x, m->lbx, m->ubx),
                                   casadi_max_viol(ng_, m->g, m->lbg, m->ubg));

      // Right-hand side of Armijo condition
      double F_sens = casadi_dot(nx_, m->dx, m->gf);
      double L1dir = F_sens - m->sigma * l1_infeas;
      double L1merit = m->f + m->sigma * l1_infeas;

      // Storing the actual merit function value in a list
      m->merit_mem[m->merit_ind] = L1merit;
      ++m->merit_ind %= merit_memsize_;

      // Calculating maximal merit function value so far
      double meritmax = m->merit_mem[0];
      for (size_t i=1; i<merit_memsize_ && i<m->iter_count; ++i) {
        if (meritmax < m->merit_mem[i]) meritmax = m->merit_mem[i];
      }

      // Stepsize
      t = 1.0;
      double fk_cand;
      // Merit function value in candidate
      double L1merit_cand = 0;

      // Reset line-search counter, success marker
      ls_iter = 0;
      ls_success = true;

      // Line-search
      if (verbose_) print("Starting line-search\n");
      if (max_iter_ls_>0) { // max_iter_ls_== 0 disables line-search

        // Line-search loop
        while (true) {
          // Increase counter
          ls_iter++;

          // Candidate step
          casadi_copy(m->x, nx_, m->x_cand);
          casadi_axpy(nx_, t, m->dx, m->x_cand);

          // Evaluating objective and constraints
          m->arg[0] = m->x_cand;
          m->arg[1] = m->p;
          m->res[0] = &fk_cand;
          m->res[1] = m->g_cand;
          if (calc_function(m, "nlp_fg")) {
            // line-search failed, skip iteration
            t = beta_ * t;
            continue;
          }

          // Calculating merit-function in candidate
          l1_infeas = std::fmax(casadi_max_viol(nx_, m->x_cand, m->lbx, m->ubx),
                                casadi_max_viol(ng_, m->g_cand, m->lbg, m->ubg));
          L1merit_cand = fk_cand + m->sigma * l1_infeas;
          if (L1merit_cand <= meritmax + t * c1_ * L1dir) {
            break;
          }

          // Line-search not successful, but we accept it.
          if (ls_iter == max_iter_ls_) {
            ls_success = false;
            break;
          }

          // Backtracking
          t = beta_ * t;
        }

        // Candidate accepted, update dual variables
        casadi_scal(ng_, 1-t, m->lam_g);
        casadi_axpy(ng_, t, m->qp_DUAL_A, m->lam_g);
        casadi_scal(nx_, 1-t, m->lam_x);
        casadi_axpy(nx_, t, m->qp_DUAL_X, m->lam_x);

        casadi_scal(nx_, t, m->dx);

      } else {
        // Full step
        casadi_copy(m->qp_DUAL_A, ng_, m->lam_g);
        casadi_copy(m->qp_DUAL_X, nx_, m->lam_x);
      }

      // Take step
      casadi_axpy(nx_, 1., m->dx, m->x);

      if (!exact_hessian_) {
        // Evaluate the gradient of the Lagrangian with the old x but new lam_g (for BFGS)
        casadi_copy(m->gf, nx_, m->gLag_old);
        casadi_mv(m->Jk, Asp_, m->lam_g, m->gLag_old, true);
        casadi_axpy(nx_, 1., m->lam_x, m->gLag_old);
      }
    }

    return 0;
  }

  void Sqpmethod::print_iteration() const {
    print("%4s %14s %9s %9s %9s %7s %2s\n", "iter", "objective", "inf_pr",
          "inf_du", "||d||", "lg(rg)", "ls");
  }

  void Sqpmethod::print_iteration(casadi_int iter, double obj,
                                  double pr_inf, double du_inf,
                                  double dx_norm, double rg,
                                  casadi_int ls_trials, bool ls_success) const {
    print("%4d %14.6e %9.2e %9.2e %9.2e ", iter, obj, pr_inf, du_inf, dx_norm);
    if (rg>0) {
      print("%7.2d ", log10(rg));
    } else {
      print("%7s ", "-");
    }
    print("%2d", ls_trials);
    if (!ls_success) print("F");
    print("\n");
  }

  void Sqpmethod::solve_QP(SqpmethodMemory* m, const double* H, const double* g,
                           const double* lbx, const double* ubx,
                           const double* A, const double* lbA, const double* ubA,
                           double* x_opt, double* lambda_x_opt, double* lambda_A_opt) const {
    // Inputs
    fill_n(m->arg, qpsol_.n_in(), nullptr);
    m->arg[CONIC_H] = H;
    m->arg[CONIC_G] = g;
    m->arg[CONIC_X0] = x_opt;
    m->arg[CONIC_LAM_X0] = lambda_x_opt;
    m->arg[CONIC_LAM_A0] = lambda_A_opt;
    m->arg[CONIC_LBX] = lbx;
    m->arg[CONIC_UBX] = ubx;
    m->arg[CONIC_A] = A;
    m->arg[CONIC_LBA] = lbA;
    m->arg[CONIC_UBA] = ubA;

    // Outputs
    fill_n(m->res, qpsol_.n_out(), nullptr);
    m->res[CONIC_X] = x_opt;
    m->res[CONIC_LAM_X] = lambda_x_opt;
    m->res[CONIC_LAM_A] = lambda_A_opt;

    // Solve the QP
    qpsol_(m->arg, m->res, m->iw, m->w, 0);
    if (verbose_) print("QP solved\n");
  }

  Dict Sqpmethod::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<SqpmethodMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["iter_count"] = m->iter_count;
    return stats;
  }
} // namespace casadi
