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


#include "bonmin_nlp.hpp"
#include "bonmin_interface.hpp"
#include "casadi/core/timing.hpp"

namespace casadi {

  BonminUserClass::BonminUserClass(const BonminInterface& solver, BonminMemory* mem)
    : solver_(solver), mem_(mem) {
    n_ = solver_.nx_;
    m_ = solver_.ng_;
  }

  BonminUserClass::~BonminUserClass() {
  }

  // returns the size of the problem
  bool BonminUserClass::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                    Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style) {
    solver_.get_nlp_info(mem_, n, m, nnz_jac_g, nnz_h_lag);

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
  }

  // returns the variable bounds
  bool BonminUserClass::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                       Index m, Number* g_l, Number* g_u) {
    casadi_assert_dev(n==solver_.nx_);
    casadi_assert_dev(m==solver_.ng_);
    return solver_.get_bounds_info(mem_, x_l, x_u, g_l, g_u);
  }

  // returns the initial point for the problem
  bool BonminUserClass::get_starting_point(Index n, bool init_x, Number* x,
                                          bool init_z, Number* z_L, Number* z_U,
                                          Index m, bool init_lambda,
                                          Number* lambda) {
    casadi_assert_dev(n==solver_.nx_);
    casadi_assert_dev(m==solver_.ng_);
    return solver_.get_starting_point(mem_, init_x, x, init_z, z_L, z_U, init_lambda, lambda);
  }

  // returns the value of the objective function
  bool BonminUserClass::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    mem_->arg[0] = x;
    mem_->arg[1] = mem_->p;
    mem_->res[0] = &obj_value;
    return solver_.calc_function(mem_, "nlp_f")==0;
  }

  // return the gradient of the objective function grad_ {x} f(x)
  bool BonminUserClass::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
    mem_->arg[0] = x;
    mem_->arg[1] = mem_->p;
    mem_->res[0] = nullptr;
    mem_->res[1] = grad_f;
    return solver_.calc_function(mem_, "nlp_grad_f")==0;
  }

  // return the value of the constraints: g(x)
  bool BonminUserClass::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
    mem_->arg[0] = x;
    mem_->arg[1] = mem_->p;
    mem_->res[0] = g;
    return solver_.calc_function(mem_, "nlp_g")==0;
  }

  // return the structure or values of the jacobian
  bool BonminUserClass::eval_jac_g(Index n, const Number* x, bool new_x,
                                  Index m, Index nele_jac, Index* iRow, Index *jCol,
                                  Number* values) {
    if (values) {
      // Evaluate numerically
      mem_->arg[0] = x;
      mem_->arg[1] = mem_->p;
      mem_->res[0] = nullptr;
      mem_->res[1] = values;
      return solver_.calc_function(mem_, "nlp_jac_g")==0;
    } else {
      // Get the sparsity pattern
      casadi_int ncol = solver_.jacg_sp_.size2();
      const casadi_int* colind = solver_.jacg_sp_.colind();
      const casadi_int* row = solver_.jacg_sp_.row();
      if (nele_jac!=colind[ncol]) return false; // consistency check

      // Pass to BONMIN
      for (casadi_int cc=0; cc<ncol; ++cc) {
        for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
          *iRow++ = row[el];
          *jCol++ = cc;
        }
      }
      return true;
    }
  }


  bool BonminUserClass::eval_h(Index n, const Number* x, bool new_x,
                              Number obj_factor, Index m, const Number* lambda,
                              bool new_lambda, Index nele_hess, Index* iRow,
                              Index* jCol, Number* values) {
    if (values) {
      // Evaluate numerically
      mem_->arg[0] = x;
      mem_->arg[1] = mem_->p;
      mem_->arg[2] = &obj_factor;
      mem_->arg[3] = lambda;
      mem_->res[0] = values;
      if (solver_.calc_function(mem_, "nlp_hess_l")) return false;
      return true;
    } else {
      // Get the sparsity pattern
      casadi_int ncol = solver_.hesslag_sp_.size2();
      const casadi_int* colind = solver_.hesslag_sp_.colind();
      const casadi_int* row = solver_.hesslag_sp_.row();

      // Pass to BONMIN
      for (casadi_int cc=0; cc<ncol; ++cc) {
        for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
          *iRow++ = row[el];
          *jCol++ = cc;
        }
      }
      return true;
    }
  }

  void BonminUserClass::finalize_solution(TMINLP::SolverReturn status,
        Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value) {
    solver_.finalize_solution(mem_, status, x, obj_value);
  }


  bool BonminUserClass::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value,
                                             Number inf_pr, Number inf_du,
                                             Number mu, Number d_norm,
                                             Number regularization_size,
                                             Number alpha_du, Number alpha_pr,
                                             Index ls_trials,
                                             const IpoptData* ip_data,
                                             IpoptCalculatedQuantities* ip_cq) {

    // Only do the callback every few iterations
    if (iter % solver_.callback_step_!=0) return true;

    /// Code copied from TNLPAdapter::FinalizeSolution
    /// See also: http://list.coin-or.org/pipermail/ipopt/2010-July/002078.html
    // http://list.coin-or.org/pipermail/ipopt/2010-April/001965.html

    bool full_callback = false;

    return solver_.intermediate_callback(mem_, x_, z_L_, z_U_, g_, lambda_, obj_value, iter,
                                         inf_pr, inf_du, mu, d_norm, regularization_size,
                                         alpha_du, alpha_pr, ls_trials, full_callback);
  }

  Index BonminUserClass::get_number_of_nonlinear_variables() {
    return solver_.get_number_of_nonlinear_variables();
  }

  bool BonminUserClass::get_list_of_nonlinear_variables(Index num_nonlin_vars,
                                                       Index* pos_nonlin_vars) {
    return solver_.get_list_of_nonlinear_variables(num_nonlin_vars, pos_nonlin_vars);
  }

 bool BonminUserClass::get_variables_types(Index n, VariableType* var_types) {
   if (solver_.discrete_.empty()) {
     std::fill_n(var_types, n, CONTINUOUS);
   } else {
     if (solver_.discrete_.size()!=n) return false;
     for (auto&& d : solver_.discrete_) {
       *var_types++ = d ? INTEGER : CONTINUOUS;
     }
   }
   return true;
 }

 bool BonminUserClass::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types) {
   casadi_assert_dev(n==solver_.nl_ex_.size());
   for (casadi_int i=0; i<n; ++i)
    var_types[i] = solver_.nl_ex_[i] ? Ipopt::TNLP::NON_LINEAR : Ipopt::TNLP::LINEAR;
   return true;
 }

 bool BonminUserClass::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types) {
   casadi_assert_dev(m==solver_.nl_g_.size());
   for (casadi_int i=0; i<m; ++i)
    const_types[i] = solver_.nl_g_[i] ? Ipopt::TNLP::NON_LINEAR : Ipopt::TNLP::LINEAR;
   return true;
 }

 const Bonmin::TMINLP::SosInfo * BonminUserClass::sosConstraints() const {
   return &solver_.sosConstraints(mem_);
 }



} // namespace casadi
