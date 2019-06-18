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


#include "qrqp.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QRQP_EXPORT
  casadi_register_conic_qrqp(Conic::Plugin* plugin) {
    plugin->creator = Qrqp::creator;
    plugin->name = "qrqp";
    plugin->doc = Qrqp::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Qrqp::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QRQP_EXPORT casadi_load_conic_qrqp() {
    Conic::registerPlugin(casadi_register_conic_qrqp);
  }

  Qrqp::Qrqp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  Qrqp::~Qrqp() {
    clear_mem();
  }

  Options Qrqp::options_
  = {{&Conic::options_},
     {{"max_iter",
       {OT_INT,
        "Maximum number of iterations [1000]."}},
      {"tol",
       {OT_DOUBLE,
        "Tolerance [1e-8]."}},
      {"du_to_pr",
       {OT_DOUBLE,
        "How much larger dual than primal error is acceptable [1000]"}},
      {"print_header",
       {OT_BOOL,
        "Print header [true]."}},
      {"print_iter",
       {OT_BOOL,
        "Print iterations [true]."}}
     }
  };

  void Qrqp::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    max_iter_ = 1000;
    print_iter_ = true;
    print_header_ = true;
    du_to_pr_ = 1000.;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="print_iter") {
        print_iter_ = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="du_to_pr") {
        du_to_pr_ = op.second;
      }
    }

    // Transpose of the Jacobian
    AT_ = A_.T();

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, true, true);

    // Symbolic QR factorization
    kkt_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Setup memory structure
    p_.du_to_pr = du_to_pr_;
    p_.print_iter = print_iter_;
    p_.sp_a = A_;
    p_.sp_h = H_;
    p_.sp_at = AT_;
    p_.sp_kkt = kkt_;
    p_.sp_v = sp_v_;
    p_.sp_r = sp_r_;
    p_.prinv = get_ptr(prinv_);
    p_.pc = get_ptr(pc_);
    p_.dmin = std::numeric_limits<double>::min();
    p_.inf = inf;
    p_.nx = nx_;
    p_.na = na_;
    p_.nz = nx_+na_;

    // Allocate memory
    casadi_int sz_w, sz_iw;
    casadi_qp_work(&p_, &sz_iw, &sz_w);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::QRQP\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Number of nonzeros in H:                   %9d\n", H_.nnz());
      print("Number of nonzeros in A:                   %9d\n", A_.nnz());
      print("Number of nonzeros in KKT:                 %9d\n", kkt_.nnz());
      print("Number of nonzeros in QR(V):               %9d\n", sp_v_.nnz());
      print("Number of nonzeros in QR(R):               %9d\n", sp_r_.nnz());
    }
  }

  int Qrqp::init_mem(void* mem) const {
    auto m = static_cast<QrqpMemory*>(mem);
    m->return_status = "";
    m->success = false;
    return 0;
  }

  int Qrqp::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QrqpMemory*>(mem);
    // Reset statistics
    for (auto&& s : m->fstats) s.second.reset();
    // Check inputs
    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }
    // Setup data structure
    casadi_qp_data<double> d;
    d.prob = &p_;
    d.nz_h = arg[CONIC_H];
    d.g = arg[CONIC_G];
    d.nz_a = arg[CONIC_A];
    casadi_qp_init(&d, iw, w);
    // Pass bounds on z
    casadi_copy(arg[CONIC_LBX], nx_, d.lbz);
    casadi_copy(arg[CONIC_LBA], na_, d.lbz+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, d.ubz);
    casadi_copy(arg[CONIC_UBA], na_, d.ubz+nx_);
    // Pass initial guess
    casadi_copy(arg[CONIC_X0], nx_, d.z);
    casadi_copy(arg[CONIC_LAM_X0], nx_, d.lam);
    casadi_copy(arg[CONIC_LAM_A0], na_, d.lam+nx_);
    // Reset solver
    if (casadi_qp_reset(&d)) return 1;
    // Return flag
    int flag = 0;
    // Constraint to be flipped, if any
    casadi_int index=-2, sign=0, r_index=-2, r_sign=0;
    // QP iterations
    casadi_int iter = 0;
    while (true) {
      // Calculate dependent quantities
      casadi_qp_calc_dependent(&d);
      // Make an active set change
      casadi_qp_flip(&d, &index, &sign, r_index, r_sign);
      // Form and factorize the KKT system
      casadi_qp_factorize(&d);
      // Termination message
      if (index==-1) {
        casadi_qp_log(&d, "QP converged");
        m->return_status = "success";
      } else if (iter>=max_iter_) {
        casadi_qp_log(&d, "QP terminated: max iter");
        m->return_status = "Maximum number of iterations reached";
        flag = 1;
      }
      // Print iteration progress:
      if (print_iter_) {
        if (iter % 10 == 0) {
          print("%5s %5s %9s %9s %5s %9s %5s %9s %5s %9s %40s\n",
                "Iter", "Sing", "fk", "|pr|", "con", "|du|", "var",
                "min_R", "con", "last_tau", "Note");
        }
        print("%5d %5d %9.2g %9.2g %5d %9.2g %5d %9.2g %5d %9.2g %40s\n",
              iter, d.sing, d.f, d.pr, d.ipr, d.du, d.idu,
              d.mina, d.imina, d.tau, d.msg);
        d.msg[0] = '\0';
      }
      // Terminate loop?
      if (index==-1 || flag!=0) break;
      // Start a new iteration
      iter++;
      // Calculate search direction
      if (casadi_qp_calc_step(&d, &r_index, &r_sign)) {
        if (print_iter_) print("QP terminated: No search direction\n");
        m->return_status = "Failed to calculate search direction";
        flag = 1;
        break;
      }
      // Line search in the calculated direction
      casadi_qp_linesearch(&d, &index, &sign);
    }
    // Get solution
    casadi_copy(&d.f, 1, res[CONIC_COST]);
    casadi_copy(d.z, nx_, res[CONIC_X]);
    casadi_copy(d.lam, nx_, res[CONIC_LAM_X]);
    casadi_copy(d.lam+nx_, na_, res[CONIC_LAM_A]);
    // Return
    if (verbose_) casadi_warning(m->return_status);
    m->success = flag ? false : true;
    return 0;
  }

  Dict Qrqp::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QrqpMemory*>(mem);
    stats["return_status"] = m->return_status;
    stats["success"] = m->success;
    return stats;
  }

} // namespace casadi
