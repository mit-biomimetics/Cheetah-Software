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


#include "fast_newton.hpp"
#include <iomanip>
#include "../core/linsol_internal.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_ROOTFINDER_FAST_NEWTON_EXPORT
  casadi_register_rootfinder_fast_newton(Rootfinder::Plugin* plugin) {
    plugin->creator = FastNewton::creator;
    plugin->name = "fast_newton";
    plugin->doc = FastNewton::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &FastNewton::options_;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_FAST_NEWTON_EXPORT casadi_load_rootfinder_fast_newton() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_fast_newton);
  }

  FastNewton::FastNewton(const std::string& name, const Function& f)
    : Rootfinder(name, f) {
  }

  FastNewton::~FastNewton() {
    clear_mem();
  }

  Options FastNewton::options_
  = {{&Rootfinder::options_},
     {{"abstol",
       {OT_DOUBLE,
        "Stopping criterion tolerance on ||g||__inf)"}},
      {"abstolStep",
       {OT_DOUBLE,
        "Stopping criterion tolerance on step size"}},
      {"max_iter",
       {OT_INT,
        "Maximum number of Newton iterations to perform before returning."}}
     }
  };

  void FastNewton::init(const Dict& opts) {

    // Call the base class initializer
    Rootfinder::init(opts);

    // Default options
    max_iter_ = 1000;
    abstol_ = 1e-12;
    abstolStep_ = 1e-12;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="abstol") {
        abstol_ = op.second;
      } else if (op.first=="abstolStep") {
        abstolStep_ = op.second;
      }
    }

    casadi_assert(oracle_.n_in()>0,
                          "Newton: the supplied f must have at least one input.");
    casadi_assert(!linsol_.is_null(),
                          "Newton::init: linear_solver must be supplied");

    jac_f_z_ = get_function("jac_f_z");

    // Symbolic factorization
    sp_jac_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Allocate memory
    alloc_w(n_, true); // x
    alloc_w(n_, true); // F
    alloc_w(sp_jac_.nnz(), true); // J

    alloc_w(sp_jac_.size1() + sp_jac_.size2(), true); // w
    alloc_w(sp_v_.nnz(), true); // v
    alloc_w(sp_r_.nnz(), true); // r
    alloc_w(sp_jac_.size2(), true); // beta

  }

 void FastNewton::set_work(void* mem, const double**& arg, double**& res,
                       casadi_int*& iw, double*& w) const {
     Rootfinder::set_work(mem, arg, res, iw, w);
     auto m = static_cast<FastNewtonMemory*>(mem);
     casadi_newton_mem<double>* M = &(m->M);

     M->n = n_;
     M->abstol = abstol_;
     M->abstol_step = abstolStep_;

     M->x = w; w += n_;
     M->g = w; w += n_;
     M->jac_g_x = w; w += sp_jac_.nnz();

     M->sp_a = sp_jac_;
     M->sp_v = sp_v_;
     M->sp_r = sp_r_;
     M->prinv = get_ptr(prinv_);
     M->pc = get_ptr(pc_);

     M->lin_w = w; w+= sp_jac_.size1()+sp_jac_.size2();
     M->lin_v = w; w+= sp_v_.nnz();
     M->lin_r = w; w+= sp_r_.nnz();
     M->lin_beta = w; w+= sp_jac_.size2();

  }

  int FastNewton::solve(void* mem) const {
    auto m = static_cast<FastNewtonMemory*>(mem);
    casadi_newton_mem<double>* M = &(m->M);

    // Get the initial guess
    casadi_copy(m->iarg[iin_], n_, M->x);

    for (m->iter=0; m->iter<max_iter_; ++m->iter) {
       /* (re)calculate f and J */
       // Use x to evaluate J
       for (casadi_int i=0;i<n_in_;++i) m->arg[i] = m->iarg[i];
       m->arg[iin_] = M->x;
       for (casadi_int i=0;i<n_out_;++i) m->res[i+1] = m->ires[i];
       m->res[0] = M->jac_g_x;
       m->res[1+iout_] = M->g;
       jac_f_z_(m->arg, m->res, m->iw, m->w);

       m->return_status = casadi_newton(M);
       if (m->return_status) break;
    }
    // Get the solution
    casadi_copy(M->x, n_, m->ires[iout_]);

    m->success = m->return_status>0;

    return 0;
  }

  void FastNewton::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_NEWTON);

    g.local("m", "struct casadi_newton_mem");

    g << "m.n = " << n_ << ";\n";
    g << "m.abstol = " << abstol_ << ";\n";
    g << "m.abstol_step = " << abstolStep_ << ";\n";

    casadi_int w_offset = 0;
    g << "m.x = w;\n"; w_offset+=n_;
    g << "m.g = w+" + str(w_offset) + ";\n"; w_offset+=n_;
    g << "m.jac_g_x = w+" + str(w_offset) + ";\n"; w_offset+=sp_jac_.nnz();

    g << "m.sp_a = " + g.sparsity(sp_jac_)+ ";\n";
    g << "m.sp_v = " + g.sparsity(sp_v_)+ ";\n";
    g << "m.sp_r = " + g.sparsity(sp_r_)+ ";\n";
    g << "m.prinv = " + g.constant(prinv_)+ ";\n";
    g << "m.pc = " + g.constant(pc_)+ ";\n";

    g << "m.lin_w = w+" + str(w_offset) + ";\n"; w_offset+=sp_jac_.size1()+sp_jac_.size2();
    g << "m.lin_v = w+" + str(w_offset) + ";\n"; w_offset+=sp_v_.nnz();
    g << "m.lin_r = w+" + str(w_offset) + ";\n"; w_offset+=sp_r_.nnz();
    g << "m.lin_beta = w+" + str(w_offset) + ";\n"; w_offset+=sp_jac_.size2();

    std::string jac_f_z = g.add_dependency(get_function("jac_f_z"));

    g.comment("Get the initial guess");
    g << g.copy("arg[" + str(iin_)+ "]", n_, "m.x") << "\n";

    g.local("iter", "casadi_int");
    g << "for (iter=0; iter<" + str(max_iter_) + "; ++iter) {\n";

    g.comment("(re)calculate f and J");
    // Use x to evaluate J
    for (casadi_int i=0;i<n_in_;++i) {
      g << "arg[" + str(i+n_in_) + "] = " << (i==iin_? "m.x" : "arg[" + str(i)+ "]") << ";\n";
    }
    g << "res[" + str(n_out_) + "] = m.jac_g_x;\n";
    for (casadi_int i=0;i<n_out_;++i) {
      g << "res[" + str(i+n_out_+1) + "] = " << (i==iout_? "m.g" : "res[" + str(i)+ "]") << ";\n";
    }
    g << jac_f_z + "(arg+" + str(n_in_) + ", res+" + str(n_out_) + ", iw, w, 0);\n";
    g << "if (casadi_newton(&m)) break;\n";
    g << "}\n";

    // Get the solution
    g.comment("Get the solution");
    g << g.copy("m.x", n_, "res[" + str(iout_)+ "]") << "\n";
  }

  void FastNewton::codegen_declarations(CodeGenerator& g) const {
    g.add_dependency(get_function("jac_f_z"));
  }

  int FastNewton::init_mem(void* mem) const {
    if (Rootfinder::init_mem(mem)) return 1;
    auto m = static_cast<FastNewtonMemory*>(mem);
    m->return_status = 0;
    m->iter = 0;
    return 0;
  }

  std::string return_code(casadi_int status) {
    switch (status) {
      case 0:
        return "max_iteration_reached";
      case 1:
         return "converged_abstol";
       case 2:
         return "converged_abstol_step";
      default:
        return "unknown";
    }
  }

  Dict FastNewton::get_stats(void* mem) const {
    Dict stats = Rootfinder::get_stats(mem);
    auto m = static_cast<FastNewtonMemory*>(mem);
    stats["return_status"] = return_code(m->return_status);
    stats["iter_count"] = m->iter;
    return stats;
  }

} // namespace casadi
