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


#include "sundials_interface.hpp"

#include "casadi/core/casadi_misc.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace casadi {

  SundialsInterface::SundialsInterface(const std::string& name, const Function& dae)
    : Integrator(name, dae) {
  }

  SundialsInterface::~SundialsInterface() {
  }

  Options SundialsInterface::options_
  = {{&Integrator::options_},
     {{"max_num_steps",
       {OT_INT,
        "Maximum number of integrator steps"}},
      {"reltol",
       {OT_DOUBLE,
        "Relative tolerence for the IVP solution"}},
      {"abstol",
       {OT_DOUBLE,
        "Absolute tolerence for the IVP solution"}},
      {"newton_scheme",
       {OT_STRING,
        "Linear solver scheme in the Newton method: DIRECT|gmres|bcgstab|tfqmr"}},
      {"max_krylov",
       {OT_INT,
        "Maximum Krylov subspace size"}},
      {"sensitivity_method",
       {OT_STRING,
        "Sensitivity method: SIMULTANEOUS|staggered"}},
      {"max_multistep_order",
       {OT_INT,
        "Maximum order for the (variable-order) multistep method"}},
      {"use_preconditioner",
       {OT_BOOL,
        "Precondition the iterative solver [default: true]"}},
      {"stop_at_end",
       {OT_BOOL,
        "Stop the integrator at the end of the interval"}},
      {"disable_internal_warnings",
       {OT_BOOL,
        "Disable SUNDIALS internal warning messages"}},
      {"quad_err_con",
       {OT_BOOL,
        "Should the quadratures affect the step size control"}},
      {"fsens_err_con",
       {OT_BOOL,
        "include the forward sensitivities in all error controls"}},
      {"steps_per_checkpoint",
       {OT_INT,
        "Number of steps between two consecutive checkpoints"}},
      {"interpolation_type",
       {OT_STRING,
        "Type of interpolation for the adjoint sensitivities"}},
      {"linear_solver",
       {OT_STRING,
        "A custom linear solver creator function [default: qr]"}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver"}},
      {"second_order_correction",
       {OT_BOOL,
        "Second order correction in the augmented system Jacobian [true]"}},
      {"step0",
       {OT_DOUBLE,
        "initial step size [default: 0/estimated]"}},
      {"max_order",
       {OT_DOUBLE,
        "Maximum order"}},
      {"nonlin_conv_coeff",
       {OT_DOUBLE,
        "Coefficient in the nonlinear convergence test"}}
     }
  };

  void SundialsInterface::init(const Dict& opts) {
    // Call the base class method
    Integrator::init(opts);

    // If sensitivity equations, make sure derivative_of_ is available
    casadi_assert(ns_==0 || !derivative_of_.is_null(),
      "Not implemented.");

    // Default options
    abstol_ = 1e-8;
    reltol_ = 1e-6;
    max_num_steps_ = 10000;
    stop_at_end_ = true;
    use_precon_ = true;
    max_krylov_ = 10;
    linear_solver_ = "qr";
    string newton_scheme = "direct";
    quad_err_con_ = false;
    string interpolation_type = "hermite";
    steps_per_checkpoint_ = 20;
    disable_internal_warnings_ = false;
    max_multistep_order_ = 5;
    second_order_correction_ = true;
    step0_ = 0;
    max_order_ = 0;
    nonlin_conv_coeff_ = 0;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="abstol") {
        abstol_ = op.second;
      } else if (op.first=="reltol") {
        reltol_ = op.second;
      } else if (op.first=="max_num_steps") {
        max_num_steps_ = op.second;
      } else if (op.first=="stop_at_end") {
        stop_at_end_ = op.second;
      } else if (op.first=="use_preconditioner") {
        use_precon_ = op.second;
      } else if (op.first=="max_krylov") {
        max_krylov_ = op.second;
      } else if (op.first=="newton_scheme") {
        newton_scheme = op.second.to_string();
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      } else if (op.first=="linear_solver_options") {
        linear_solver_options_ = op.second;
      } else if (op.first=="quad_err_con") {
        quad_err_con_ = op.second;
      } else if (op.first=="interpolation_type") {
        interpolation_type = op.second.to_string();
      } else if (op.first=="steps_per_checkpoint") {
        steps_per_checkpoint_ = op.second;
      } else if (op.first=="disable_internal_warnings") {
        disable_internal_warnings_ = op.second;
      } else if (op.first=="max_multistep_order") {
        max_multistep_order_ = op.second;
      } else if (op.first=="second_order_correction") {
        second_order_correction_ = op.second;
      } else if (op.first=="step0") {
        step0_ = op.second;
      } else if (op.first=="max_order") {
        max_order_ = op.second;
      } else if (op.first=="nonlin_conv_coeff") {
        nonlin_conv_coeff_ = op.second;
      }
    }

    // Type of Newton scheme
    if (newton_scheme=="direct") {
      newton_scheme_ = SD_DIRECT;
    } else if (newton_scheme=="gmres") {
      newton_scheme_ = SD_GMRES;
    } else if (newton_scheme=="bcgstab") {
      newton_scheme_ = SD_BCGSTAB;
    } else if (newton_scheme=="tfqmr") {
      newton_scheme_ = SD_TFQMR;
    } else {
      casadi_error("Unknown Newton scheme: " + newton_scheme);
    }

    // Interpolation_type
    if (interpolation_type=="hermite") {
      interp_ = SD_HERMITE;
    } else if (interpolation_type=="polynomial") {
      interp_ = SD_POLYNOMIAL;
    } else {
      casadi_error("Unknown interpolation type: " + interpolation_type);
    }

    // Get or create Jacobians and linear system solvers
    for (bool backward : {false, true}) {
      // Skip backward?
      if (backward && nrx_==0) continue;

      // Get Jacobian function
      Function J;
      if (ns_==0) {
        J = getJ(backward);
      } else {
        SundialsInterface* d = derivative_of_.get<SundialsInterface>();
        casadi_assert_dev(d!=nullptr);
        if (d->ns_==0) {
          J = d->get_function(backward ? "jacB" : "jacF");
        } else {
          J = d->getJ(backward);
        }
      }
      set_function(J, J.name(), true);
      alloc_w(J.nnz_out(0), true);
    }

    // Allocate work vectors
    alloc_w(np_, true); // p
    alloc_w(nrp_, true); // rp
    alloc_w(2*max(nx_+nz_, nrx_+nrz_), true); // v1, v2

    // Allocate linear solvers
    linsolF_ = Linsol("linsolF", linear_solver_,
      get_function("jacF").sparsity_out(0), linear_solver_options_);
    if (nrx_>0) {
      linsolB_ = Linsol("linsolB", linear_solver_,
        get_function("jacB").sparsity_out(0), linear_solver_options_);
    }
  }

  int SundialsInterface::init_mem(void* mem) const {
    if (Integrator::init_mem(mem)) return 1;
    auto m = static_cast<SundialsMemory*>(mem);

    // Allocate n-vectors
    m->xz = N_VNew_Serial(nx_+nz_);
    m->q = N_VNew_Serial(nq_);
    m->rxz = N_VNew_Serial(nrx_+nrz_);
    m->rq = N_VNew_Serial(nrq_);

    m->mem_linsolF = linsolF_.checkout();
    if (!linsolB_.is_null()) m->mem_linsolB = linsolB_.checkout();

    return 0;
  }

  void SundialsInterface::free_mem(void *mem) const {
    Integrator::free_mem(mem);
    auto m = static_cast<SundialsMemory*>(mem);

    linsolF_.release(m->mem_linsolF);
    if (!linsolB_.is_null()) linsolB_.release(m->mem_linsolB);
  }

  void SundialsInterface::reset(IntegratorMemory* mem, double t, const double* x,
                                const double* z, const double* p) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(p, np_, m->p);

    // Set the state
    casadi_copy(x, nx_, NV_DATA_S(m->xz));
    casadi_copy(z, nz_, NV_DATA_S(m->xz)+nx_);

    // Reset summation states
    N_VConst(0., m->q);
  }

  void SundialsInterface::resetB(IntegratorMemory* mem, double t, const double* rx,
                                 const double* rz, const double* rp) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(rp, nrp_, m->rp);

    // Get the backward state
    casadi_copy(rx, nrx_, NV_DATA_S(m->rxz));

    // Reset summation states
    N_VConst(0., m->rq);
  }

  SundialsMemory::SundialsMemory() {
    this->xz  = nullptr;
    this->q = nullptr;
    this->rxz = nullptr;
    this->rq = nullptr;
    this->first_callB = true;
  }

  SundialsMemory::~SundialsMemory() {
    if (this->xz) N_VDestroy_Serial(this->xz);
    if (this->q) N_VDestroy_Serial(this->q);
    if (this->rxz) N_VDestroy_Serial(this->rxz);
    if (this->rq) N_VDestroy_Serial(this->rq);
  }

  Dict SundialsInterface::get_stats(void* mem) const {
    Dict stats = Integrator::get_stats(mem);
    auto m = static_cast<SundialsMemory*>(mem);

    // Counters, forward problem
    stats["nsteps"] = static_cast<casadi_int>(m->nsteps);
    stats["nfevals"] = static_cast<casadi_int>(m->nfevals);
    stats["nlinsetups"] = static_cast<casadi_int>(m->nlinsetups);
    stats["netfails"] = static_cast<casadi_int>(m->netfails);
    stats["qlast"] = m->qlast;
    stats["qcur"] = m->qcur;
    stats["hinused"] = m->hinused;
    stats["hlast"] = m->hlast;
    stats["hcur"] = m->hcur;
    stats["tcur"] = m->tcur;
    stats["nniters"] = static_cast<casadi_int>(m->nniters);
    stats["nncfails"] = static_cast<casadi_int>(m->nncfails);

    // Counters, backward problem
    stats["nstepsB"] = static_cast<casadi_int>(m->nstepsB);
    stats["nfevalsB"] = static_cast<casadi_int>(m->nfevalsB);
    stats["nlinsetupsB"] = static_cast<casadi_int>(m->nlinsetupsB);
    stats["netfailsB"] = static_cast<casadi_int>(m->netfailsB);
    stats["qlastB"] = m->qlastB;
    stats["qcurB"] = m->qcurB;
    stats["hinusedB"] = m->hinusedB;
    stats["hlastB"] = m->hlastB;
    stats["hcurB"] = m->hcurB;
    stats["tcurB"] = m->tcurB;
    stats["nnitersB"] = static_cast<casadi_int>(m->nnitersB);
    stats["nncfailsB"] = static_cast<casadi_int>(m->nncfailsB);
    return stats;
  }

  void SundialsInterface::print_stats(IntegratorMemory* mem) const {
    auto m = to_mem(mem);
    print("FORWARD INTEGRATION:\n");
    print("Number of steps taken by SUNDIALS: %ld\n", m->nsteps);
    print("Number of calls to the user’s f function: %ld\n", m->nfevals);
    print("Number of calls made to the linear solver setup function: %ld\n", m->nlinsetups);
    print("Number of error test failures: %ld\n", m->netfails);
    print("Method order used on the last internal step: %d\n", m->qlast);
    print("Method order to be used on the next internal step: %d\n", m->qcur);
    print("Actual value of initial step size: %g\n", m->hinused);
    print("Step size taken on the last internal step: %g\n", m->hlast);
    print("Step size to be attempted on the next internal step: %g\n", m->hcur);
    print("Current internal time reached: %g\n");
    print("Number of nonlinear iterations performed: %ld\n", m->nniters);
    print("Number of nonlinear convergence failures: %ld\n", m->nncfails);
    if (nrx_>0) {
      print("BACKWARD INTEGRATION:\n");
      print("Number of steps taken by SUNDIALS: %ld\n", m->nstepsB);
      print("Number of calls to the user’s f function: %ld\n", m->nfevalsB);
      print("Number of calls made to the linear solver setup function: %ld\n", m->nlinsetupsB);
      print("Number of error test failures: %ld\n", m->netfailsB);
      print("Method order used on the last internal step: %d\n" , m->qlastB);
      print("Method order to be used on the next internal step: %d\n", m->qcurB);
      print("Actual value of initial step size: %g\n", m->hinusedB);
      print("Step size taken on the last internal step: %g\n", m->hlastB);
      print("Step size to be attempted on the next internal step: %g\n", m->hcurB);
      print("Current internal time reached: %g\n", m->tcurB);
      print("Number of nonlinear iterations performed: %ld\n", m->nnitersB);
      print("Number of nonlinear convergence failures: %ld\n", m->nncfailsB);
    }
    print("\n");
  }

  void SundialsInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<SundialsMemory*>(mem);

    // Set work in base classes
    Integrator::set_work(mem, arg, res, iw, w);

    // Work vectors
    m->p = w; w += np_;
    m->rp = w; w += nrp_;
    m->v1 = w; w += max(nx_+nz_, nrx_+nrz_);
    m->v2 = w; w += max(nx_+nz_, nrx_+nrz_);
    m->jac = w; w += get_function("jacF").nnz_out(0);
    if (nrx_>0) {
      m->jacB = w; w += get_function("jacB").nnz_out(0);
    }
  }

} // namespace casadi
