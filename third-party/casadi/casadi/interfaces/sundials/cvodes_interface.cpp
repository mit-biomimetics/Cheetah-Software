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


#include "cvodes_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

#define THROWING(fcn, ...) \
cvodes_error(CASADI_STR(fcn), fcn(__VA_ARGS__))

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_CVODES_EXPORT
  casadi_register_integrator_cvodes(Integrator::Plugin* plugin) {
    plugin->creator = CvodesInterface::creator;
    plugin->name = "cvodes";
    plugin->doc = CvodesInterface::meta_doc.c_str();;
    plugin->version = CASADI_VERSION;
    plugin->options = &CvodesInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_CVODES_EXPORT casadi_load_integrator_cvodes() {
    Integrator::registerPlugin(casadi_register_integrator_cvodes);
  }

  CvodesInterface::CvodesInterface(const std::string& name, const Function& dae)
    : SundialsInterface(name, dae) {
  }

  CvodesInterface::~CvodesInterface() {
    clear_mem();
  }

  Options CvodesInterface::options_
  = {{&SundialsInterface::options_},
     {{"linear_multistep_method",
       {OT_STRING,
        "Integrator scheme: BDF|adams"}},
      {"nonlinear_solver_iteration",
       {OT_STRING,
        "Nonlinear solver type: NEWTON|functional"}},
      {"fsens_all_at_once",
       {OT_BOOL,
        "Calculate all right hand sides of the sensitivity equations at once"}}
     }
  };

  void CvodesInterface::init(const Dict& opts) {
    if (verbose_) casadi_message(name_ + "::init");

    // Initialize the base classes
    SundialsInterface::init(opts);

    // Default options
    string linear_multistep_method = "bdf";
    string nonlinear_solver_iteration = "newton";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="linear_multistep_method") {
        linear_multistep_method = op.second.to_string();
      } else if (op.first=="nonlinear_solver_iteration") {
        nonlinear_solver_iteration = op.second.to_string();
      }
    }

    // Create function
    create_function("odeF", {"x", "p", "t"}, {"ode"});
    create_function("quadF", {"x", "p", "t"}, {"quad"});
    create_function("odeB", {"rx", "rp", "x", "p", "t"}, {"rode"});
    create_function("quadB", {"rx", "rp", "x", "p", "t"}, {"rquad"});

    // Algebraic variables not supported
    casadi_assert(nz_==0 && nrz_==0,
      "CVODES does not support algebraic variables");

    if (linear_multistep_method=="adams") {
      lmm_ = CV_ADAMS;
    } else if (linear_multistep_method=="bdf") {
      lmm_ = CV_BDF;
    } else {
      casadi_error("Unknown linear multistep method: " + linear_multistep_method);
    }

    if (nonlinear_solver_iteration=="newton") {
      iter_ = CV_NEWTON;
    } else if (nonlinear_solver_iteration=="functional") {
      iter_ = CV_FUNCTIONAL;
    } else {
      casadi_error("Unknown nonlinear solver iteration: " + nonlinear_solver_iteration);
    }

    // Attach functions for jacobian information
    if (newton_scheme_!=SD_DIRECT || (ns_>0 && second_order_correction_)) {
      create_function("jtimesF", {"t", "x", "p", "fwd:x"}, {"fwd:ode"});
      if (nrx_>0) {
        create_function("jtimesB",
                        {"t", "x", "p", "rx", "rp", "fwd:rx"}, {"fwd:rode"});
      }
    }
  }

  int CvodesInterface::init_mem(void* mem) const {
    if (SundialsInterface::init_mem(mem)) return 1;
    auto m = to_mem(mem);

    // Create CVodes memory block
    m->mem = CVodeCreate(lmm_, iter_);
    casadi_assert(m->mem!=nullptr, "CVodeCreate: Creation failed");

    // Set error handler function
    THROWING(CVodeSetErrHandlerFn, m->mem, ehfun, m);

    // Set user data
    THROWING(CVodeSetUserData, m->mem, m);

    // Initialize CVodes
    double t0 = 0;
    THROWING(CVodeInit, m->mem, rhs, t0, m->xz);

    // Set tolerances
    THROWING(CVodeSStolerances, m->mem, reltol_, abstol_);

    // Maximum number of steps
    THROWING(CVodeSetMaxNumSteps, m->mem, max_num_steps_);

    // Initial step size
    if (step0_) THROWING(CVodeSetInitStep, m->mem, step0_);

    // Maximum order of method
    if (max_order_) THROWING(CVodeSetMaxOrd, m->mem, max_order_);

    // Coeff. in the nonlinear convergence test
    if (nonlin_conv_coeff_) THROWING(CVodeSetNonlinConvCoef, m->mem, nonlin_conv_coeff_);

    // attach a linear solver
    if (newton_scheme_==SD_DIRECT) {
      // Direct scheme
      CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
      cv_mem->cv_lmem   = m;
      cv_mem->cv_lsetup = lsetup;
      cv_mem->cv_lsolve = lsolve;
      cv_mem->cv_setupNonNull = TRUE;
    } else {
      // Iterative scheme
      casadi_int pretype = use_precon_ ? PREC_LEFT : PREC_NONE;
      switch (newton_scheme_) {
      case SD_DIRECT: casadi_assert_dev(0);
      case SD_GMRES: THROWING(CVSpgmr, m->mem, pretype, max_krylov_); break;
      case SD_BCGSTAB: THROWING(CVSpbcg, m->mem, pretype, max_krylov_); break;
      case SD_TFQMR: THROWING(CVSptfqmr, m->mem, pretype, max_krylov_); break;
      }
      THROWING(CVSpilsSetJacTimesVecFn, m->mem, jtimes);
      if (use_precon_) THROWING(CVSpilsSetPreconditioner, m->mem, psetup, psolve);
    }

    // Quadrature equations
    if (nq_>0) {
      // Initialize quadratures in CVodes
      THROWING(CVodeQuadInit, m->mem, rhsQ, m->q);

      // Should the quadrature errors be used for step size control?
      if (quad_err_con_) {
        THROWING(CVodeSetQuadErrCon, m->mem, true);

        // Quadrature error tolerances
        // TODO(Joel): vector absolute tolerances
        THROWING(CVodeQuadSStolerances, m->mem, reltol_, abstol_);
      }
    }

    // Initialize adjoint sensitivities
    if (nrx_>0) {
      casadi_int interpType = interp_==SD_HERMITE ? CV_HERMITE : CV_POLYNOMIAL;
      THROWING(CVodeAdjInit, m->mem, steps_per_checkpoint_, interpType);
    }

    m->first_callB = true;
    return 0;
  }

  int CvodesInterface::rhs(double t, N_Vector x, N_Vector xdot, void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = NV_DATA_S(x);
      m->arg[1] = m->p;
      m->arg[2] = &t;
      m->res[0] = NV_DATA_S(xdot);
      s.calc_function(m, "odeF");
      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "rhs failed: " << e.what() << endl;
      return -1;
    }
  }

  void CvodesInterface::reset(IntegratorMemory* mem, double t, const double* x,
                              const double* z, const double* _p) const {
    if (verbose_) casadi_message(name_ + "::reset");
    auto m = to_mem(mem);

    // Reset the base classes
    SundialsInterface::reset(mem, t, x, z, _p);

    // Re-initialize
    THROWING(CVodeReInit, m->mem, t, m->xz);

    // Re-initialize quadratures
    if (nq_>0) {
      N_VConst(0.0, m->q);
      THROWING(CVodeQuadReInit, m->mem, m->q);
    }

    // Re-initialize backward integration
    if (nrx_>0) {
      THROWING(CVodeAdjReInit, m->mem);
    }

    // Set the stop time of the integration -- don't integrate past this point
    if (stop_at_end_) setStopTime(m, grid_.back());
  }

  void CvodesInterface::advance(IntegratorMemory* mem, double t, double* x,
                                double* z, double* q) const {
    auto m = to_mem(mem);

    casadi_assert(t>=grid_.front(),
      "CvodesInterface::integrate(" + str(t) + "): "
      "Cannot integrate to a time earlier than t0 (" + str(grid_.front()) + ")");
    casadi_assert(t<=grid_.back() || !stop_at_end_,
      "CvodesInterface::integrate(" + str(t) + "): "
      "Cannot integrate past a time later than tf (" + str(grid_.back()) + ") "
      "unless stop_at_end is set to False.");

    // Integrate, unless already at desired time
    const double ttol = 1e-9;
    if (fabs(m->t-t)>=ttol) {
      // Integrate forward ...
      if (nrx_>0) {
        // ... with taping
        THROWING(CVodeF, m->mem, t, m->xz, &m->t, CV_NORMAL, &m->ncheck);
      } else {
        // ... without taping
        THROWING(CVode, m->mem, t, m->xz, &m->t, CV_NORMAL);
      }

      // Get quadratures
      if (nq_>0) {
        double tret;
        THROWING(CVodeGetQuad, m->mem, &tret, m->q);
      }
    }

    // Set function outputs
    casadi_copy(NV_DATA_S(m->xz), nx_, x);
    casadi_copy(NV_DATA_S(m->q), nq_, q);

    // Get stats
    THROWING(CVodeGetIntegratorStats, m->mem, &m->nsteps, &m->nfevals, &m->nlinsetups,
             &m->netfails, &m->qlast, &m->qcur, &m->hinused,
             &m->hlast, &m->hcur, &m->tcur);
    THROWING(CVodeGetNonlinSolvStats, m->mem, &m->nniters, &m->nncfails);
  }

  void CvodesInterface::resetB(IntegratorMemory* mem, double t, const double* rx,
                               const double* rz, const double* rp) const {
    auto m = to_mem(mem);

    // Reset the base classes
    SundialsInterface::resetB(mem, t, rx, rz, rp);

    if (m->first_callB) {
      // Create backward problem
      THROWING(CVodeCreateB, m->mem, lmm_, iter_, &m->whichB);
      THROWING(CVodeInitB, m->mem, m->whichB, rhsB, grid_.back(), m->rxz);
      THROWING(CVodeSStolerancesB, m->mem, m->whichB, reltol_, abstol_);
      THROWING(CVodeSetUserDataB, m->mem, m->whichB, m);
      if (newton_scheme_==SD_DIRECT) {
        // Direct scheme
        CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
        CVadjMem ca_mem = cv_mem->cv_adj_mem;
        CVodeBMem cvB_mem = ca_mem->cvB_mem;
        cvB_mem->cv_lmem   = m;
        cvB_mem->cv_mem->cv_lmem = m;
        cvB_mem->cv_mem->cv_lsetup = lsetupB;
        cvB_mem->cv_mem->cv_lsolve = lsolveB;
        cvB_mem->cv_mem->cv_setupNonNull = TRUE;
      } else {
        // Iterative scheme
        casadi_int pretype = use_precon_ ? PREC_LEFT : PREC_NONE;
        switch (newton_scheme_) {
        case SD_DIRECT: casadi_assert_dev(0);
        case SD_GMRES: THROWING(CVSpgmrB, m->mem, m->whichB, pretype, max_krylov_); break;
        case SD_BCGSTAB: THROWING(CVSpbcgB, m->mem, m->whichB, pretype, max_krylov_); break;
        case SD_TFQMR: THROWING(CVSptfqmrB, m->mem, m->whichB, pretype, max_krylov_); break;
        }
        THROWING(CVSpilsSetJacTimesVecFnB, m->mem, m->whichB, jtimesB);
        if (use_precon_) THROWING(CVSpilsSetPreconditionerB, m->mem, m->whichB, psetupB, psolveB);
      }

      // Quadratures for the backward problem
      THROWING(CVodeQuadInitB, m->mem, m->whichB, rhsQB, m->rq);
      if (quad_err_con_) {
        THROWING(CVodeSetQuadErrConB, m->mem, m->whichB, true);
        THROWING(CVodeQuadSStolerancesB, m->mem, m->whichB, reltol_, abstol_);
      }

      // Mark initialized
      m->first_callB = false;
    } else {
      THROWING(CVodeReInitB, m->mem, m->whichB, grid_.back(), m->rxz);
      THROWING(CVodeQuadReInitB, m->mem, m->whichB, m->rq);
    }
  }

  void CvodesInterface::retreat(IntegratorMemory* mem, double t,
                                double* rx, double* rz, double* rq) const {
    auto m = to_mem(mem);
    // Integrate, unless already at desired time
    if (t<m->t) {
      THROWING(CVodeB, m->mem, t, CV_NORMAL);
      THROWING(CVodeGetB, m->mem, m->whichB, &m->t, m->rxz);
      if (nrq_>0) {
        THROWING(CVodeGetQuadB, m->mem, m->whichB, &m->t, m->rq);
      }
    }

    // Save outputs
    casadi_copy(NV_DATA_S(m->rxz), nrx_, rx);
    casadi_copy(NV_DATA_S(m->rq), nrq_, rq);

    // Get stats
    CVodeMem cv_mem = static_cast<CVodeMem>(m->mem);
    CVadjMem ca_mem = cv_mem->cv_adj_mem;
    CVodeBMem cvB_mem = ca_mem->cvB_mem;
    THROWING(CVodeGetIntegratorStats, cvB_mem->cv_mem, &m->nstepsB,
           &m->nfevalsB, &m->nlinsetupsB, &m->netfailsB, &m->qlastB,
           &m->qcurB, &m->hinusedB, &m->hlastB, &m->hcurB, &m->tcurB);
    THROWING(CVodeGetNonlinSolvStats, cvB_mem->cv_mem, &m->nnitersB, &m->nncfailsB);
  }

  void CvodesInterface::cvodes_error(const char* module, int flag) {
    // Successfull return or warning
    if (flag>=CV_SUCCESS) return;
    // Construct error message
    char* flagname = CVodeGetReturnFlagName(flag);
    stringstream ss;
    ss << module << " returned \"" << flagname << "\". Consult CVODES documentation.";
    free(flagname);
    casadi_error(ss.str());
  }

  void CvodesInterface::ehfun(int error_code, const char *module, const char *function,
                              char *msg, void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      if (!s.disable_internal_warnings_) {
        uerr() << msg << endl;
      }
    } catch(exception& e) {
      uerr() << "ehfun failed: " << e.what() << endl;
    }
  }

  int CvodesInterface::rhsQ(double t, N_Vector x, N_Vector qdot, void *user_data) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = NV_DATA_S(x);
      m->arg[1] = m->p;
      m->arg[2] = &t;
      m->res[0] = NV_DATA_S(qdot);
      s.calc_function(m, "quadF");
      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "rhsQ failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::rhsB(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                            void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = NV_DATA_S(rx);
      m->arg[1] = m->rp;
      m->arg[2] = NV_DATA_S(x);
      m->arg[3] = m->p;
      m->arg[4] = &t;
      m->res[0] = NV_DATA_S(rxdot);
      s.calc_function(m, "odeB");

      // Negate (note definition of g)
      casadi_scal(s.nrx_, -1., NV_DATA_S(rxdot));

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "rhsB failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::rhsQB(double t, N_Vector x, N_Vector rx,
                             N_Vector rqdot, void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = NV_DATA_S(rx);
      m->arg[1] = m->rp;
      m->arg[2] = NV_DATA_S(x);
      m->arg[3] = m->p;
      m->arg[4] = &t;
      m->res[0] = NV_DATA_S(rqdot);
      s.calc_function(m, "quadB");

      // Negate (note definition of g)
      casadi_scal(s.nrq_, -1., NV_DATA_S(rqdot));

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "rhsQB failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x,
                              N_Vector xdot, void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = &t;
      m->arg[1] = NV_DATA_S(x);
      m->arg[2] = m->p;
      m->arg[3] = NV_DATA_S(v);
      m->res[0] = NV_DATA_S(Jv);
      s.calc_function(m, "jtimesF");
      return 0;
    } catch(casadi_int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "jtimes failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::jtimesB(N_Vector v, N_Vector Jv, double t, N_Vector x,
                               N_Vector rx, N_Vector rxdot, void *user_data ,
                               N_Vector tmpB) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      m->arg[0] = &t;
      m->arg[1] = NV_DATA_S(x);
      m->arg[2] = m->p;
      m->arg[3] = NV_DATA_S(rx);
      m->arg[4] = m->rp;
      m->arg[5] = NV_DATA_S(v);
      m->res[0] = NV_DATA_S(Jv);
      s.calc_function(m, "jtimesB");
      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "jtimes failed: " << e.what() << endl;
      return -1;
    }
  }

  void CvodesInterface::setStopTime(IntegratorMemory* mem, double tf) const {
    // Set the stop time of the integration -- don't integrate past this point
    auto m = to_mem(mem);
    THROWING(CVodeSetStopTime, m->mem, tf);
  }

  int CvodesInterface::psolve(double t, N_Vector x, N_Vector xdot, N_Vector r,
                              N_Vector z, double gamma, double delta, int lr,
                              void *user_data, N_Vector tmp) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;

      // Get right-hand sides in m->v1
      double* v = NV_DATA_S(r);
      casadi_copy(v, s.nx_, m->v1);

      // Solve for undifferentiated right-hand-side, save to output
      if (s.linsolF_.solve(m->jac, m->v1, 1, false, m->mem_linsolF))
        casadi_error("Linear system solve failed");
      v = NV_DATA_S(z); // possibly different from r
      casadi_copy(m->v1, s.nx1_, v);

      // Sensitivity equations
      if (s.ns_>0) {
        // Second order correction
        if (s.second_order_correction_) {
          // The outputs will double as seeds for jtimesF
          casadi_fill(v + s.nx1_, s.nx_ - s.nx1_, 0.);
          m->arg[0] = &t; // t
          m->arg[1] = NV_DATA_S(x); // x
          m->arg[2] = m->p; // p
          m->arg[3] = v; // fwd:x
          m->res[0] = m->v2; // fwd:ode
          s.calc_function(m, "jtimesF");

          // Subtract m->v2 from m->v1, scaled with -gamma
          casadi_axpy(s.nx_ - s.nx1_, m->gamma, m->v2 + s.nx1_, m->v1 + s.nx1_);
        }

        // Solve for sensitivity right-hand-sides
        if (s.linsolF_.solve(m->jac, m->v1 + s.nx1_, s.ns_, false, m->mem_linsolF))
          casadi_error("Linear solve failed");

        // Save to output, reordered
        casadi_copy(m->v1 + s.nx1_, s.nx_-s.nx1_, v+s.nx1_);
      }

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "psolve failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                               N_Vector rvecB, N_Vector zvecB, double gammaB,
                               double deltaB, int lr, void *user_data, N_Vector tmpB) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;

      // Get right-hand sides in m->v1
      double* v = NV_DATA_S(rvecB);
      casadi_copy(v, s.nrx_, m->v1);

      // Solve for undifferentiated right-hand-side, save to output
      if (s.linsolB_.solve(m->jacB, m->v1, 1, false, m->mem_linsolB))
        casadi_error("Linear solve failed");
      v = NV_DATA_S(zvecB); // possibly different from rvecB
      casadi_copy(m->v1, s.nrx1_, v);

      // Sensitivity equations
      if (s.ns_>0) {
        // Second order correction
        if (s.second_order_correction_) {
          // The outputs will double as seeds for jtimesF
          casadi_fill(v + s.nrx1_, s.nrx_ - s.nrx1_, 0.);
          m->arg[0] = &t; // t
          m->arg[1] = NV_DATA_S(x); // x
          m->arg[2] = m->p; // p
          m->arg[3] = NV_DATA_S(xB); // rx
          m->arg[4] = m->rp; // rp
          m->arg[5] = v; // fwd:rx
          m->res[0] = m->v2; // fwd:rode
          s.calc_function(m, "jtimesB");

          // Subtract m->v2 from m->v1, scaled with gammaB
          casadi_axpy(s.nrx_-s.nrx1_, -m->gammaB, m->v2 + s.nrx1_, m->v1 + s.nrx1_);
        }

        // Solve for sensitivity right-hand-sides
        if (s.linsolB_.solve(m->jacB, m->v1 + s.nx1_, s.ns_, false, m->mem_linsolB)) {
          casadi_error("Linear solve failed");
        }

        // Save to output, reordered
        casadi_copy(m->v1 + s.nx1_, s.nx_-s.nx1_, v+s.nx1_);
      }

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "psolveB failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::psetup(double t, N_Vector x, N_Vector xdot, booleantype jok,
                              booleantype *jcurPtr, double gamma, void *user_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Store gamma for later
      m->gamma = gamma;

      // Calculate Jacobian
      double d1 = -gamma, d2 = 1.;
      m->arg[0] = &t;
      m->arg[1] = NV_DATA_S(x);
      m->arg[2] = m->p;
      m->arg[3] = &d1;
      m->arg[4] = &d2;
      m->res[0] = m->jac;
      if (s.calc_function(m, "jacF")) casadi_error("'jacF' calculation failed");

      // Prepare the solution of the linear system (e.g. factorize)
      if (s.linsolF_.nfact(m->jac, m->mem_linsolF)) casadi_error("'jacF' factorization failed");

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "psetup failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::psetupB(double t, N_Vector x, N_Vector rx, N_Vector rxdot,
                               booleantype jokB, booleantype *jcurPtrB, double gammaB,
                               void *user_data, N_Vector tmp1B, N_Vector tmp2B,
                               N_Vector tmp3B) {
    try {
      auto m = to_mem(user_data);
      auto& s = m->self;
      // Store gamma for later
      m->gammaB = gammaB;
      // Calculate Jacobian
      double one=1;
      m->arg[0] = &t;
      m->arg[1] = NV_DATA_S(rx);
      m->arg[2] = m->rp;
      m->arg[3] = NV_DATA_S(x);
      m->arg[4] = m->p;
      m->arg[5] = &gammaB;
      m->arg[6] = &one;
      m->res[0] = m->jacB;
      if (s.calc_function(m, "jacB")) casadi_error("'jacB' calculation failed");

      // Prepare the solution of the linear system (e.g. factorize)
      if (s.linsolB_.nfact(m->jacB, m->mem_linsolB)) casadi_error("'jacB' factorization failed");

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "psetupB failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                              booleantype *jcurPtr,
                              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      //auto& s = m->self;

      // Current time
      double t = cv_mem->cv_tn;

      // Scaling factor before J
      double gamma = cv_mem->cv_gamma;

      // Call the preconditioner setup function (which sets up the linear solver)
      if (psetup(t, x, xdot, FALSE, jcurPtr,
                 gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "lsetup failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::lsetupB(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                               booleantype *jcurPtr,
                               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      CVadjMem ca_mem;
      //CVodeBMem cvB_mem;

      int flag;

      // Current time
      double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
      double gamma = cv_mem->cv_gamma;

      cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);

      ca_mem = cv_mem->cv_adj_mem;
      //cvB_mem = ca_mem->ca_bckpbCrt;

      // Get FORWARD solution from interpolation.
      flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, nullptr);
      if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");

      // Call the preconditioner setup function (which sets up the linear solver)
      if (psetupB(t, ca_mem->ca_ytmp, x, xdot, FALSE, jcurPtr,
                  gamma, static_cast<void*>(m), vtemp1, vtemp2, vtemp3)) return 1;

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "lsetupB failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector x, N_Vector xdot) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      //auto& s = m->self;

      // Current time
      double t = cv_mem->cv_tn;

      // Scaling factor before J
      double gamma = cv_mem->cv_gamma;

      // Accuracy
      double delta = 0.0;

      // Left/right preconditioner
      casadi_int lr = 1;

      // Call the preconditioner solve function (which solves the linear system)
      if (psolve(t, x, xdot, b, b, gamma, delta,
                 lr, static_cast<void*>(m), nullptr)) return 1;

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "lsolve failed: " << e.what() << endl;
      return -1;
    }
  }

  int CvodesInterface::lsolveB(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                               N_Vector x, N_Vector xdot) {
    try {
      auto m = to_mem(cv_mem->cv_lmem);
      CVadjMem ca_mem;
      //CVodeBMem cvB_mem;

      int flag;

      // Current time
      double t = cv_mem->cv_tn; // TODO(Joel): is this correct?
      double gamma = cv_mem->cv_gamma;

      cv_mem = static_cast<CVodeMem>(cv_mem->cv_user_data);

      ca_mem = cv_mem->cv_adj_mem;
      //cvB_mem = ca_mem->ca_bckpbCrt;

      // Get FORWARD solution from interpolation.
      flag = ca_mem->ca_IMget(cv_mem, t, ca_mem->ca_ytmp, nullptr);
      if (flag != CV_SUCCESS) casadi_error("Could not interpolate forward states");



      // Accuracy
      double delta = 0.0;

      // Left/right preconditioner
      int lr = 1;

      // Call the preconditioner solve function (which solves the linear system)
      if (psolveB(t, ca_mem->ca_ytmp, x, xdot, b, b, gamma, delta, lr,
                  static_cast<void*>(m), nullptr)) return 1;

      return 0;
    } catch(int flag) { // recoverable error
      return flag;
    } catch(exception& e) { // non-recoverable error
      uerr() << "lsolveB failed: " << e.what() << endl;
      return -1;
    }
  }

  Function CvodesInterface::getJ(bool b) const {
    return oracle_.is_a("SXFunction") ? getJ<SX>(b) : getJ<MX>(b);
  }

  template<typename MatType>
  Function CvodesInterface::getJ(bool backward) const {
    vector<MatType> a = MatType::get_input(oracle_);
    vector<MatType> r = const_cast<Function&>(oracle_)(a);
    MatType c_x = MatType::sym("c_x");
    MatType c_xdot = MatType::sym("c_xdot");

    // Get the Jacobian in the Newton iteration
    if (backward) {
      MatType jac = c_x*MatType::jacobian(r[DE_RODE], a[DE_RX])
                  + c_xdot*MatType::eye(nrx_);
      return Function("jacB",
                      {a[DE_T], a[DE_RX], a[DE_RP],
                       a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
     } else {
      MatType jac = c_x*MatType::jacobian(r[DE_ODE], a[DE_X])
                  + c_xdot*MatType::eye(nx_);
      return Function("jacF", {a[DE_T], a[DE_X], a[DE_P], c_x, c_xdot}, {jac});
    }
  }

  CvodesMemory::CvodesMemory(const CvodesInterface& s) : self(s) {
    this->mem = nullptr;

    // Reset checkpoints counter
    this->ncheck = 0;
  }

  CvodesMemory::~CvodesMemory() {
    if (this->mem) CVodeFree(&this->mem);
  }

} // namespace casadi
