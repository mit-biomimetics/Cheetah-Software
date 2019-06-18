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


#include "runge_kutta.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_RK_EXPORT
      casadi_register_integrator_rk(Integrator::Plugin* plugin) {
    plugin->creator = RungeKutta::creator;
    plugin->name = "rk";
    plugin->doc = RungeKutta::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &RungeKutta::options_;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_RK_EXPORT casadi_load_integrator_rk() {
    Integrator::registerPlugin(casadi_register_integrator_rk);
  }

  RungeKutta::RungeKutta(const std::string& name, const Function& dae)
    : FixedStepIntegrator(name, dae) {
  }

  RungeKutta::~RungeKutta() {
  }

  void RungeKutta::init(const Dict& opts) {
    // Call the base class init
    FixedStepIntegrator::init(opts);

    // Algebraic variables not (yet?) supported
    casadi_assert(nz_==0 && nrz_==0,
                          "Explicit Runge-Kutta integrators do not support algebraic variables");
  }

  void RungeKutta::setupFG() {
    f_ = create_function("f", {"x", "z", "p", "t"}, {"ode", "alg", "quad"});
    g_ = create_function("g", {"rx", "rz", "rp", "x", "z", "p", "t"},
                              {"rode", "ralg", "rquad"});

    // Symbolic inputs
    MX x0 = MX::sym("x0", this->x());
    MX p = MX::sym("p", this->p());
    MX t = MX::sym("t", this->t());

    // Intermediate variables (does not enter in F_, only in G_)
    MX v = MX::sym("v", x0.size1(), x0.size2()*3);
    vector<MX> x = horzsplit(v, x0.size2());
    casadi_assert_dev(x.size()==3);

    // Definitions of x
    vector<MX> x_def(3);

    // Time points
    vector<MX> tt(3);

    // Forward integration
    {
      // Arguments when calling f
      vector<MX> f_arg(DAE_NUM_IN);
      vector<MX> f_res;
      f_arg[DAE_P] = p;

      // k1
      f_arg[DAE_T] = t;
      f_arg[DAE_X] = x0;
      f_res = f_(f_arg);
      MX k1 = f_res[DAE_ODE];
      MX k1q = f_res[DAE_QUAD];

      // k2
      tt[0] = f_arg[DAE_T] = t + h_/2;
      x_def[0] = f_arg[DAE_X] = x0 + (h_/2) * k1;
      f_res = f_(f_arg);
      MX k2 = f_res[DAE_ODE];
      MX k2q = f_res[DAE_QUAD];

      // k3
      tt[1] = tt[0];
      x_def[1] = f_arg[DAE_X] = x0 + (h_/2) * k2;
      f_res = f_(f_arg);
      MX k3 = f_res[DAE_ODE];
      MX k3q = f_res[DAE_QUAD];

      // k4
      tt[2] = f_arg[DAE_T] = t + h_;
      x_def[2] = f_arg[DAE_X] = x0 + h_ * k3;
      f_res = f_(f_arg);
      MX k4 = f_res[DAE_ODE];
      MX k4q = f_res[DAE_QUAD];

      // Take step
      MX xf = x0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
      MX qf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      // TODO(Joel): Change this and make x_k1, x_k2, x_k3 and x_k4 algebraic outputs
      f_arg[DAE_T] = t;
      f_arg[DAE_X] = x0;
      f_arg[DAE_P] = p;
      f_arg[DAE_Z] = v;
      f_res[DAE_ODE] = xf;
      f_res[DAE_QUAD] = qf;
      f_res[DAE_ALG] = horzcat(x_def);
      F_ = Function("dae", f_arg, f_res);
      alloc(F_);
    }

    // Backward integration
    if (!g_.is_null()) {
      // Symbolic inputs
      MX rx0 = MX::sym("rx0", this->rx());
      MX rp = MX::sym("rp", this->rp());

      // Intermediate variables (do not enter in G_)
      MX rv = MX::sym("rv", rx0.size1(), 3*rx0.size2());
      vector<MX> rx_def(3);

      // Arguments when calling g
      vector<MX> g_arg(RDAE_NUM_IN);
      vector<MX> g_res;
      g_arg[RDAE_P] = p;
      g_arg[RDAE_RP] = rp;

      // k1
      g_arg[RDAE_T] = tt[2];
      g_arg[RDAE_X] = x[2];
      g_arg[RDAE_RX] = rx0;
      g_res = g_(g_arg);
      MX k1 = g_res[RDAE_ODE];
      MX k1q = g_res[RDAE_QUAD];

      // k2
      g_arg[RDAE_T] = tt[1];
      g_arg[RDAE_X] = x[1];
      g_arg[RDAE_RX] = rx_def[2] = rx0 + (h_/2) * k1;
      g_res = g_(g_arg);
      MX  k2 = g_res[RDAE_ODE];
      MX  k2q = g_res[RDAE_QUAD];

      // k3
      g_arg[RDAE_T] = tt[0];
      g_arg[RDAE_X] = x[0];
      g_arg[RDAE_RX] = rx_def[1] = rx0 + (h_/2) * k2;
      g_res = g_(g_arg);
      MX k3 = g_res[RDAE_ODE];
      MX k3q = g_res[RDAE_QUAD];

      // k4
      g_arg[RDAE_T] = t;
      g_arg[RDAE_X] = x0;
      g_arg[RDAE_RX] = rx_def[0] = rx0 + h_ * k3;
      g_res = g_(g_arg);
      MX k4 = g_res[RDAE_ODE];
      MX  k4q = g_res[RDAE_QUAD];

      // Take step
      MX rxf = rx0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
      MX rqf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      g_arg[RDAE_T] = t;
      g_arg[RDAE_X] = x0;
      g_arg[RDAE_P] = p;
      g_arg[RDAE_Z] = v;
      g_arg[RDAE_RX] = rx0;
      g_arg[RDAE_RP] = rp;
      g_arg[RDAE_RZ] = rv;
      g_res[RDAE_ODE] = rxf;
      g_res[RDAE_QUAD] = rqf;
      g_res[RDAE_ALG] = horzcat(rx_def);
      G_ = Function("rdae", g_arg, g_res);
      alloc(G_);
    }
  }

} // namespace casadi
