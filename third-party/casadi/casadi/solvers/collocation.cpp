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


#include "collocation.hpp"
#include "casadi/core/polynomial.hpp"
#include "casadi/core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_COLLOCATION_EXPORT
      casadi_register_integrator_collocation(Integrator::Plugin* plugin) {
    plugin->creator = Collocation::creator;
    plugin->name = "collocation";
    plugin->doc = Collocation::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Collocation::options_;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_COLLOCATION_EXPORT casadi_load_integrator_collocation() {
    Integrator::registerPlugin(casadi_register_integrator_collocation);
  }

  Collocation::Collocation(const std::string& name, const Function& dae)
    : ImplicitFixedStepIntegrator(name, dae) {
  }

  Collocation::~Collocation() {
  }

  Options Collocation::options_
  = {{&ImplicitFixedStepIntegrator::options_},
     {{"interpolation_order",
       {OT_INT,
        "Order of the interpolating polynomials"}},
      {"collocation_scheme",
       {OT_STRING,
        "Collocation scheme: radau|legendre"}}
     }
  };

  void Collocation::init(const Dict& opts) {
    // Default options
    deg_ = 3;
    collocation_scheme_ = "radau";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="interpolation_order") {
        deg_ = op.second;
      } else if (op.first=="collocation_scheme") {
        collocation_scheme_ = op.second.to_string();
      }
    }

    // Call the base class init
    ImplicitFixedStepIntegrator::init(opts);
  }

  void Collocation::setupFG() {
    f_ = create_function("f", {"x", "z", "p", "t"}, {"ode", "alg", "quad"});
    g_ = create_function("g", {"rx", "rz", "rp", "x", "z", "p", "t"},
                              {"rode", "ralg", "rquad"});

    // All collocation time points
    std::vector<double> tau_root = collocation_points(deg_, collocation_scheme_);
    tau_root.insert(tau_root.begin(), 0);

    // Coefficients of the collocation equation
    vector<vector<double> > C(deg_+1, vector<double>(deg_+1, 0));

    // Coefficients of the continuity equation
    vector<double> D(deg_+1, 0);

    // Coefficients of the quadratures
    vector<double> B(deg_+1, 0);

    // For all collocation points
    for (casadi_int j=0; j<deg_+1; ++j) {

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for (casadi_int r=0; r<deg_+1; ++r) {
        if (r!=j) {
          p *= Polynomial(-tau_root[r], 1)/(tau_root[j]-tau_root[r]);
        }
      }

      // Evaluate the polynomial at the final time to get the
      // coefficients of the continuity equation
      if (collocation_scheme_=="radau") {
        D[j] = j==deg_ ? 1 : 0;
      } else {
        D[j] = p(1.0);
      }

      // Evaluate the time derivative of the polynomial at all collocation points to
      // get the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for (casadi_int r=0; r<deg_+1; ++r) {
        C[j][r] = dp(tau_root[r]);
      }

      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      B[j] = ip(1.0);
    }

    // Symbolic inputs
    MX x0 = MX::sym("x0", this->x());
    MX p = MX::sym("p", this->p());
    MX t = MX::sym("t", this->t());

    // Implicitly defined variables (z and x)
    MX v = MX::sym("v", deg_*(nx_+nz_));
    vector<casadi_int> v_offset(1, 0);
    for (casadi_int d=0; d<deg_; ++d) {
      v_offset.push_back(v_offset.back()+nx_);
      v_offset.push_back(v_offset.back()+nz_);
    }
    vector<MX> vv = vertsplit(v, v_offset);
    vector<MX>::const_iterator vv_it = vv.begin();

    // Collocated states
    vector<MX> x(deg_+1), z(deg_+1);
    for (casadi_int d=1; d<=deg_; ++d) {
      x[d] = reshape(*vv_it++, size_in(INTEGRATOR_X0));
      z[d] = reshape(*vv_it++, size_in(INTEGRATOR_Z0));
    }
    casadi_assert_dev(vv_it==vv.end());

    // Collocation time points
    vector<MX> tt(deg_+1);
    for (casadi_int d=0; d<=deg_; ++d) {
      tt[d] = t + h_*tau_root[d];
    }

    // Equations that implicitly define v
    vector<MX> eq;

    // Quadratures
    MX qf = MX::zeros(this->q());

    // End state
    MX xf = D[0]*x0;

    // For all collocation points
    for (casadi_int j=1; j<deg_+1; ++j) {
      //for (casadi_int j=deg_; j>=1; --j) {

      // Evaluate the DAE
      vector<MX> f_arg(DAE_NUM_IN);
      f_arg[DAE_T] = tt[j];
      f_arg[DAE_P] = p;
      f_arg[DAE_X] = x[j];
      f_arg[DAE_Z] = z[j];
      vector<MX> f_res = f_(f_arg);

      // Get an expression for the state derivative at the collocation point
      MX xp_j = C[0][j] * x0;
      for (casadi_int r=1; r<deg_+1; ++r) {
        xp_j += C[r][j] * x[r];
      }

      // Add collocation equation
      eq.push_back(vec(h_*f_res[DAE_ODE] - xp_j));

      // Add the algebraic conditions
      eq.push_back(vec(f_res[DAE_ALG]));

      // Add contribution to the final state
      xf += D[j]*x[j];

      // Add contribution to quadratures
      qf += (B[j]*h_)*f_res[DAE_QUAD];
    }

    // Form forward discrete time dynamics
    vector<MX> F_in(DAE_NUM_IN);
    F_in[DAE_T] = t;
    F_in[DAE_X] = x0;
    F_in[DAE_P] = p;
    F_in[DAE_Z] = v;
    vector<MX> F_out(DAE_NUM_OUT);
    F_out[DAE_ODE] = xf;
    F_out[DAE_ALG] = vertcat(eq);
    F_out[DAE_QUAD] = qf;
    F_ = Function("dae", F_in, F_out);
    alloc(F_);

    // Backwards dynamics
    // NOTE: The following is derived so that it will give the exact adjoint
    // sensitivities whenever g is the reverse mode derivative of f.
    if (!g_.is_null()) {

      // Symbolic inputs
      MX rx0 = MX::sym("rx0", this->rx());
      MX rp = MX::sym("rp", this->rp());

      // Implicitly defined variables (rz and rx)
      MX rv = MX::sym("v", deg_*(nrx_+nrz_));
      vector<casadi_int> rv_offset(1, 0);
      for (casadi_int d=0; d<deg_; ++d) {
        rv_offset.push_back(rv_offset.back()+nrx_);
        rv_offset.push_back(rv_offset.back()+nrz_);
      }
      vector<MX> rvv = vertsplit(rv, rv_offset);
      vector<MX>::const_iterator rvv_it = rvv.begin();

      // Collocated states
      vector<MX> rx(deg_+1), rz(deg_+1);
      for (casadi_int d=1; d<=deg_; ++d) {
        rx[d] = reshape(*rvv_it++, this->rx().size());
        rz[d] = reshape(*rvv_it++, this->rz().size());
      }
      casadi_assert_dev(rvv_it==rvv.end());

      // Equations that implicitly define v
      eq.clear();

      // Quadratures
      MX rqf = MX::zeros(this->rq());

      // End state
      MX rxf = D[0]*rx0;

      // For all collocation points
      for (casadi_int j=1; j<deg_+1; ++j) {

        // Evaluate the backward DAE
        vector<MX> g_arg(RDAE_NUM_IN);
        g_arg[RDAE_T] = tt[j];
        g_arg[RDAE_P] = p;
        g_arg[RDAE_X] = x[j];
        g_arg[RDAE_Z] = z[j];
        g_arg[RDAE_RX] = rx[j];
        g_arg[RDAE_RZ] = rz[j];
        g_arg[RDAE_RP] = rp;
        vector<MX> g_res = g_(g_arg);

        // Get an expression for the state derivative at the collocation point
        MX rxp_j = -D[j]*rx0;
        for (casadi_int r=1; r<deg_+1; ++r) {
          rxp_j += (B[r]*C[j][r]) * rx[r];
        }

        // Add collocation equation
        eq.push_back(vec(h_*B[j]*g_res[RDAE_ODE] - rxp_j));

        // Add the algebraic conditions
        eq.push_back(vec(g_res[RDAE_ALG]));

        // Add contribution to the final state
        rxf += -B[j]*C[0][j]*rx[j];

        // Add contribution to quadratures
        rqf += h_*B[j]*g_res[RDAE_QUAD];
      }

      // Form backward discrete time dynamics
      vector<MX> G_in(RDAE_NUM_IN);
      G_in[RDAE_T] = t;
      G_in[RDAE_X] = x0;
      G_in[RDAE_P] = p;
      G_in[RDAE_Z] = v;
      G_in[RDAE_RX] = rx0;
      G_in[RDAE_RP] = rp;
      G_in[RDAE_RZ] = rv;
      vector<MX> G_out(RDAE_NUM_OUT);
      G_out[RDAE_ODE] = rxf;
      G_out[RDAE_ALG] = vertcat(eq);
      G_out[RDAE_QUAD] = rqf;
      G_ = Function("rdae", G_in, G_out);
      alloc(G_);
    }
  }

  void Collocation::reset(IntegratorMemory* mem, double t, const double* x,
                                const double* z, const double* p) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Reset the base classes
    ImplicitFixedStepIntegrator::reset(mem, t, x, z, p);

    // Initial guess for Z
    double* Z = get_ptr(m->Z);
    for (casadi_int d=0; d<deg_; ++d) {
      casadi_copy(x, nx_, Z);
      Z += nx_;
      casadi_copy(z, nz_, Z);
      Z += nz_;
    }
  }

  void Collocation::resetB(IntegratorMemory* mem, double t, const double* rx,
                               const double* rz, const double* rp) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Reset the base classes
    ImplicitFixedStepIntegrator::resetB(mem, t, rx, rz, rp);

    // Initial guess for RZ
    double* RZ = get_ptr(m->RZ);
    for (casadi_int d=0; d<deg_; ++d) {
      casadi_copy(rx, nrx_, RZ);
      RZ += nrx_;
      casadi_copy(rz, nrz_, RZ);
      RZ += nrz_;
    }
  }

} // namespace casadi
