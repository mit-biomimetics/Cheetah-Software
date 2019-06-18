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


#ifndef CASADI_INTEGRATOR_HPP
#define CASADI_INTEGRATOR_HPP

#include "function.hpp"
#include "linsol.hpp"
#include "rootfinder.hpp"

namespace casadi {

  /** \defgroup main_integrator
      Create an ODE/DAE integrator
      Solves an initial value problem (IVP) coupled to a terminal value problem
      with differential equation given as an implicit ODE coupled to an algebraic
      equation and a set of quadratures:

      \verbatim
      Initial conditions at t=t0
      x(t0)  = x0
      q(t0)  = 0

      Forward integration from t=t0 to t=tf
      der(x) = function(x, z, p, t)                  Forward ODE
      0 = fz(x, z, p, t)                  Forward algebraic equations
      der(q) = fq(x, z, p, t)                  Forward quadratures

      Terminal conditions at t=tf
      rx(tf)  = rx0
      rq(tf)  = 0

      Backward integration from t=tf to t=t0
      der(rx) = gx(rx, rz, rp, x, z, p, t)        Backward ODE
      0 = gz(rx, rz, rp, x, z, p, t)        Backward algebraic equations
      der(rq) = gq(rx, rz, rp, x, z, p, t)        Backward quadratures

      where we assume that both the forward and backwards integrations are index-1
      (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
      gx, gz and gq have a linear dependency on rx, rz and rp.

      \endverbatim

      \generalsection{Integrator}
      \pluginssection{Integrator}

      \author Joel Andersson
      \date 2011-2015
  */
  /** \defgroup integrator
  * @copydoc main_integrator
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_integrator
  * \endif
  */
  ///@{
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const SXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const MXDict& dae, const Dict& opts=Dict());
#ifndef SWIG
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const Function& dae, const Dict& opts=Dict());
#endif // SWIG
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_integrator(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_integrator(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_integrator(const std::string& name);

  /** \brief Get input scheme of integrators */
  CASADI_EXPORT std::vector<std::string> integrator_in();

  /** \brief Get integrator output scheme of integrators */
  CASADI_EXPORT std::vector<std::string> integrator_out();

  /** \brief Get integrator input scheme name by index */
  CASADI_EXPORT std::string integrator_in(casadi_int ind);

  /** \brief Get output scheme name by index */
  CASADI_EXPORT std::string integrator_out(casadi_int ind);

  /** \brief Get the number of integrator inputs */
  CASADI_EXPORT casadi_int integrator_n_in();

  /** \brief Get the number of integrator outputs */
  CASADI_EXPORT casadi_int integrator_n_out();
  /** @} */

#ifndef SWIG
/// Inputs of the symbolic representation of the DAE
enum DeIn {
  DE_T,
  DE_X,
  DE_Z,
  DE_P,
  DE_RX,
  DE_RZ,
  DE_RP,
  DE_NUM_IN};

/// Shortnames for DAE symbolic representation inputs
const std::vector<std::string> DE_INPUTS = {"t", "x", "z", "p", "rx", "rz", "rp"};

/// Inputs of the symbolic representation of the DAE
enum DeOut {
  DE_ODE,
  DE_ALG,
  DE_QUAD,
  DE_RODE,
  DE_RALG,
  DE_RQUAD,
  DE_NUM_OUT};

/// Shortnames for DAE symbolic representation outputs
const std::vector<std::string> DE_OUTPUTS = {"ode", "alg", "quad", "rode", "ralg", "rquad"};

/// Input arguments of an ODE/DAE function
enum DAEInput {
  /// Differential state
  DAE_X,
  /// Algebraic state
  DAE_Z,
  /// Parameter
  DAE_P,
  /// Explicit time dependence
  DAE_T,
  /// Number of arguments
  DAE_NUM_IN
};

/// Output arguments of an DAE function
enum DAEOutput {
  /// Right hand side of the implicit ODE
  DAE_ODE,
  /// Right hand side of algebraic equations
  DAE_ALG,
  /// Right hand side of quadratures equations
  DAE_QUAD,
  /// Number of arguments
  DAE_NUM_OUT
};

/// Input arguments of an ODE/DAE backward integration function
enum RDAEInput {
  /// Backward differential state
  RDAE_RX,
  /// Backward algebraic state
  RDAE_RZ,
  /// Backward  parameter vector
  RDAE_RP,
  /// Forward differential state
  RDAE_X,
  /// Forward algebraic state
  RDAE_Z,
  /// Parameter vector
  RDAE_P,
  /// Explicit time dependence
  RDAE_T,
  /// Number of arguments
  RDAE_NUM_IN
};

/// Output arguments of an ODE/DAE backward integration function
enum RDAEOutput {
  /// Right hand side of ODE
  RDAE_ODE,
  /// Right hand side of algebraic equations
  RDAE_ALG,
  /// Right hand side of quadratures
  RDAE_QUAD,
  /// Number of arguments
  RDAE_NUM_OUT
};

/// Input arguments of an integrator
enum IntegratorInput {
  /// Differential state at the initial time
  INTEGRATOR_X0,
  /// Parameters
  INTEGRATOR_P,
  /// Initial guess for the algebraic variable
  INTEGRATOR_Z0,
  /// Backward differential state at the final time
  INTEGRATOR_RX0,
  /// Backward parameter vector
  INTEGRATOR_RP,
  /// Initial guess for the backwards algebraic variable
  INTEGRATOR_RZ0,
  /// Number of input arguments of an integrator
  INTEGRATOR_NUM_IN
};

/// Output arguments of an integrator
enum IntegratorOutput {
  /// Differential state at the final time
  INTEGRATOR_XF,
  /// Quadrature state at the final time
  INTEGRATOR_QF,
  /// Algebraic variable at the final time
  INTEGRATOR_ZF,
  /// Backward differential state at the initial time
  INTEGRATOR_RXF,
  /// Backward quadrature state at the initial time
  INTEGRATOR_RQF,
  /// Backward algebraic variable at the initial time
  INTEGRATOR_RZF,
  /// Number of output arguments of an integrator
  INTEGRATOR_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
