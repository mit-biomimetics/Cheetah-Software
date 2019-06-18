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


#ifndef CASADI_DPLE_HPP
#define CASADI_DPLE_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_dple

      Discrete periodic Lyapunov Equation solver
      Given matrices \f$A_k\f$ and symmetric \f$V_k,  k = 0..K-1\f$
      \verbatim
      A_k in R^(n x n)
      V_k in R^n
      \endverbatim
      provides all of \f$P_k\f$ that satisfy:
      \verbatim
      P_0 = A_(K-1)*P_(K-1)*A_(K-1)' + V_k
      P_k+1 = A_k*P_k*A_k' + V_k  for k = 1..K-1
      \endverbatim

      \generalsection{Dple}
      \pluginssection{Dple}

      \author Joris Gillis
      \date 2013-2016
  */

  /** \defgroup dple
  * @copydoc main_dple
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_dple
  * \endif
  */
  ///@{
  CASADI_EXPORT Function dplesol(const std::string& name, const std::string& solver,
                           const SpDict& qp, const Dict& opts=Dict());
  CASADI_EXPORT MX dplesol(const MX& A, const MX& V, const std::string& solver,
    const Dict& opts=Dict());
  CASADI_EXPORT MXVector dplesol(const MXVector& A, const MXVector& V, const std::string& solver,
    const Dict& opts=Dict());
  CASADI_EXPORT DMVector dplesol(const DMVector& A, const DMVector& V, const std::string& solver,
    const Dict& opts=Dict());
  ///@}

  /** \brief Get input scheme of DPLE solvers */
  CASADI_EXPORT std::vector<std::string> dple_in();

  /** \brief Get output scheme of DPLE solvers */
  CASADI_EXPORT std::vector<std::string> dple_out();

  /** \brief Get DPLE input scheme name by index */
  CASADI_EXPORT std::string dple_in(casadi_int ind);

  /** \brief Get DPLE output scheme name by index */
  CASADI_EXPORT std::string dple_out(casadi_int ind);

  /** \brief Get the number of QP solver inputs */
  CASADI_EXPORT casadi_int dple_n_in();

  /** \brief Get the number of QP solver outputs */
  CASADI_EXPORT casadi_int dple_n_out();

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_dple(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_dple(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_dple(const std::string& name);

  /** @} */

#ifndef SWIG
/// Input arguments of a \e dple solver [dpleIn]
enum DpleInput {
  /// A matrices (horzcat when const_dim, diagcat otherwise) [a]
  DPLE_A,
  /// V matrices (horzcat when const_dim, diagcat otherwise) [v]
  DPLE_V,
  DPLE_NUM_IN
};

/// Output arguments of a \e dple solver [dpleOut]
enum DpleOutput {
  /// Lyapunov matrix (horzcat when const_dim, diagcat otherwise) (Cholesky of P if pos_def) [p]
  DPLE_P,
  /// Number of arguments.
  DPLE_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_DPLE_HPP
