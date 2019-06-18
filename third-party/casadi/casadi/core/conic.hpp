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


#ifndef CASADI_CONIC_HPP
#define CASADI_CONIC_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_conic
      Create a QP solver
      Solves the following strictly convex problem:

      \verbatim
      min          1/2 x' H x + g' x
      x

      subject to
      LBA <= A x <= UBA
      LBX <= x   <= UBX

      resize(Q x, np, np) + P >= 0 (psd)

      with :
      H sparse (n x n) positive definite
      g dense  (n x 1)
      A sparse (nc x n)
      Q sparse symmetric (np^2 x n)
      P sparse symmetric (np x nq)

      n: number of decision variables (x)
      nc: number of constraints (A)
      nq: shape of psd constraint matrix

      \endverbatim

      If H is not positive-definite, the solver should throw an error.

      Second-order cone constraints can be added as psd constraints
      through a helper function 'soc':

      x in R^n
      y in R

      || x ||_2 <= y

        <=>

      soc(x, y) psd

      This can be proven with soc(x, y)=[y*I   x; x'  y]
      using the Shur complement.


      \generalsection{Conic}
      \pluginssection{Conic}

      \author Joel Andersson
      \date 2011-2015
  */

  /** \defgroup conic
  * @copydoc main_conic
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_conic
  * \endif
  */
  ///@{
  CASADI_EXPORT Function conic(const std::string& name, const std::string& solver,
                               const SpDict& qp, const Dict& opts=Dict());
  CASADI_EXPORT Function qpsol(const std::string& name, const std::string& solver,
                               const SXDict& qp, const Dict& opts=Dict());
  CASADI_EXPORT Function qpsol(const std::string& name, const std::string& solver,
                               const MXDict& qp, const Dict& opts=Dict());
  ///@}

  /** \brief Get input scheme of QP solvers */
  CASADI_EXPORT std::vector<std::string> conic_in();

  /** \brief Get QP solver output scheme of QP solvers */
  CASADI_EXPORT std::vector<std::string> conic_out();

  /** \brief Get QP solver input scheme name by index */
  CASADI_EXPORT std::string conic_in(casadi_int ind);

  /** \brief Get output scheme name by index */
  CASADI_EXPORT std::string conic_out(casadi_int ind);

  /** \brief Get the number of QP solver inputs */
  CASADI_EXPORT casadi_int conic_n_in();

  /** \brief Get the number of QP solver outputs */
  CASADI_EXPORT casadi_int conic_n_out();

  /** \brief Get all options for a plugin */
  CASADI_EXPORT std::vector<std::string> conic_options(const std::string& name);

  /** \brief Get type info for a particular option */
  CASADI_EXPORT std::string conic_option_type(const std::string& name, const std::string& op);

  /** \brief Get documentation for a particular option */
  CASADI_EXPORT std::string conic_option_info(const std::string& name, const std::string& op);

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_conic(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_conic(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_conic(const std::string& name);

  /** Generate native code in the interfaced language for debugging */
  CASADI_EXPORT void conic_debug(const Function& f, const std::string &filename);

  /** Generate native code in the interfaced language for debugging */
  CASADI_EXPORT void conic_debug(const Function& f, std::ostream &file);

  /** @} */

#ifndef SWIG
/// Input arguments of a QP problem
enum ConicInput {
  /// The square matrix H: sparse, (n x n). Only the lower triangular part is actually used.
  /// The matrix is assumed to be symmetrical.
  CONIC_H,
  /// The vector g: dense,  (n x 1)
  CONIC_G,
  /// The matrix A: sparse, (nc x n) - product with x must be dense.
  CONIC_A,
  /// dense, (nc x 1)
  CONIC_LBA,
  /// dense, (nc x 1)
  CONIC_UBA,
  /// dense, (n x 1)
  CONIC_LBX,
  /// dense, (n x 1)
  CONIC_UBX,
  /// dense, (n x 1)
  CONIC_X0,
  /// dense
  CONIC_LAM_X0,
  /// dense
  CONIC_LAM_A0,
  /// The matrix Q: sparse symmetric, (np^2 x n)
  CONIC_Q,
  /// The matrix P: sparse symmetric, (np x np)
  CONIC_P,
  CONIC_NUM_IN};

/// Output arguments of an QP Solver
enum ConicOutput {
  /// The primal solution
  CONIC_X,
  /// The optimal cost
  CONIC_COST,
  /// The dual solution corresponding to linear bounds
  CONIC_LAM_A,
  /// The dual solution corresponding to simple bounds
  CONIC_LAM_X,
  CONIC_NUM_OUT};
#endif // SWIG

} // namespace casadi

#endif // CASADI_CONIC_HPP
