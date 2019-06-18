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


#ifndef CASADI_NLPSOL_HPP
#define CASADI_NLPSOL_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_nlpsol
      Create an NLP solver
      Creates a solver for the following parametric nonlinear program (NLP):
      \verbatim

      min          F(x, p)
      x

      subject to
      LBX <=   x    <= UBX
      LBG <= G(x, p) <= UBG
      p  == P

      nx: number of decision variables
      ng: number of constraints
      np: number of parameters

      \endverbatim

      \generalsection{Nlpsol}
      \pluginssection{Nlpsol}

      \author Joel Andersson
      \date 2011-2015
  */

  /** \defgroup nlpsol
  * @copydoc main_nlpsol
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_nlpsol
  * \endif
  */
  ///@{
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const SXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const MXDict& nlp, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const std::string& fname, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Importer& compiler, const Dict& opts=Dict());
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const NlpBuilder& nl, const Dict& opts=Dict());
#ifndef SWIG
  CASADI_EXPORT Function nlpsol(const std::string& name, const std::string& solver,
                                const Function& nlp, const Dict& opts=Dict());
#endif // SWIG
  ///@}

  /** \brief Get input scheme of NLP solvers
  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif
  */
  CASADI_EXPORT std::vector<std::string> nlpsol_in();

  /** \brief Get NLP solver output scheme of NLP solvers
  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif
  */
  CASADI_EXPORT std::vector<std::string> nlpsol_out();

  /** \brief Get NLP solver input scheme name by index
  * \if EXPANDED
  * @copydoc scheme_NlpsolInput
  * \endif
  */
  CASADI_EXPORT std::string nlpsol_in(casadi_int ind);

  /** \brief Get output scheme name by index
  * \if EXPANDED
  * @copydoc scheme_NlpsolOutput
  * \endif
  */
  CASADI_EXPORT std::string nlpsol_out(casadi_int ind);

  /** \brief Number of NLP solver inputs */
  CASADI_EXPORT casadi_int nlpsol_n_in();

  /** \brief Number of NLP solver outputs */
  CASADI_EXPORT casadi_int nlpsol_n_out();

  ///@{
  /** \brief Default input for an NLP solver */
  CASADI_EXPORT double nlpsol_default_in(casadi_int ind);
  CASADI_EXPORT std::vector<double> nlpsol_default_in();
  ///@}

  /** \brief Get all options for a plugin */
  CASADI_EXPORT std::vector<std::string> nlpsol_options(const std::string& name);

  /** \brief Get type info for a particular option */
  CASADI_EXPORT std::string nlpsol_option_type(const std::string& name, const std::string& op);

  /** \brief Get documentation for a particular option */
  CASADI_EXPORT std::string nlpsol_option_info(const std::string& name, const std::string& op);

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_nlpsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_nlpsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_nlpsol(const std::string& name);

  /** @} */

#ifndef SWIG
/// Input arguments of an NLP function
enum NLPInput {
  /// Decision variable
  NL_X,
  /// Fixed parameter
  NL_P,
  /// Number of NLP inputs
  NL_NUM_IN
};

/// Shortname for onput arguments of an NLP function
const std::vector<std::string> NL_INPUTS = {"x", "p"};

/// Output arguments of an NLP function
enum NLPOutput {
  /// Objective function
  NL_F,
  /// Constraint function
  NL_G,
  /// Number of NLP outputs
  NL_NUM_OUT
};

/// Shortname for output arguments of an NLP function
const std::vector<std::string> NL_OUTPUTS = {"f", "g"};

/// Input arguments of an NLP Solver
enum NlpsolInput {
  /// Decision variables, initial guess (nx x 1)
  NLPSOL_X0,
  /// Value of fixed parameters (np x 1)
  NLPSOL_P,
  /// Decision variables lower bound (nx x 1), default -inf
  NLPSOL_LBX,
  /// Decision variables upper bound (nx x 1), default +inf
  NLPSOL_UBX,
  /// Constraints lower bound (ng x 1), default -inf
  NLPSOL_LBG,
  /// Constraints upper bound (ng x 1), default +inf
  NLPSOL_UBG,
  /// Lagrange multipliers for bounds on X, initial guess (nx x 1)
  NLPSOL_LAM_X0,
  /// Lagrange multipliers for bounds on G, initial guess (ng x 1)
  NLPSOL_LAM_G0,
  NLPSOL_NUM_IN
};

/// Output arguments of an NLP Solver
enum NlpsolOutput {
  /// Decision variables at the optimal solution (nx x 1)
  NLPSOL_X,
  /// Cost function value at the optimal solution (1 x 1)
  NLPSOL_F,
  /// Constraints function at the optimal solution (ng x 1)
  NLPSOL_G,
  /// Lagrange multipliers for bounds on X at the solution (nx x 1)
  NLPSOL_LAM_X,
  /// Lagrange multipliers for bounds on G at the solution (ng x 1)
  NLPSOL_LAM_G,
  /// Lagrange multipliers for bounds on P at the solution (np x 1)
  NLPSOL_LAM_P,
  NLPSOL_NUM_OUT
};
#endif // SWIG

} // namespace casadi

#endif // CASADI_NLPSOL_HPP
