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


#ifndef CASADI_INTERPOLANT_HPP
#define CASADI_INTERPOLANT_HPP

#include "function.hpp"

namespace casadi {


  /** \defgroup main_interpolant
   * An interpolant function for lookup table data
   *
   * \param[in] name label for the resulting Function
   * \param[in] solver name of the plugin
   * \param[in] grid collection of 1D grids whose outer product
   *            defines the full N-D rectangular grid
   * \param[in] values flattened vector of all values
   *            for all gridpoints
   *
   * Syntax 1D
   * \verbatim
   * # Python
   * xgrid = np.linspace(1,6,6)
   * V = [-1,-1,-2,-3,0,2]
   * LUT = casadi.interpolant("LUT","bspline",[xgrid],V)
   * print(LUT(2.5))
   * \endverbatim
   * \verbatim
   * % Matlab
   * xgrid = 1:6;
   * V = [-1 -1 -2 -3 0 2];
   * LUT = casadi.interpolant('LUT','bspline',{xgrid},V);
   * LUT(2.5)
   * \endverbatim
   *
   * Syntax 2D
   * \verbatim
   * # Python
   * xgrid = np.linspace(-5,5,11)
   * ygrid = np.linspace(-4,4,9)
   * X,Y = np.meshgrid(xgrid,ygrid,indexing='ij')
   * R = np.sqrt(5*X**2 + Y**2)+ 1
   * data = np.sin(R)/R
   * data_flat = data.ravel(order='F')
   * LUT = casadi.interpolant('name','bspline',[xgrid,ygrid],data_flat)
   * print(LUT([0.5,1]))
   * \enverbatim
   * \verbatim
   * % Matlab
   * xgrid = -5:1:5;
   * ygrid = -4:1:4;
   * R = sqrt(5*X.^2 + Y.^2)+ 1;
   * V = sin(R)./(R);
   * LUT = interpolant('LUT','bspline',{xgrid, ygrid},V(:));
   * LUT([0.5 1])
   * \endverbatim
   *
   *  \generalsection{Interpolant}
   *  \pluginssection{Interpolant}
   * \author Joel Andersson
   * \date 2016
   */

  /** \defgroup interpolant
  * @copydoc main_interpolant
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_interpolant
  * \endif
  */
  ///@{
  CASADI_EXPORT Function interpolant(const std::string& name,
                                     const std::string& solver,
                                     const std::vector<std::vector<double> >& grid,
                                     const std::vector<double>& values,
                                     const Dict& opts=Dict());
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_interpolant(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_interpolant(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_interpolant(const std::string& name);

  /** @} */

} // namespace casadi

#endif // CASADI_INTERPOLANT_HPP
