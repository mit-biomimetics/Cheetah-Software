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

#include "generic_matrix.hpp"

namespace casadi {


  double index_interp1d(const std::vector<double>& x, double xq, bool equidistant) {

    if (equidistant) {
      double delta = x[1]-x[0];
      double i = (xq-x[0])/delta;
      double imax = static_cast<double>(x.size()-1);
      return std::max(std::min(i, imax), 0.0);

    } else {
      std::vector<double>::const_iterator it = std::lower_bound(x.begin(), x.end(), xq);

      // End of x
      if (it==x.end()) return static_cast<double>(x.size()-1);

      // Start of x
      if (it==x.begin()) return 0.0;

      casadi_int i = std::distance(x.begin(), it);

      // Exactly on an entry
      if (*it == xq) return static_cast<double>(i);

      double b = *it;
      double a = *(it-1);

      return static_cast<double>(i)+(xq-b)/(b-a);
    }
  }

} // namespace casadi
