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


#ifndef CASADI_CASADI_LIMITS_HPP
#define CASADI_CASADI_LIMITS_HPP

#include <cmath>
#include <limits>

/** \brief The casadi namespace */
namespace casadi {

  /** \brief casadi_limits class

  The following class, which acts as a complements to the standard numeric_limits class, allows
  specifying certain properties of scalar objects. The template can be specialized for
  e.g. symbolic scalars
  \author Joel Andersson
  \date 2011
  */

  template<class T>
  class CASADI_EXPORT casadi_limits {
    public:
      static bool is_zero(const T& val) {
        return val==0;
      }
      static bool is_equal(const T& x, const T& y, casadi_int depth) {
        return x==y;
      }
      static bool is_almost_zero(const T& val, double tol) {
        return val<=tol && val>=-tol;
      }
      static bool is_one(const T& val) {
        return val==1;
      }
      static bool is_minus_one(const T& val) {
        return val==-1;
      }
      static bool is_constant(const T& val) {
        return true;
      }
      static bool is_integer(const T& val) {
        return val==static_cast<casadi_int>(val);
      }
      static bool is_inf(const T& val) {
        return std::numeric_limits<T>::has_infinity && val==std::numeric_limits<T>::infinity();
      }
      static bool is_minus_inf(const T& val) {
        return std::numeric_limits<T>::has_infinity && val==-std::numeric_limits<T>::infinity();
      }
      static bool is_nan(const T& val) {
        return std::numeric_limits<T>::has_quiet_NaN && val!=val;
      }
      static const T zero;
      static const T one;
      static const T two;
      static const T minus_one;
  };

  template<class T>
  inline bool is_zero(const T& x) {
    return casadi_limits<T>::is_zero(x);
  }

  template<class T>
  const T casadi_limits<T>::zero = T(0);

  template<class T>
  const T casadi_limits<T>::one = 1;

  template<class T>
  const T casadi_limits<T>::two = 2;

  template<class T>
  const T casadi_limits<T>::minus_one = -1;

} // namespace casadi
#endif // CASADI_CASADI_LIMITS_HPP
