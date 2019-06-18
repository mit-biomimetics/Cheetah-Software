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


#ifndef CASADI_NONZEROS_HPP
#define CASADI_NONZEROS_HPP

namespace casadi {


/** NonZeros class for Matrix
  NonZeros is the return type for operator[] of the Matrix class, it allows access to the value as well as changing the parent object
  \author Joel Andersson
  \date 2011
*/

/// Access to a set of nonzeros
template<typename M, typename K>
class NonZeros : public M {
  public:
    /// Constructor
    NonZeros(M& mat, const K& k) : mat_(mat), k_(k) { mat.get_nz(*this, false, k); }

    ///@{
    /// Methods that modify a part of the parent object (A[k] = ?, A[k] += ?, etc.)
    const M& operator=(const NonZeros<M, K> &y);
    const M& operator=(const M &y);
    M operator+=(const M &y);
    M operator-=(const M &y);
    M operator*=(const M &y);
    M operator/=(const M &y);
    ///@}

  private:
    /// A reference to the matrix that is allowed to be modified
    M& mat_;

    /// The element of the matrix that is allowed to be modified
    K k_;
};

// Implementation
template<typename M, typename K>
const M& NonZeros<M, K>::operator=(const NonZeros<M, K> &y) {
  mat_.set_nz(y, false, k_);
  return y;
}

// Implementation
template<typename M, typename K>
const M& NonZeros<M, K>::operator=(const M &y) {
  mat_.set_nz(y, false, k_);
  return y;
}

template<typename M, typename K>
M NonZeros<M, K>::operator+=(const M &y) {
  M s = *this+y;
  mat_.set_nz(s, false, k_);
  return s;
}

template<typename M, typename K>
M NonZeros<M, K>::operator-=(const M &y) {
  M s = *this-y;
  mat_.set_nz(s, false, k_);
  return s;
}

template<typename M, typename K>
M NonZeros<M, K>::operator*=(const M &y) {
   M s = *this*y;
   mat_.set_nz(s, false, k_);
   return s;
}

template<typename M, typename K>
M NonZeros<M, K>::operator/=(const M &y) {
  M s = *this/y;
  mat_.set_nz(s, false, k_);
  return s;
}

} // namespace casadi


#endif // CASADI_NONZEROS_HPP
