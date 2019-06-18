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


#ifndef CASADI_SUBMATRIX_HPP
#define CASADI_SUBMATRIX_HPP

namespace casadi {


  /** SubMatrix class for Matrix
      SubMatrix is the return type for operator() of the Matrix class, it allows access to the value as well as changing the parent object
      \author Joel Andersson
      \date 2011-2014
  */
  template<typename M, typename I, typename J>
  class SubMatrix : public M {
  private:
    /// A reference to the matrix that is allowed to be modified
    M& mat_;

    /// The element of the matrix that is allowed to be modified
    I i_;
    J j_;
  public:
    /// Constructor
    SubMatrix(M& mat, const I& i, const J& j) : mat_(mat), i_(i), j_(j) {
      mat.get(*this, false, i, j);
    }

    ///@{
    /// Methods that modify a part of the parent object (A(i, j) = ?, A(i, j) += ?, etc.)
    inline const M& operator=(const SubMatrix<M, I, J> &y) {
      mat_.set(y, false, i_, j_);
      return y;
    }

    inline const M& operator=(const M &y) {
      mat_.set(y, false, i_, j_);
      return y;
    }

    inline M operator+=(const M &y) {
      M s = *this+y;
      mat_.set(s, false, i_, j_);
      return s;
    }

    inline M operator-=(const M &y) {
      M s = *this-y;
      mat_.set(s, false, i_, j_);
      return s;
    }

    inline M operator*=(const M &y) {
      M s = *this*y;
      mat_.set(s, false, i_, j_);
      return s;
    }

    inline M operator/=(const M &y) {
      M s = *this/y;
      mat_.set(s, false, i_, j_);
      return s;
    }
    ///@}
  };

  /** SubIndex class for Matrix
      Same as the above class but for single argument return for operator()
      \author Joel Andersson
      \date 2011-2014
  */
  template<typename M, typename I>
  class SubIndex : public M {
  private:
    /// A reference to the matrix that is allowed to be modified
    M& mat_;

    /// The element of the matrix that is allowed to be modified
    I i_;
  public:
    /// Constructor
    SubIndex(M& mat, const I& i) : mat_(mat), i_(i) {
      mat.get(*this, false, i);
    }

    ///@{
    /// Methods that modify a part of the parent object (A(i) = ?, A(i) += ?, etc.)
    inline const M& operator=(const SubIndex<M, I> &y) {
      mat_.set(y, false, i_);
      return y;
    }

    inline const M& operator=(const M &y) {
      mat_.set(y, false, i_);
      return y;
    }

    inline M operator+=(const M &y) {
      M s = *this+y;
      mat_.set(s, false, i_);
      return s;
    }

    inline M operator-=(const M &y) {
      M s = *this-y;
      mat_.set(s, false, i_);
      return s;
    }

    inline M operator*=(const M &y) {
      M s = *this*y;
      mat_.set(s, false, i_);
      return s;
    }

    inline M operator/=(const M &y) {
      M s = *this/y;
      mat_.set(s, false, i_);
      return s;
    }
    ///@}
  };

} // namespace casadi


#endif // CASADI_SUBMATRIX_HPP
