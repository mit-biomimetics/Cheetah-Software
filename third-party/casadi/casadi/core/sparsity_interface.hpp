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


#ifndef CASADI_SPARSITY_INTERFACE_HPP
#define CASADI_SPARSITY_INTERFACE_HPP

#include "casadi_misc.hpp"

namespace casadi {
  /** \brief Empty Base
      This class is extended in SWIG.
   */
  struct CASADI_EXPORT SparsityInterfaceCommon {};

#ifndef SWIG
  /** \brief Sparsity interface class

      This is a common base class for GenericMatrix (i.e. MX and Matrix<>) and Sparsity, introducing a
      uniform syntax and implementing common functionality using the curiously recurring template pattern
      (CRTP) idiom.\n

      \author Joel Andersson
      \date 2014
  */
  template<typename MatType>
  class SparsityInterface : public SparsityInterfaceCommon {
  protected:
    // Helper functions
    inline const MatType& self() const { return static_cast<const MatType&>(*this); }
    inline MatType& self() { return static_cast<MatType&>(*this); }
  public:

    /// \cond CLUTTER
    static std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, const std::vector<casadi_int>& vert_offset,
                 const std::vector<casadi_int>& horz_offset);
    static std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, casadi_int vert_incr, casadi_int horz_incr);
    static MatType veccat(const std::vector< MatType >& x);
    static MatType vec(const MatType& x);
    static MatType repmat(const MatType& x, casadi_int n, casadi_int m=1);
    static std::vector<casadi_int> offset(const std::vector< MatType > &v, bool vert=true);
    static std::vector< MatType > diagsplit(const MatType& x,
                                            const std::vector<casadi_int>& output_offset);
    static std::vector< MatType > diagsplit(const MatType& x, casadi_int incr);
    static std::vector< MatType > diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2);
    static MatType mtimes(const std::vector<MatType> &args);
    static std::vector<MatType > horzsplit(const MatType& x, casadi_int incr);
    static std::vector<MatType > vertsplit(const MatType& x, casadi_int incr);
    static MatType repmat(const MatType &A, const std::pair<casadi_int, casadi_int>& rc) {
      return MatType::repmat(A, rc.first, rc.second);
    }
    /// \endcond

/** \defgroup expression_tools Expression tools
* Functions for manipulating DM, SX, MX or Sparsity
*
*/
/**
\ingroup expression_tools
@{
*/
    /** \brief Concatenate a list of matrices horizontally
     * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
     *
     *   horzcat(horzsplit(x, ...)) = x
     */
    inline friend MatType horzcat(const std::vector<MatType> &v) {
      return MatType::horzcat(v);
    }

    /** \brief Concatenate a list of matrices vertically
     * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
     *
     *   vertcat(vertsplit(x, ...)) = x
     */
    inline friend MatType vertcat(const std::vector<MatType> &v) {
      return MatType::vertcat(v);
    }

    /** \brief  split horizontally, retaining groups of columns
     * \param offset List of all start columns for each group
     *      the last column group will run to the end.
     *
     *   horzcat(horzsplit(x, ...)) = x
     */
    inline friend std::vector<MatType >
      horzsplit(const MatType &x, const std::vector<casadi_int>& offset) {
      return MatType::horzsplit(x, offset);
    }

    /** \brief  split horizontally, retaining fixed-sized groups of columns
     * \param incr Size of each group of columns
     *
     *   horzcat(horzsplit(x, ...)) = x
     */
    inline friend std::vector<MatType > horzsplit(const MatType& x, casadi_int incr=1) {
      return MatType::horzsplit(x, incr);
    }

    /** * \brief  split vertically, retaining groups of rows
     * \param output_offset List of all start rows for each group
     *      the last row group will run to the end.
     *
     *   vertcat(vertsplit(x, ...)) = x
     */
    inline friend std::vector<MatType >
      vertsplit(const MatType& x, const std::vector<casadi_int>& offset) {
      return MatType::vertsplit(x, offset);
    }

    /** \brief Helper function, get offsets corresponding to a vector of matrices
     */
    inline friend std::vector<casadi_int > offset(const std::vector<MatType> &v, bool vert=true) {
      return MatType::offset(v, vert);
    }

    /** \brief  split vertically, retaining fixed-sized groups of rows
     * \param incr Size of each group of rows
     *
     *   vertcat(vertsplit(x, ...)) = x

     \doctest
     print vertsplit(SX.sym("a",4))
     \doctestout
     [SX(a_0), SX(a_1), SX(a_2), SX(a_3)]
     \enddoctest

     \doctest
     print vertsplit(SX.sym("a",4),2)
     \doctestout
     [SX([a_0, a_1]), SX([a_2, a_3])]
     \enddoctest

     If the number of rows is not a multiple of \p incr,
     the last entry returned will have a size smaller than \p incr.

     \doctest
     print vertsplit(DM([0,1,2,3,4]),2)
     \doctestout
     [DM([0, 1]), DM([2, 3]), DM(4)]
     \enddoctest
     *
     */
    inline friend std::vector<MatType > vertsplit(const MatType &x, casadi_int incr=1) {
      return MatType::vertsplit(x, incr);
    }

    /** \brief Construct a matrix from a list of list of blocks.
     */
    inline friend MatType blockcat(const std::vector< std::vector<MatType > > &v) {
      return MatType::blockcat(v);
    }

    /** \brief Construct a matrix from 4 blocks
     */
    inline friend MatType
      blockcat(const MatType &A, const MatType &B, const MatType &C, const MatType &D) {
      return vertcat(horzcat(A, B), horzcat(C, D));
    }

    /** \brief  chop up into blocks
     * \param vert_offset Defines the boundaries of the block rows
     * \param horz_offset Defines the boundaries of the block columns
     *
     *   blockcat(blocksplit(x,..., ...)) = x
     */
    inline friend std::vector< std::vector< MatType > >
      blocksplit(const MatType& x,
                 const std::vector<casadi_int>& vert_offset,
                 const std::vector<casadi_int>& horz_offset) {
      return MatType::blocksplit(x, vert_offset, horz_offset);
    }

    /** \brief  chop up into blocks
     * \param vert_incr Defines the increment for block boundaries in row dimension
     * \param horz_incr Defines the increment for block boundaries in column dimension
     *
     *   blockcat(blocksplit(x,..., ...)) = x
     */
    inline friend std::vector< std::vector< MatType > >
      blocksplit(const MatType& x, casadi_int vert_incr=1, casadi_int horz_incr=1) {
      return MatType::blocksplit(x, vert_incr, horz_incr);
    }

    /** \brief Construct a matrix with given block on the diagonal
     */
    inline friend MatType diagcat(const std::vector<MatType> &A) {
      return MatType::diagcat(A);
    }

    /** \brief  split diagonally, retaining square matrices
     * \param output_offset1 List of all start locations (row) for each group
     *      the last matrix will run to the end.
     * \param output_offset2 List of all start locations (row) for each group
     *      the last matrix will run to the end.
     *
     *   diagcat(diagsplit(x, ...)) = x
     */
    friend std::vector< MatType >
      diagsplit(const MatType& x,
                const std::vector<casadi_int>& output_offset1,
                const std::vector<casadi_int>& output_offset2) {
      return MatType::diagsplit(x, output_offset1, output_offset2);
    }

    /** \brief  split diagonally, retaining square matrices
     * \param output_offset List of all start locations for each group
     *      the last matrix will run to the end.
     *
     *   diagcat(diagsplit(x, ...)) = x
     */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, const std::vector<casadi_int>& output_offset) {
      return MatType::diagsplit(x, output_offset);
    }

    /** \brief  split diagonally, retaining groups of square matrices
     * \param incr Size of each matrix
     *
     *  diagsplit(diagsplit(x, ...)) = x
     */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, casadi_int incr=1) {
      return MatType::diagsplit(x, incr);
    }

    /** \brief  split diagonally, retaining fixed-sized matrices
     * \param incr1 Row dimension of each matrix
     * \param incr2 Column dimension of each matrix
     *
     *  diagsplit(diagsplit(x, ...)) = x
     */
    inline friend std::vector< MatType >
      diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2) {
      return MatType::diagsplit(x, incr1, incr2);
    }

    /** \brief  concatenate vertically while vectorizing all arguments with vec
     */
    inline friend MatType veccat(const std::vector< MatType >& x) {
      return MatType::veccat(x);
    }

    /** \brief Matrix product of two matrices
    */
    inline friend MatType mtimes(const MatType &x, const MatType &y) {
      return MatType::mtimes(x, y);
    }

    /** \brief Matrix product of n matrices
     */
    inline friend MatType mtimes(const std::vector<MatType> &args) {
      return MatType::mtimes(args);
    }

    /** \brief Multiply-accumulate operation
        Matrix product of two matrices (x and y), adding the result to
        a third matrix z. The result has the same sparsity pattern as
        C meaning that other entries of (x*y) are ignored.
        The operation is equivalent to: z+mtimes(x,y).project(z.sparsity()).
    */
    inline friend MatType
      mac(const MatType &x, const MatType &y, const MatType &z) {
      return MatType::mac(x, y, z);
    }

    /** \brief Transpose
     */
    inline friend MatType transpose(const MatType& X) {
      return X.T();
    }

    /** \brief  make a vector
        Reshapes/vectorizes the matrix such that the shape becomes (expr.numel(), 1).
        Columns are stacked on top of each other.
        Same as reshape(expr, expr.numel(), 1)

        a c \n
        b d \n

        turns into

        a \n
        b \n
        c \n
        d \n
    */
    inline friend MatType vec(const MatType& x) {
      return MatType::vec(x);
    }

    /** \brief Returns a reshaped version of the matrix
     */
    inline friend MatType reshape(const MatType& x, casadi_int nrow, casadi_int ncol) {
      return MatType::reshape(x, nrow, ncol);
    }

    /** \brief Returns a reshaped version of the matrix, dimensions as a vector
    */
    inline friend MatType reshape(const MatType& x, std::pair<casadi_int, casadi_int> rc) {
      return MatType::reshape(x, rc.first, rc.second);
    }

    /** \brief Reshape the matrix
    */
    inline friend MatType reshape(const MatType& x, const Sparsity& sp) {
      return MatType::reshape(x, sp);
    }

    /** \brief Obtain the structural rank of a sparsity-pattern
    */
    inline friend casadi_int sprank(const MatType& x) {
      return MatType::sprank(x);
    }

    /** \brief 0-norm (nonzero count) of a Matrix-matrix product
    */
    inline friend casadi_int norm_0_mul(const MatType &x, const MatType &y) {
      return MatType::norm_0_mul(x, y);
    }

    /** \brief Get the upper triangular part of a matrix
    */
    inline friend MatType triu(const MatType& x, bool includeDiagonal=true) {
      return MatType::triu(x, includeDiagonal);
    }

    /** \brief Get the lower triangular part of a matrix
    */
    inline friend MatType tril(const MatType& x, bool includeDiagonal=true) {
      return MatType::tril(x, includeDiagonal);
    }

    /** \brief Kronecker tensor product
     *
     * Creates a block matrix in which each element (i, j) is a_ij*b
     */
    inline friend MatType kron(const MatType& a, const MatType& b) {
      return MatType::kron(a, b);
    }

    /** \brief Repeat matrix A n times vertically and m times horizontally
     */
    inline friend MatType repmat(const MatType &A, casadi_int n, casadi_int m=1) {
      return MatType::repmat(A, n, m);
    }

    /** \brief Repeat matrix A n times vertically and m times horizontally
     */
    inline friend MatType repmat(const MatType &A, const std::pair<casadi_int, casadi_int>& rc) {
      return MatType::repmat(A, rc);
    }

    /** \brief Concatenate horizontally, two matrices */
    inline friend MatType horzcat(const MatType &x, const MatType &y) {
      return horzcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate horizontally, three matrices */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z) {
      return horzcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate horizontally, four matrices */
    inline friend MatType horzcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return horzcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Concatenate vertically, two matrices */
    inline friend MatType vertcat(const MatType &x, const MatType &y) {
      return vertcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate vertically, three matrices */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z) {
      return vertcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate vertically, four matrices */
    inline friend MatType vertcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return vertcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Concatenate along diagonal, two matrices */
    inline friend MatType diagcat(const MatType &x, const MatType &y) {
      return diagcat(std::vector<MatType>{x, y});
    }

    /** \brief Concatenate along diagonal, three matrices */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z) {
      return diagcat(std::vector<MatType>{x, y, z});
    }

    /** \brief Concatenate along diagonal, four matrices */
    inline friend MatType diagcat(const MatType &x, const MatType &y, const MatType &z,
                                  const MatType &w) {
      return diagcat(std::vector<MatType>{x, y, z, w});
    }

    /** \brief Return a row-wise summation of elements */
    inline friend MatType sum1(const MatType &x) { return MatType::sum1(x);}

    /** \brief Return a column-wise summation of elements  */
    inline friend MatType sum2(const MatType &x) { return MatType::sum2(x);}

/** \@} */
  };
#endif // SWIG

#ifndef SWIG
  template<typename MatType>
  MatType SparsityInterface<MatType>::vec(const MatType& x) {
    if (x.size2()==1) return x;
    return reshape(x, x.numel(), 1);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::repmat(const MatType& x, casadi_int n, casadi_int m) {
    if (n==1 && m==1) return x;
    MatType allrows = vertcat(std::vector<MatType>(n, x));
    if (n==0) allrows = MatType(0, x.size2());
    MatType ret = horzcat(std::vector<MatType>(m, allrows));
    if (m==0) ret = MatType(allrows.size1(), 0);
    return ret;
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::blocksplit(const MatType& x,
                                         const std::vector<casadi_int>& vert_offset,
                                         const std::vector<casadi_int>& horz_offset) {
    std::vector<MatType> rows = MatType::vertsplit(x, vert_offset);
    std::vector< std::vector< MatType > > ret;
    for (auto&& r : rows) ret.push_back(MatType::horzsplit(r, horz_offset));
    return ret;
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::blocksplit(const MatType& x,
      casadi_int vert_incr, casadi_int horz_incr) {
    casadi_assert_dev(horz_incr>=1);
    casadi_assert_dev(vert_incr>=1);
    casadi_int sz1 = x.size1();
    std::vector<casadi_int> offset1 = range(0, sz1, vert_incr);
    offset1.push_back(sz1);
    casadi_int sz2 = x.size2();
    std::vector<casadi_int> offset2 = range(0, sz2, horz_incr);
    offset2.push_back(sz2);
    return blocksplit(x, offset1, offset2);
  }

  template<typename MatType>
  std::vector<casadi_int>
  SparsityInterface<MatType>::offset(const std::vector< MatType > &v, bool vert) {
    std::vector<casadi_int> ret(v.size()+1);
    ret[0]=0;
    for (casadi_int i=0; i<v.size(); ++i) {
      ret[i+1] = ret[i] + (vert ? v[i].size1() : v[i].size2());
    }
    return ret;
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::veccat(const std::vector< MatType >& x) {
    std::vector< MatType > x_vec = x;
    for (typename std::vector< MatType >::iterator it=x_vec.begin();
         it!=x_vec.end(); ++it) {
      *it = vec(*it);
    }
    return vertcat(x_vec);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x,
      const std::vector<casadi_int>& output_offset) {
    casadi_assert(x.is_square(), "diagsplit(x,incr)::input must be square but got "
                          + x.dim()  + ".");
    return MatType::diagsplit(x, output_offset, output_offset);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x, casadi_int incr) {
    casadi_assert_dev(incr>=1);
    casadi_assert(x.is_square(), "diagsplit(x,incr)::input must be square but got "
                          + x.dim()  + ".");
    std::vector<casadi_int> offset2 = range(0, x.size2(), incr);
    offset2.push_back(x.size2());
    return MatType::diagsplit(x, offset2);
  }

  template<typename MatType>
  std::vector< MatType >
  SparsityInterface<MatType>::diagsplit(const MatType& x, casadi_int incr1, casadi_int incr2) {
    casadi_assert_dev(incr1>=1);
    casadi_assert_dev(incr2>=1);
    std::vector<casadi_int> offset1 = range(0, x.size1(), incr1);
    offset1.push_back(x.size1());
    std::vector<casadi_int> offset2 = range(0, x.size2(), incr2);
    offset2.push_back(x.size2());
    return MatType::diagsplit(x, offset1, offset2);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::mtimes(const std::vector<MatType> &args) {
    casadi_assert(args.size()>=1,
                          "mul(std::vector<MatType> &args): "
                          "supplied list must not be empty.");
    MatType ret = args[0];
    for (casadi_int i=1; i<args.size(); ++i) ret = MatType::mtimes(ret, args[i]);
    return ret;
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::horzsplit(const MatType& x, casadi_int incr) {
    casadi_assert_dev(incr>=1);
    casadi_int sz2 = x.size2();
    std::vector<casadi_int> offset2 = range(0, sz2, incr);
    offset2.push_back(sz2);
    return MatType::horzsplit(x, offset2);
  }

  template<typename MatType>
  std::vector<MatType > SparsityInterface<MatType>::vertsplit(const MatType& x, casadi_int incr) {
    casadi_assert_dev(incr>=1);
    casadi_int sz1 = x.size1();
    std::vector<casadi_int> offset1 = range(0, sz1, incr);
    offset1.push_back(sz1);
    return MatType::vertsplit(x, offset1);
  }
#endif // SWIG

} // namespace casadi

#endif // CASADI_SPARSITY_INTERFACE_HPP
