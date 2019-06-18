/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *    Copyright (C) 2005-2013 Timothy A. Davis
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


#ifndef CASADI_SPARSITY_HPP
#define CASADI_SPARSITY_HPP

#include "shared_object.hpp"
#include "printable.hpp"
#include "casadi_common.hpp"
#include "sparsity_interface.hpp"
#include "generic_type.hpp"
#include <vector>
#include <list>
#include <limits>
#include <unordered_map>

namespace casadi {
  // Forward declaration
  class SparsityInternal;

  #ifndef SWIG
    /** \brief Compact representation of a sparsity pattern */
    struct CASADI_EXPORT SparsityStruct {
      casadi_int nrow;
      casadi_int ncol;
      const casadi_int* colind;
      const casadi_int* row;
    };
  #endif // SWIG

  /** \brief General sparsity class
   *
   * The storage format is a compressed column storage (CCS) format.\n
   *
   In this format, the structural non-zero elements are stored in column-major order, starting from
   the upper left corner of the matrix and ending in the lower right corner.

   In addition to the dimension (size1(), size2()), (i.e. the number of rows and the number of
   columns respectively), there are also two vectors of integers:

   1. "colind" [length size2()+1], which contains the index to the first non-zero element on or
   after the corresponding column. All the non-zero elements of a particular i are thus the elements
   with index el that fulfills: colind[i] <= el < colind[i+1].

   2. "row" [same length as the number of non-zero elements, nnz()] The rows for each of the
   structural non-zeros.

   Note that with this format, it is cheap to loop over all the non-zero elements of a particular
   column, at constant time per element, but expensive to jump to access a location (i, j).

   If the matrix is dense, i.e. length(row) == size1()*size2(), the format reduces to standard dense
   column major format, which allows access to an arbitrary element in constant time.

   Since the object is reference counted (it inherits from SharedObject), several matrices are
   allowed to share the same sparsity pattern.

   The implementations of methods marked as such in this class has been taken from the
   CSparse package and modified to fit CasADi data structures and separation of
   sparsity pattern calculation and numerical evaluation.
   These functions are Copyright(c) Timothy A. Davis, 2006-2009
   and licensed as a derivative work under the GNU LGPL

   * \see Matrix
   *
   * \author Joel Andersson
   * \date 2010-2015
   */
  class CASADI_EXPORT Sparsity
    : public SharedObject,
      public SWIG_IF_ELSE(SparsityInterfaceCommon, SparsityInterface<Sparsity>),
      public SWIG_IF_ELSE(PrintableCommon, Printable<Sparsity>) {
  public:

    /// Default constructor
    explicit Sparsity(casadi_int dummy=0);

    /** \brief Pattern with all structural zeros */
    Sparsity(casadi_int nrow, casadi_int ncol);

    /// Construct from sparsity pattern vectors given in compressed column storage format
    Sparsity(casadi_int nrow, casadi_int ncol,
             const std::vector<casadi_int>& colind, const std::vector<casadi_int>& row,
             bool order_rows=false);

    /** \brief Create a sparse matrix with all structural zeros */
    explicit Sparsity(const std::pair<casadi_int, casadi_int>& rc);

#ifndef SWIG
    /// Construct from sparsity pattern vectors given in compressed column storage format
    Sparsity(casadi_int nrow, casadi_int ncol, const casadi_int* colind, const casadi_int* row,
            bool order_rows=false);

    /** \brief  Create from node */
    static Sparsity create(SparsityInternal *node);

    /// Base class
    typedef SparsityInterface<Sparsity> B;

    /// Expose base class functions
    using B::horzsplit;
    using B::diagsplit;
    using B::vertsplit;
    using B::mtimes;
#endif

    /** \brief Create a scalar sparsity pattern **/
    ///@{
    static Sparsity scalar(bool dense_scalar=true)
    { return dense_scalar ? dense(1, 1) : Sparsity(1, 1); }
    ///@}

    /** \brief Create a dense rectangular sparsity pattern **/
    ///@{
    static Sparsity dense(casadi_int nrow, casadi_int ncol=1);
    static Sparsity dense(const std::pair<casadi_int, casadi_int> &rc) {
      return dense(rc.first, rc.second);
    }
    ///@}

    /** \brief Create the sparsity pattern for a unit vector of length n and a nonzero on
     * position el **/
    ///@{
    static Sparsity unit(casadi_int n, casadi_int el);
    ///@}

    /** \brief Create a upper triangular square sparsity pattern **/
    static Sparsity upper(casadi_int n);

    /** \brief Create a lower triangular square sparsity pattern **/
    static Sparsity lower(casadi_int n);

    /** \brief Create diagonal sparsity pattern **/
    ///@{
    static Sparsity diag(casadi_int nrow) { return diag(nrow, nrow);}
    static Sparsity diag(casadi_int nrow, casadi_int ncol);
    static Sparsity diag(const std::pair<casadi_int, casadi_int> &rc) {
      return diag(rc.first, rc.second);
    }
    ///@}

    /** \brief Create a single band in a square sparsity pattern
     *
     * band(n, 0) is equivalent to diag(n) \n
     * band(n, -1) has a band below the diagonal \n
     * \param p indicate
     **/
    static Sparsity band(casadi_int n, casadi_int p);

    /** \brief Create banded square sparsity pattern
     *
     * banded(n, 0) is equivalent to diag(n) \n
     * banded(n, 1) is tri-diagonal matrix \n
     **/
    static Sparsity banded(casadi_int n, casadi_int p);

    /** \brief Construct a block sparsity pattern from (row, col) vectors */
    static Sparsity rowcol(const std::vector<casadi_int>& row,
                           const std::vector<casadi_int>& col,
                           casadi_int nrow, casadi_int ncol);

    /** \brief Create a sparsity pattern given the nonzeros in sparse triplet form
    **/
    static Sparsity triplet(casadi_int nrow, casadi_int ncol,
                            const std::vector<casadi_int>& row, const std::vector<casadi_int>& col,
                            std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool invert_mapping);

    /** \brief Create a sparsity pattern given the nonzeros in sparse triplet form
        (no nonzero mapping)
        rows_are_sorted==true means that the row entries already in increasing order
        for each col and without any duplicates
    **/
    static Sparsity triplet(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& row,
                            const std::vector<casadi_int>& col);

    /** \brrief Create a sparsity from nonzeros
    *
    * Inverse of `find()`
    */
    static Sparsity nonzeros(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& nz,
                              bool ind1=SWIG_IND1);

    /** Create from a single vector containing the pattern in compressed column storage format:
     * The format:
     * The first two entries are the number of rows (nrow) and columns (ncol)
     * The next ncol+1 entries are the column offsets (colind).
       Note that the last element, colind[ncol], gives the number of nonzeros
     * The last colind[ncol] entries are the row indices
     **/
    ///@{
    static Sparsity compressed(const std::vector<casadi_int>& v, bool order_rows=false);
#ifndef SWIG
    static Sparsity compressed(const casadi_int* v, bool order_rows=false);
#endif // SWIG
    ///@}

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Correctness of sparsity patterns are checked during
              construction */
    void sanity_check(bool complete=false) const {}
#endif // WITH_DEPRECATED_FEATURES

    /** Get the diagonal of the matrix/create a diagonal matrix
        (mapping will contain the nonzero mapping)
        When the input is square, the diagonal elements are returned.
        If the input is vector-like, a diagonal matrix is constructed with it.
    */
    Sparsity get_diag(std::vector<casadi_int>& SWIG_OUTPUT(mapping)) const;

    /// Compress a sparsity pattern
    std::vector<casadi_int> compress() const;

#ifndef SWIG
    /// Access a member function or object
    const SparsityInternal* operator->() const;

    /// Reference to internal structure
    const SparsityInternal& operator*() const;
#endif // SWIG
    /// \name Check if two sparsity patterns are identical
    /// @{
    bool is_equal(const Sparsity& y) const;
    bool is_equal(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                 const std::vector<casadi_int>& row) const;
#ifndef SWIG
    bool is_equal(casadi_int nrow, casadi_int ncol,
                  const casadi_int* colind, const casadi_int* row) const;
#endif // SWIG

    bool operator==(const Sparsity& y) const { return is_equal(y);}
    /// @}

    /// Check if two sparsity patterns are difference
    bool operator!=(const Sparsity& y) const {return !is_equal(y);}

    /// Check if pattern is horizontal repeat of another
    bool is_stacked(const Sparsity& y, casadi_int n) const;

#ifndef SWIG
    /** \brief Implicit or explicit type conversion to C representation
        In the C runtime, sparsity patterns are represented as a "const casadi_int*".
        This enables using the C runtime functions using a natural syntax.
    */
    operator const casadi_int*() const;

    /** \brief Implicit or explicit type conversion to compact representation */
    operator const std::vector<casadi_int>&() const;

    /** \brief Implicit or explicit type conversion to C representation */
    operator SparsityStruct() const;
#endif // SWIG

    /// \name Size and element counting
    /// @{

    /// Get the number of rows
    casadi_int size1() const;

    /// Get the number of rows, Octave-style syntax
    casadi_int rows() const {return size1();}

    /// Get the number of columns
    casadi_int size2() const;

    /// Get the number of columns, Octave-style syntax
    casadi_int columns() const {return size2();}

    /** \brief The total number of elements, including structural zeros, i.e. size2()*size1()
      * Beware of overflow
      * \see nnz()
      */
    casadi_int numel() const;

    /** \brief The percentage of nonzero
     * Equivalent to (100.0 * nnz())/numel(), but avoids overflow
     */
    double density() const;

    /** \brief Check if the sparsity is empty
     *
     * A sparsity is considered empty if one of the dimensions is zero
     * (or optionally both dimensions)
     */
    bool is_empty(bool both=false) const;

    /** \brief Get the number of (structural) non-zeros
      * \see numel()
      */
    casadi_int nnz() const;

    /** \brief Number of non-zeros in the upper triangular half,
     * i.e. the number of elements (i, j) with j>=i */
    casadi_int nnz_upper(bool strictly=false) const;

    /** \brief Number of non-zeros in the lower triangular half,
     * i.e. the number of elements (i, j) with j<=i */
    casadi_int nnz_lower(bool strictly=false) const;

    /** \brief Number of non-zeros on the diagonal, i.e. the number of elements (i, j) with j==i */
    casadi_int nnz_diag() const;

    /** \brief Upper half-bandwidth */
    casadi_int bw_upper() const;

    /** \brief Lower half-bandwidth */
    casadi_int bw_lower() const;

    /** \brief  Get the shape */
    std::pair<casadi_int, casadi_int> size() const;

    /** \brief  Get the size along a particular dimensions */
    casadi_int size(casadi_int axis) const;
    /// @}

    /** Obtain information about sparsity */
    Dict info() const;

    /** Construct instance from info */
    static Sparsity from_info(const Dict& info);

    /** Export sparsity pattern to file
    *
    * Supported formats:
    *   - .mtx   Matrix Market
    */
    void to_file(const std::string& filename, const std::string& format_hint="") const;

    static Sparsity from_file(const std::string& filename, const std::string& format_hint="");

#ifndef SWIG
    /** \brief Serialize */
    void serialize(std::ostream &stream) const;
#endif

    /** \brief Serialize */
    std::string serialize() const;

    /** \brief Build Sparsity from serialization */
    static Sparsity deserialize(std::istream& istream);

    /** \brief Build Sparsity from serialization */
    static Sparsity deserialize(const std::string& s);

#ifndef SWIG
    /** \brief Get a reference to row-vector,
     * containing rows for all non-zero elements (see class description) */
    const casadi_int* row() const;

    /** \brief Get a reference to the colindex of all column element (see class description) */
    const casadi_int* colind() const;
#endif

    /** \brief Get the row for each non-zero entry
        Together with the column-vector, this vector gives the sparsity of the matrix in
        sparse triplet format, and together with the colind vector, one obtains the sparsity
        in column compressed format. */
    std::vector<casadi_int> get_row() const;

    /** \brief Get the column index for each column
        Together with the row-vector, one obtains the sparsity pattern in the
        column compressed format. */
    std::vector<casadi_int> get_colind() const;

    /** \brief  Get a reference to the colindex of column cc (see class description) */
    casadi_int colind(casadi_int cc) const;

    /** \brief Get the row of a non-zero element */
    casadi_int row(casadi_int el) const;

    /** \brief Get the column for each non-zero entry
        Together with the row-vector, this vector gives the sparsity of the matrix in
        sparse triplet format, i.e. the column and row for each non-zero elements  */
    std::vector<casadi_int> get_col() const;

    /// Resize
    void resize(casadi_int nrow, casadi_int ncol);

    /** \brief Get the index of a non-zero element
        Add the element if it does not exist and copy object if it's not unique */
    casadi_int add_nz(casadi_int rr, casadi_int cc);

    /** \brief Get the index of an existing non-zero element
        return -1 if the element does not exist */
    casadi_int get_nz(casadi_int rr, casadi_int cc) const;

    /// Returns true if the pattern has a non-zero at location rr, cc
    bool has_nz(casadi_int rr, casadi_int cc) const;

    /** \brief Get a set of non-zero element
        return -1 if the element does not exist */
    std::vector<casadi_int> get_nz(const std::vector<casadi_int>& rr,
                                    const std::vector<casadi_int>& cc) const;

    /** \brief Get the nonzero index for a set of elements
        The index vector is used both for input and outputs and must be sorted by increasing
        nonzero index, i.e. column-wise.
        Elements not found in the sparsity pattern are set to -1.
    */
    void get_nz(std::vector<casadi_int>& SWIG_INOUT(indices)) const;

    /// Get nonzeros in lower triangular part
    std::vector<casadi_int> get_lower() const;

    /// Get nonzeros in upper triangular part
    std::vector<casadi_int> get_upper() const;

    /// Get the sparsity in compressed column storage (CCS) format
    void get_ccs(std::vector<casadi_int>& SWIG_OUTPUT(colind),
                std::vector<casadi_int>& SWIG_OUTPUT(row)) const;

    /// Get the sparsity in compressed row storage (CRS) format
    void get_crs(std::vector<casadi_int>& SWIG_OUTPUT(rowind),
                std::vector<casadi_int>& SWIG_OUTPUT(col)) const;

    /// Get the sparsity in sparse triplet format
    void get_triplet(std::vector<casadi_int>& SWIG_OUTPUT(row),
                    std::vector<casadi_int>& SWIG_OUTPUT(col)) const;

    /** \brief Get a submatrix
     *
     * Returns the sparsity of the submatrix, with a mapping such that
     *   submatrix[k] = originalmatrix[mapping[k]]
     */
    Sparsity sub(const std::vector<casadi_int>& rr,
                 const std::vector<casadi_int>& cc,
                 std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool ind1=false) const;

    /** \brief Get a set of elements
     *
     * Returns the sparsity of the corresponding elements, with a mapping such that
     *   submatrix[k] = originalmatrix[mapping[k]]
     */
    Sparsity sub(const std::vector<casadi_int>& rr, const Sparsity& sp,
                 std::vector<casadi_int>& SWIG_OUTPUT(mapping), bool ind1=false) const;

    /// Transpose the matrix
    Sparsity T() const;

    /** \brief Transpose the matrix and get the reordering of the non-zero entries
    *
    *  \param[out] mapping the non-zeros of the original matrix
    *              for each non-zero of the new matrix
    */
    Sparsity transpose(std::vector<casadi_int>& SWIG_OUTPUT(mapping),
                        bool invert_mapping=false) const;

    /// Check if the sparsity is the transpose of another
    bool is_transpose(const Sparsity& y) const;

    /// Check if the sparsity is a reshape of another
    bool is_reshape(const Sparsity& y) const;

    /// @{
    /** \brief Combine two sparsity patterns
        Returns the new sparsity pattern as well as a mapping with the same length as
        the number of non-zero elements
        The mapping matrix contains the arguments for each nonzero, the first bit indicates
        if the first argument is nonzero,
        the second bit indicates if the second argument is nonzero (note that none of,
        one of or both of the arguments can be nonzero) */
#ifndef SWIG
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool fx0_is_zero,
                            std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool fx0_is_zero) const;
    /// @}

    /// @{
    /** \brief Union of two sparsity patterns */
#ifndef SWIG
    Sparsity unite(const Sparsity& y, std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity unite(const Sparsity& y) const;
    Sparsity operator+(const Sparsity& b) const;
    /// @}

    /// @{
    /** \brief Intersection of two sparsity patterns
        Returns the new sparsity pattern as well as a mapping with the same length as the
        number of non-zero elements
        The value is 1 if the non-zero comes from the first (i.e. this) object, 2 if it is from
        the second and 3 (i.e. 1 | 2) if from both */
#ifndef SWIG
    Sparsity intersect(const Sparsity& y,
                       std::vector<unsigned char>& mapping) const;
#endif // SWIG
    Sparsity intersect(const Sparsity& y) const;
    Sparsity operator*(const Sparsity& b) const;
    /// @}

    /// Take the inverse of a sparsity pattern; flip zeros and non-zeros
    Sparsity pattern_inverse() const;

#ifndef SWIG
    /** \brief Propagate sparsity using 0-1 logic through a matrix product,
     * no memory allocation: <tt>z = mul(x, y)</tt> with work vector
     * Forward mode.
     */
    static void mul_sparsityF(const bvec_t* x, const Sparsity& x_sp,
                              const bvec_t* y, const Sparsity& y_sp,
                              bvec_t* z, const Sparsity& z_sp,
                              bvec_t* w);

    /** \brief Propagate sparsity using 0-1 logic through a matrix product,
     * no memory allocation: <tt>z = mul(x, y)</tt> with work vector
     * Reverse mode.
     */
    static void mul_sparsityR(bvec_t* x, const Sparsity& x_sp,
                              bvec_t* y, const Sparsity& y_sp,
                              bvec_t* z, const Sparsity& z_sp,
                              bvec_t* w);

    /// \cond INTERNAL
    /// @{
    /** \brief Accessed by SparsityInterface */
    static Sparsity horzcat(const std::vector<Sparsity> & sp);
    static Sparsity vertcat(const std::vector<Sparsity> & sp);
    static Sparsity blockcat(const std::vector< std::vector< Sparsity > > &v);
    static Sparsity diagcat(const std::vector< Sparsity > &v);
    static std::vector<Sparsity>
      horzsplit(const Sparsity& x, const std::vector<casadi_int>& output_offset);
    static std::vector<Sparsity>
      vertsplit(const Sparsity& x, const std::vector<casadi_int>& output_offset);
    static std::vector<Sparsity>
      diagsplit(const Sparsity& x,
                const std::vector<casadi_int>& offset1,
                const std::vector<casadi_int>& offset2);
    static Sparsity mtimes(const Sparsity& x, const Sparsity& y);
    static Sparsity mac(const Sparsity& x, const Sparsity& y, const Sparsity& z) { return z;}
    static Sparsity reshape(const Sparsity& x, casadi_int nrow, casadi_int ncol);
    static Sparsity reshape(const Sparsity& x, const Sparsity& sp);
    static casadi_int sprank(const Sparsity& x);
    static casadi_int norm_0_mul(const Sparsity& x, const Sparsity& B);
    static Sparsity kron(const Sparsity& x, const Sparsity& b);
    static Sparsity triu(const Sparsity& x, bool includeDiagonal=true);
    static Sparsity tril(const Sparsity& x, bool includeDiagonal=true);

    static Sparsity sum2(const Sparsity &x);
    static Sparsity sum1(const Sparsity &x);
#endif //SWIG

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros

        For the matrices A to B
        A(m, n)
        length(jj)=m , length(ii)=n
        B(nrow, ncol)

        A=enlarge(m, n, ii, jj) makes sure that

        B[jj, ii] == A
    */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr,
                  const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Enlarge the matrix along the first dimension (i.e. insert rows) */
    void enlargeRows(casadi_int nrow, const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Enlarge the matrix along the second dimension (i.e. insert columns) */
    void enlargeColumns(casadi_int ncol, const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Make a patten dense */
    Sparsity makeDense(std::vector<casadi_int>& SWIG_OUTPUT(mapping)) const;

    /** \brief Erase rows and/or columns of a matrix */
    std::vector<casadi_int> erase(const std::vector<casadi_int>& rr,
                                  const std::vector<casadi_int>& cc, bool ind1=false);

    /** \brief Erase elements of a matrix */
    std::vector<casadi_int> erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /// Append another sparsity patten vertically (NOTE: only efficient if vector)
    void append(const Sparsity& sp);

    /// Append another sparsity patten horizontally
    void appendColumns(const Sparsity& sp);

    /// Is scalar?
    bool is_scalar(bool scalar_and_dense=false) const;

    /// Is dense?
    bool is_dense() const;

    /** \brief  Check if the pattern is a row vector (i.e. size1()==1) */
    bool is_row() const;

    /** \brief  Check if the pattern is a column vector (i.e. size2()==1) */
    bool is_column() const;

    /** \brief  Check if the pattern is a row or column vector */
    bool is_vector() const;

    /// Is diagonal?
    bool is_diag() const;

    /// Is square?
    bool is_square() const;

    /// Is symmetric?
    bool is_symmetric() const;

    /// Is upper triangular?
    bool is_triu() const;

    /// Is lower triangular?
    bool is_tril() const;

    /// Check whether the sparsity-pattern indicates structural singularity
    bool is_singular() const;

    /** \brief Do the rows appear sequentially on each column
    *
    * \param[in] strictly if true, then do not allow multiple entries
    */
    bool rowsSequential(bool strictly=true) const;

    /** \brief Remove duplicate entries
    *
    * The same indices will be removed from the \a mapping vector,
    * which must have the same length as the number of nonzeros
    */
    void removeDuplicates(std::vector<casadi_int>& SWIG_INOUT(mapping));

#ifndef SWIG
    typedef std::unordered_multimap<std::size_t, WeakRef> CachingMap;

    /// Cached sparsity patterns
    static CachingMap& getCache();

    /// (Dense) scalar
    static const Sparsity& getScalar();

    /// (Sparse) scalar
    static const Sparsity& getScalarSparse();

    /// Empty zero-by-zero
    static const Sparsity& getEmpty();

#endif //SWIG

    /** \brief Calculate the elimination tree
        See Direct Methods for Sparse Linear Systems by Davis (2006).
        If the parameter ata is false, the algorithm is equivalent to MATLAB's etree(A), except that
        the indices are zero-based. If ata is true, the algorithm is equivalent to MATLAB's
        etree(A, 'col').

        The implementation is a modified version of cs_etree in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL
    */
    std::vector<casadi_int> etree(bool ata=false) const;

    /** \brief Symbolic LDL factorization
        Returns the sparsity pattern of L^T

        The implementation is a modified version of LDL
        Copyright(c) Timothy A. Davis, 2005-2013
        Licensed as a derivative work under the GNU LGPL
    */
    Sparsity ldl(std::vector<casadi_int>& SWIG_OUTPUT(p), bool amd=true) const;

    /** \brief Symbolic QR factorization
        Returns the sparsity pattern of V (compact representation of Q) and R
        as well as vectors needed for the numerical factorization and solution.
        The implementation is a modified version of CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL
    */
    void qr_sparse(Sparsity& SWIG_OUTPUT(V), Sparsity& SWIG_OUTPUT(R),
                   std::vector<casadi_int>& SWIG_OUTPUT(prinv),
                   std::vector<casadi_int>& SWIG_OUTPUT(pc), bool amd=true) const;

    /** \brief Depth-first search on the adjacency graph of the sparsity
        See Direct Methods for Sparse Linear Systems by Davis (2006).
    */
    casadi_int dfs(casadi_int j, casadi_int top, std::vector<casadi_int>& SWIG_INOUT(xi),
            std::vector<casadi_int>& SWIG_INOUT(pstack),
            const std::vector<casadi_int>& pinv, std::vector<bool>& SWIG_INOUT(marked)) const;

    /** \brief Find the strongly connected components of the bigraph defined by the sparsity pattern
        of a square matrix

        See Direct Methods for Sparse Linear Systems by Davis (2006).
        Returns:
        - Number of components
        - Offset for each components (length: 1 + number of components)
        - Indices for each components, component i has indices
          index[offset[i]], ..., index[offset[i+1]]

        In the case that the matrix is symmetric, the result has a particular interpretation:
        Given a symmetric matrix A and
        n = A.scc(p, r)

        => A[p, p] will appear block-diagonal with n blocks and
        with the indices of the block boundaries to be found in r.

        The implementation is a modified version of cs_scc in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL
    */
    casadi_int scc(std::vector<casadi_int>& SWIG_OUTPUT(index),
            std::vector<casadi_int>& SWIG_OUTPUT(offset)) const;

    /** \brief Calculate the block triangular form (BTF)
        See Direct Methods for Sparse Linear Systems by Davis (2006).

        The function computes the Dulmage-Mendelsohn decomposition, which allows you to reorder
        the rows and columns of a matrix to bring it into block triangular form (BTF).

        It will not consider the distance of off-diagonal elements to the diagonal:
        there is no guarantee you will get a block-diagonal matrix if you supply a randomly
        permuted block-diagonal matrix.

        If your matrix is symmetrical, this method is of limited use; permutation can make it
        non-symmetric.

        \sa scc

        The implementation is a modified version of cs_dmperm in CSparse
        Copyright(c) Timothy A. Davis, 2006-2009
        Licensed as a derivative work under the GNU LGPL
    */
    casadi_int btf(std::vector<casadi_int>& SWIG_OUTPUT(rowperm),
            std::vector<casadi_int>& SWIG_OUTPUT(colperm),
            std::vector<casadi_int>& SWIG_OUTPUT(rowblock),
            std::vector<casadi_int>& SWIG_OUTPUT(colblock),
            std::vector<casadi_int>& SWIG_OUTPUT(coarse_rowblock),
            std::vector<casadi_int>& SWIG_OUTPUT(coarse_colblock)) const;

    /** \brief Approximate minimal degree preordering
      Fill-reducing ordering applied to the sparsity pattern of a linear system
      prior to factorization.
      The system must be symmetric, for an unsymmetric matrix A, first form the square
      of the pattern, A'*A.

      The implementation is a modified version of cs_amd in CSparse
      Copyright(c) Timothy A. Davis, 2006-2009
      Licensed as a derivative work under the GNU LGPL
    */
    std::vector<casadi_int> amd() const;

#ifndef SWIG
    /** \brief Propagate sparsity through a linear solve
     */
    void spsolve(bvec_t* X, const bvec_t* B, bool tr) const;
#endif // SWIG

    /** \brief Get the location of all non-zero elements as they would appear in a Dense matrix
        A : DenseMatrix  4 x 3
        B : SparseMatrix 4 x 3 , 5 structural non-zeros

        k = A.find()
        A[k] will contain the elements of A that are non-zero in B

        Inverse of `nonzeros`.
    */
    std::vector<casadi_int> find(bool ind1=SWIG_IND1) const;

#ifndef SWIG
    /// Get the location of all nonzero elements (inplace version)
    void find(std::vector<casadi_int>& loc, bool ind1=false) const;
#endif // SWIG

    /** \brief Perform a unidirectional coloring: A greedy distance-2 coloring algorithm
        (Algorithm 3.1 in A. H. GEBREMEDHIN, F. MANNE, A. POTHEN) */
    Sparsity uni_coloring(const Sparsity& AT=Sparsity(),
                          casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Perform a star coloring of a symmetric matrix:
        A greedy distance-2 coloring algorithm
        Algorithm 4.1 in
          What Color Is Your Jacobian? Graph Coloring for Computing Derivatives
          A. H. GEBREMEDHIN, F. MANNE, A. POTHEN
          SIAM Rev., 47(4), 629–705 (2006)

        Ordering options: None (0), largest first (1)
    */
    Sparsity star_coloring(casadi_int ordering = 1,
                            casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Perform a star coloring of a symmetric matrix:
        A new greedy distance-2 coloring algorithm
        Algorithm 4.1 in
          NEW ACYCLIC AND STAR COLORING ALGORITHMS WITH APPLICATION TO COMPUTING HESSIANS
          A. H. GEBREMEDHIN, A. TARAFDAR, F. MANNE, A. POTHEN
          SIAM J. SCI. COMPUT. Vol. 29, No. 3, pp. 1042–1072 (2007)

        Ordering options: None (0), largest first (1)
    */
    Sparsity star_coloring2(casadi_int ordering = 1,
                            casadi_int cutoff = std::numeric_limits<casadi_int>::max()) const;

    /** \brief Order the columns by decreasing degree */
    std::vector<casadi_int> largest_first() const;

    /** \brief Permute rows and/or columns
        Multiply the sparsity with a permutation matrix from the left and/or from the right
        P * A * trans(P), A * trans(P) or A * trans(P) with P defined by an index vector
        containing the row for each col. As an alternative, P can be transposed (inverted).
    */
    Sparsity pmult(const std::vector<casadi_int>& p,
                    bool permute_rows=true, bool permute_columns=true,
                    bool invert_permutation=false) const;

    /// Get the dimension as a string
    std::string dim(bool with_nz=false) const;

    /** \brief Dimension string as a postfix to a name
      Rules:
      1. Dense and scalar: ""
      2. 0-by-0: "[]"
      3. Dense column vector: "[5]"
      4. Dense matrix: "[5x10]"
      5. Otherwise: "[5x10,3nz]"
    */
    std::string postfix_dim() const;

    /// Describe the nonzero location k as a string
    std::string repr_el(casadi_int k) const;

    /** \brief Print a textual representation of sparsity
     */
    void spy(std::ostream &stream=casadi::uout()) const;

    /** \brief Generate a script for Matlab or Octave which visualizes
     * the sparsity using the spy command  */
    void spy_matlab(const std::string& mfile) const;

    /** \brief Export matrix in specific language
     *
     * lang: only 'matlab' supported for now
     * \verbatim
     * options:
     *   inline: Indicates if you want everything on a single line (default: False)
     *   name: Name of exported variable (default: 'sp')
     *   as_matrix: Matlab does not have a sparsity object. (default: false)
    *               With this option true, a numeric matrix will be constructed
     * \endverbatim
     */
    void export_code(const std::string& lang, std::ostream &stream=casadi::uout(),
       const Dict& options=Dict()) const;

    /// Readable name of the public class
    static std::string type_name() {return "Sparsity";}

    // Hash the sparsity pattern
    std::size_t hash() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /** \brief Get KKT system sparsity
     * [H + I1, J'; J, I2] where I1 and I2 are optional
     */
    static Sparsity kkt(const Sparsity& H, const Sparsity& J,
                        bool with_x_diag=true, bool with_lam_g_diag=true);

#ifndef SWIG
    /** \brief Assign the nonzero entries of one sparsity pattern to the nonzero
     * entries of another sparsity pattern */
    template<typename T>
    void set(T* data, const T* val_data, const Sparsity& val_sp) const;

    /** \brief Add the nonzero entries of one sparsity pattern to the nonzero entries
     * of another sparsity pattern */
    template<typename T>
    void add(T* data, const T* val_data, const Sparsity& val_sp) const;

    /** \brief Bitwise or of the nonzero entries of one sparsity pattern and the nonzero
     * entries of another sparsity pattern */
    template<typename T>
    void bor(T* data, const T* val_data, const Sparsity& val_sp) const;

    static std::string file_format(const std::string& filename, const std::string& format_hint);
    static std::set<std::string> file_formats;
  private:

    /// Construct a sparsity pattern from vectors, reuse cached pattern if possible
    void assign_cached(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                      const std::vector<casadi_int>& row, bool order_rows=false);

    /// Construct a sparsity pattern from vectors, reuse cached pattern if possible
    void assign_cached(casadi_int nrow, casadi_int ncol,
                        const casadi_int* colind, const casadi_int* row, bool order_rows=false);

#endif //SWIG
  };

  /** \brief Hash value of an integer */
  template<typename T>
  inline size_t hash_value(T v) { return size_t(v);}

  /** \brief Generate a hash value incrementally (function taken from boost) */
  template<typename T>
  inline void hash_combine(std::size_t& seed, T v) {
    seed ^= hash_value(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /** \brief Generate a hash value incrementally, array */
  inline void hash_combine(std::size_t& seed, const casadi_int* v, std::size_t sz) {
    for (casadi_int i=0; i<sz; ++i) hash_combine(seed, v[i]);
  }

  /** \brief Generate a hash value incrementally (function taken from boost) */
  inline void hash_combine(std::size_t& seed, const std::vector<casadi_int>& v) {
    hash_combine(seed, get_ptr(v), v.size());
  }

  /** \brief Hash a sparsity pattern */
  CASADI_EXPORT std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol,
                                          const std::vector<casadi_int>& colind,
                                          const std::vector<casadi_int>& row);

  CASADI_EXPORT std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol,
                                          const casadi_int* colind,
                                          const casadi_int* row);

#ifndef SWIG
  // Template instantiations
  template<typename DataType>
  void Sparsity::set(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      std::copy(val_data, val_data+sz, data);
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    } else if (val_nel==1) { // if scalar
      std::fill(data, data+sz, val_sz==0 ? DataType(0) : val_data[0]);
    } else if (sz2==val_sz2 && sz1==val_sz1) {
      // Matching dimensions
      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] = val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          } else {
            data[el] = 0;
          }
        }
      }
    } else if (sz1==val_sz2 && sz2==val_sz1 && sz2 == 1) {
      // Assign transposed (this is column)
      const casadi_int* v_cind = val_sp.colind();
      const casadi_int* r = row();
      for (casadi_int el=0; el<sz; ++el) {
        casadi_int rr=r[el];
        data[el] = v_cind[rr]==v_cind[rr+1] ? 0 : val_data[v_cind[rr]];
      }
    } else if (sz1==val_sz2 && sz2==val_sz1 && sz1 == 1) {
      // Assign transposed (this is row)
      for (casadi_int el=0; el<sz; ++el) data[el] = 0;
      const casadi_int* cind = colind();
      const casadi_int* v_r = val_sp.row();
      for (casadi_int el=0; el<val_sz; ++el) {
        casadi_int rr=v_r[el];
        if (cind[rr]!=cind[rr+1]) {
          data[cind[rr]] = val_data[el];
        }
      }
    } else {
      // Make sure that dimension matches
      casadi_error("Sparsity::set<DataType>: shape mismatch. lhs is "
                   + dim() + ", while rhs is " + val_sp.dim() + ".");
    }
  }

  template<typename DataType>
  void Sparsity::add(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();
    const casadi_int nel = sz1*sz2;

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      for (casadi_int k=0; k<sz; ++k) {
        data[k] += val_data[k];
      }
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    }  else if (val_nel==1) { // if scalar
      if (val_sz!=0) {
        for (casadi_int k=0; k<sz; ++k) {
          data[k] += val_data[0];
        }
      }
    } else {
      // Quick return if empty
      if (nel==0 && val_nel==0) return;

      // Make sure that dimension matches
      casadi_assert(sz2==val_sz2 && sz1==val_sz1,
                            "Sparsity::add<DataType>: shape mismatch. lhs is "
                            + dim() + ", while rhs is " + val_sp.dim() + ".");

      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] += val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Sparsity::bor(DataType* data, const DataType* val_data, const Sparsity& val_sp) const {
    // Get dimensions of this
    const casadi_int sz = nnz();
    const casadi_int sz1 = size1();
    const casadi_int sz2 = size2();
    const casadi_int nel = sz1*sz2;

    // Get dimensions of assigning matrix
    const casadi_int val_sz = val_sp.nnz();
    const casadi_int val_sz1 = val_sp.size1();
    const casadi_int val_sz2 = val_sp.size2();
    const casadi_int val_nel = val_sz1*val_sz2;

    // Check if sparsity matches
    if (val_sp==*this) {
      for (casadi_int k=0; k<sz; ++k) {
        data[k] |= val_data[k];
      }
    } else if (this->is_empty()) {
      // Quick return
      return;
    } else if (val_sp.is_empty()) {
      // Quick return
      return;
    }  else if (val_nel==1) { // if scalar
      if (val_sz!=0) {
        for (casadi_int k=0; k<sz; ++k) {
          data[k] |= val_data[0];
        }
      }
    } else {
      // Quick return if empty
      if (nel==0 && val_nel==0) return;

      // Make sure that dimension matches
      casadi_assert(sz2==val_sz2 && sz1==val_sz1,
                            "Sparsity::add<DataType>: shape mismatch. lhs is "
                            + dim() + ", while rhs is " + val_sp.dim() + ".");

      // Sparsity
      const casadi_int* c = row();
      const casadi_int* rind = colind();
      const casadi_int* v_c = val_sp.row();
      const casadi_int* v_rind = val_sp.colind();

      // For all columns
      for (casadi_int i=0; i<sz2; ++i) {

        // Nonzero of the assigning matrix
        casadi_int v_el = v_rind[i];

        // First nonzero of the following column
        casadi_int v_el_end = v_rind[i+1];

        // Next row of the assigning matrix
        casadi_int v_j = v_el<v_el_end ? v_c[v_el] : sz1;

        // Assign all nonzeros
        for (casadi_int el=rind[i]; el!=rind[i+1]; ++el) {

          //  Get row
          casadi_int j=c[el];

          // Forward the assigning nonzero
          while (v_j<j) {
            v_el++;
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }

          // Assign nonzero
          if (v_j==j) {
            data[el] |= val_data[v_el++];
            v_j = v_el<v_el_end ? v_c[v_el] : sz1;
          }
        }
      }
    }
  }

#endif //SWIG

  ///@{
  /// Readability typedefs
  typedef std::map<std::string, Sparsity> SpDict;
  ///@}

} // namespace casadi

#endif // CASADI_SPARSITY_HPP
