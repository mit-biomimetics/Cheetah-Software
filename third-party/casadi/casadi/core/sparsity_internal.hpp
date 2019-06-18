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


#ifndef CASADI_SPARSITY_INTERNAL_HPP
#define CASADI_SPARSITY_INTERNAL_HPP

#include "sparsity.hpp"
#include "shared_object_internal.hpp"
/// \cond INTERNAL

namespace casadi {

  class CASADI_EXPORT SparsityInternal : public SharedObjectInternal {
  private:
    /* \brief Sparsity pattern in compressed column storage (CCS) format
       The first two entries are the number of rows (nrow) and columns (ncol).
       The next (ncol+1) entries are the column offsets (colind). This means that
       the number of nonzeros (nnz) is given as sp_[sp_[1]+2].
       The last nnz entries are the rows of the nonzeros (row). See public class
       for more info about the CCS format used in CasADi. */
    std::vector<casadi_int> sp_;

    /** \brief Structure to hold the block triangular form */
    struct Btf {
      casadi_int nb;
      std::vector<casadi_int> rowperm, colperm;
      std::vector<casadi_int> rowblock, colblock;
      std::vector<casadi_int> coarse_rowblock, coarse_colblock;
    };

    /* \brief The block-triangular factorization for the sparsity
      Calculated on first call, then cached
    */
    mutable Btf* btf_;

  public:
    /// Construct a sparsity pattern from arrays
    SparsityInternal(casadi_int nrow, casadi_int ncol,
                     const casadi_int* colind, const casadi_int* row);

    /// Destructor
    ~SparsityInternal() override;

    /** \brief Get number of rows (see public class) */
    inline const std::vector<casadi_int>& sp() const { return sp_;}

    /** \brief Get number of rows (see public class) */
    inline casadi_int size1() const { return sp_[0];}

    /** \brief Get number of columns (see public class) */
    inline casadi_int size2() const { return sp_[1];}

    /** \brief Get column offsets (see public class) */
    inline const casadi_int* colind() const { return &sp_.front()+2;}

    /** \brief Get row indices (see public class) */
    inline const casadi_int* row() const { return colind()+size2()+1;}

    /// Number of structural non-zeros
    inline casadi_int nnz() const { return colind()[size2()];}

    /** \brief Get the diagonal of the matrix/create a diagonal matrix
     *
     * \param[out] mapping will contain the nonzero mapping
     */
    Sparsity get_diag(std::vector<casadi_int>& mapping) const;

    /// has diagonal entries?
    bool has_diag() const;

    /// Drop diagonal entries
    Sparsity drop_diag() const;

    /** \brief Depth-first search
      * The implementation is a modified version of cs_dfs in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    casadi_int dfs(casadi_int j, casadi_int top, std::vector<casadi_int>& xi,
                   std::vector<casadi_int>& pstack,
                   const std::vector<casadi_int>& pinv, std::vector<bool>& marked) const;

    /** \brief Find the strongly connected components of a square matrix
      * The implementation is a modified version of cs_scc in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    casadi_int scc(std::vector<casadi_int>& p, std::vector<casadi_int>& r) const;

    /** \brief Approximate minimal degree preordering
      * The implementation is a modified version of cs_amd in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    std::vector<casadi_int> amd() const;

    /** \brief Calculate the elimination tree for a matrix
      * len[w] >= ata ? ncol + nrow : ncol
      * len[parent] == ncol
      * Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
      * Modified version of cs_etree in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static void etree(const casadi_int* sp, casadi_int* parent, casadi_int *w, casadi_int ata);

    /** \brief Traverse an elimination tree using depth first search
      * Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
      * Modified version of cs_tdfs in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int postorder_dfs(casadi_int j, casadi_int k, casadi_int* head, casadi_int* next,
                             casadi_int* post, casadi_int* stack);

    /** \brief Calculate the postorder permuation
      * Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
      * len[w] >= 3*n
      * len[post] == n
      * Modified version of cs_post in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static void postorder(const casadi_int* parent, casadi_int n, casadi_int* post, casadi_int* w);

    /** \brief Needed by casadi_qr_colind
      * Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
      * Modified version of cs_leaf in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int leaf(casadi_int i, casadi_int j, const casadi_int* first,
                    casadi_int* maxfirst,
                    casadi_int* prevleaf, casadi_int* ancestor, casadi_int* jleaf);

    /** \brief Calculate the column offsets for the QR R matrix
      * Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
      * len[counts] = ncol
      * len[w] >= 5*ncol + nrow + 1
      * Modified version of cs_counts in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int qr_counts(const casadi_int* tr_sp, const casadi_int* parent,
                         const casadi_int* post, casadi_int* counts, casadi_int* w);

    /** \brief Calculate the number of nonzeros in the QR V matrix
      * Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
      * len[w] >= nrow + 3*ncol
      * len[pinv] == nrow + ncol
      * len[leftmost] == nrow
      * Modified version of cs_sqr in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int qr_nnz(const casadi_int* sp, casadi_int* pinv, casadi_int* leftmost,
                      const casadi_int* parent, casadi_int* nrow_ext, casadi_int* w);

    /** \brief Setup QP solver
      * Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
      * len[w] >= nrow + 7*ncol + 1
      * len[pinv] == nrow + ncol
      * len[leftmost] == nrow
      */
    static void qr_init(const casadi_int* sp, const casadi_int* sp_tr,
                        casadi_int* leftmost, casadi_int* parent, casadi_int* pinv,
                        casadi_int* nrow_ext, casadi_int* v_nnz, casadi_int* r_nnz, casadi_int* w);

    /** \brief Get the row indices for V and R in QR factorization
      * Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
      * Note: nrow <= nrow_ext <= nrow+ncol
      * len[iw] = nrow_ext + ncol
      * len[x] = nrow_ext
      * sp_v = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_v]
      * len[v] nnz_v
      * sp_r = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_r]
      * len[r] nnz_r
      * len[beta] ncol
      * Modified version of cs_qr in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static void qr_sparsities(
                  const casadi_int* sp_a, casadi_int nrow_ext, casadi_int* sp_v, casadi_int* sp_r,
                  const casadi_int* leftmost, const casadi_int* parent, const casadi_int* pinv,
                  casadi_int* iw);

    /** \brief Calculate the column offsets for the L factor of an LDL^T factorization
      * Strictly lower entries only.
      * Ref: User Guide for LDL by Tim Davis
      * len[colind] = ncol+1
      * len[parent] = ncol
      * len[w] >= ncol
      * Modified version of LDL
      * Copyright(c) Timothy A. Davis, 2005-2013
      * Licensed as a derivative work under the GNU LGPL
      */
    static void ldl_colind(const casadi_int* sp, casadi_int* parent,
      casadi_int* l_colind, casadi_int* w);

    /** \brief Calculate the row indices for the L factor of an LDL^T factorization
      * Strictly lower entries only.
      * Ref: User Guide for LDL by Tim Davis
      * len[w] >= n
      * Modified version of LDL
      * Copyright(c) Timothy A. Davis, 2005-2013
      * Licensed as a derivative work under the GNU LGPL
      */
    static void ldl_row(const casadi_int* sp, const casadi_int* parent,
      casadi_int* l_colind, casadi_int* l_row, casadi_int *w);

    /// Transpose the matrix
    Sparsity T() const;

    /** \brief Transpose the matrix and get the reordering of the non-zero entries,
     *
     * \param[out] mapping the non-zeros of the original matrix for each non-zero of the new matrix
     */
    Sparsity transpose(std::vector<casadi_int>& mapping, bool invert_mapping=false) const;

    /// Check if the sparsity is the transpose of another
    bool is_transpose(const SparsityInternal& y) const;

    /// Check if the sparsity is a reshape of another
    bool is_reshape(const SparsityInternal& y) const;

    /** \brief Breadth-first search for coarse decomposition
      * The implementation is a modified version of cs_bfs in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    void bfs(casadi_int n, std::vector<casadi_int>& wi, std::vector<casadi_int>& wj,
             std::vector<casadi_int>& queue, const std::vector<casadi_int>& imatch,
             const std::vector<casadi_int>& jmatch, casadi_int mark) const;

    /** \brief Collect matched columns and rows into p and q
      * The implementation is a modified version of cs_matched in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static void matched(
      casadi_int n, const std::vector<casadi_int>& wj, const std::vector<casadi_int>& imatch,
      std::vector<casadi_int>& p, std::vector<casadi_int>& q, std::vector<casadi_int>& cc,
      std::vector<casadi_int>& rr, casadi_int set, casadi_int mark);

    /** \brief Collect unmatched columns into the permutation vector p
      * The implementation is a modified version of cs_unmatched in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static void unmatched(casadi_int m, const std::vector<casadi_int>& wi,
      std::vector<casadi_int>& p, std::vector<casadi_int>& rr, casadi_int set);

    /** \brief return 1 if column i is in R2
      * The implementation is a modified version of cs_rprune in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int rprune(casadi_int i, casadi_int j, double aij, void *other);

     /** \brief drop entries for which fkeep(A(i, j)) is false; return nz if OK, else -1
       * The implementation is a modified version of cs_fkeep in CSparse
       * Copyright(c) Timothy A. Davis, 2006-2009
       * Licensed as a derivative work under the GNU LGPL
       */
    static casadi_int drop(casadi_int (*fkeep)(casadi_int, casadi_int, double, void *), void *other,
                    casadi_int nrow, casadi_int ncol,
                    std::vector<casadi_int>& colind, std::vector<casadi_int>& row);

    /// Compute the Dulmage-Mendelsohn decomposition
    casadi_int btf(std::vector<casadi_int>& rowperm, std::vector<casadi_int>& colperm,
                          std::vector<casadi_int>& rowblock, std::vector<casadi_int>& colblock,
                          std::vector<casadi_int>& coarse_rowblock,
                          std::vector<casadi_int>& coarse_colblock) const {
      T()->dmperm(colperm, rowperm, colblock, rowblock,
                  coarse_colblock, coarse_rowblock);
      return rowblock.size()-1;
    }

    /// Get cached block triangular form
    const Btf& btf() const;

     /** \brief Compute the Dulmage-Mendelsohn decomposition
       * The implementation is a modified version of cs_dmperm in CSparse
       * Copyright(c) Timothy A. Davis, 2006-2009
       * Licensed as a derivative work under the GNU LGPL
       */
    void dmperm(std::vector<casadi_int>& rowperm, std::vector<casadi_int>& colperm,
                std::vector<casadi_int>& rowblock, std::vector<casadi_int>& colblock,
                std::vector<casadi_int>& coarse_rowblock,
                std::vector<casadi_int>& coarse_colblock) const;

    /** \brief Compute the maximum transversal (maximum matching)
      * The implementation is a modified version of cs_maxtrans in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    void maxtrans(std::vector<casadi_int>& imatch,
                  std::vector<casadi_int>& jmatch, Sparsity& trans, casadi_int seed) const;

    /** \brief Find an augmenting path
      * The implementation is a modified version of cs_augment in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    void augment(casadi_int k, std::vector<casadi_int>& jmatch,
                 casadi_int *cheap, std::vector<casadi_int>& w, casadi_int *js,
                 casadi_int *is, casadi_int *ps) const;

     /** \brief return a random permutation vector
       * return a random permutation vector, the identity perm, or p = n-1:-1:0.
       * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.
       * The implementation is a modified version of cs_randperm in CSparse
       * Copyright(c) Timothy A. Davis, 2006-2009
       * Licensed as a derivative work under the GNU LGPL
       */
    static std::vector<casadi_int> randperm(casadi_int n, casadi_int seed);

    /** \brief Invert a permutation vector */
    static std::vector<casadi_int> invertPermutation(const std::vector<casadi_int>& p);

    /// C = A(p, q) where p and q are permutations of 0..m-1 and 0..n-1.
    Sparsity permute(const std::vector<casadi_int>& pinv,
      const std::vector<casadi_int>& q, casadi_int values) const;

    /** \brief C = A(p, q) where p and q are permutations of 0..m-1 and 0..n-1
      * The implementation is a modified version of cs_permute in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    void permute(const std::vector<casadi_int>& pinv,
                 const std::vector<casadi_int>& q, casadi_int values,
                 std::vector<casadi_int>& colind_C,
                 std::vector<casadi_int>& row_C) const;

    /** \brief clear w
      * The implementation is a modified version of cs_wclear in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int wclear(casadi_int mark, casadi_int lemax, casadi_int *w, casadi_int n);

    /** \brief keep off-diagonal entries; drop diagonal entries
      * The implementation is a modified version of cs_diag in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    static casadi_int diag(casadi_int i, casadi_int j, double aij, void *other);

    /** \brief C = A*B
      * The implementation is a modified version of cs_multiply in CSparse
      * Copyright(c) Timothy A. Davis, 2006-2009
      * Licensed as a derivative work under the GNU LGPL
      */
    Sparsity multiply(const Sparsity& B) const;

     /** \brief x = x + beta * A(:, j), where x is a dense vector and A(:, j) is sparse
       * The implementation is a modified version of cs_scatter in CSparse
       * Copyright(c) Timothy A. Davis, 2006-2009
       * Licensed as a derivative work under the GNU LGPL
       */
    casadi_int scatter(casadi_int j, std::vector<casadi_int>& w, casadi_int mark,
                       casadi_int* Ci, casadi_int nz) const;

    /// Get row() as a vector
    std::vector<casadi_int> get_row() const;

    /// Get colind() as a vector
    std::vector<casadi_int> get_colind() const;

    /// Get the column for each nonzero
    std::vector<casadi_int> get_col() const;

    /// Resize
    Sparsity _resize(casadi_int nrow, casadi_int ncol) const;

    /// Reshape a sparsity, order of nonzeros remains the same
    Sparsity _reshape(casadi_int nrow, casadi_int ncol) const;

    /// Number of elements
    casadi_int numel() const;

    /// Number of non-zeros in the lower triangular half
    casadi_int nnz_lower(bool strictly=false) const;

    /// Number of non-zeros in the upper triangular half
    casadi_int nnz_upper(bool strictly=false) const;

    /// Number of non-zeros on the diagonal
    casadi_int nnz_diag() const;

    /** \brief Upper half-bandwidth */
    casadi_int bw_upper() const;

    /** \brief Lower half-bandwidth */
    casadi_int bw_lower() const;

    /// Shape
    std::pair<casadi_int, casadi_int> size() const;

    /// Is scalar?
    bool is_scalar(bool scalar_and_dense) const;

    /** \brief Check if the sparsity is empty
     *
     * A sparsity is considered empty if one of the dimensions is zero
     * (or optionally both dimensions)
     */
    bool is_empty(bool both=false) const;

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

    /// Is lower triangular?
    bool is_tril() const;

    /// is upper triangular?
    bool is_triu() const;

    /// Get upper triangular part
    Sparsity _triu(bool includeDiagonal) const;

    /// Get lower triangular part
    Sparsity _tril(bool includeDiagonal) const;

    /// Get nonzeros in lower triangular part
    std::vector<casadi_int> get_lower() const;

    /// Get nonzeros in upper triangular part
    std::vector<casadi_int> get_upper() const;

    /// Get the dimension as a string
    std::string dim(bool with_nz=false) const;

    /// Describe the nonzero location k as a string
    std::string repr_el(casadi_int k) const;

    /// Sparsity pattern for a matrix-matrix product (details in public class)
    Sparsity _mtimes(const Sparsity& y) const;

    ///@{
    /// Union of two sparsity patterns
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool function0_is_zero,
                            std::vector<unsigned char>& mapping) const;
    Sparsity combine(const Sparsity& y, bool f0x_is_zero, bool function0_is_zero) const;

    template<bool with_mapping>
    Sparsity combineGen1(const Sparsity& y, bool f0x_is_zero, bool function0_is_zero,
                                std::vector<unsigned char>& mapping) const;

    template<bool with_mapping, bool f0x_is_zero, bool function0_is_zero>
    Sparsity combineGen(const Sparsity& y, std::vector<unsigned char>& mapping) const;
    ///@}

    /// Take the inverse of a sparsity pattern; flip zeros and non-zeros
    Sparsity pattern_inverse() const;

    /// Check if two sparsity patterns are the same
    bool is_equal(const Sparsity& y) const;

    /// Check if two sparsity patterns are the same
    bool is_equal(casadi_int y_nrow, casadi_int y_ncol, const std::vector<casadi_int>& y_colind,
                 const std::vector<casadi_int>& y_row) const;

    /// Check if two sparsity patterns are the same
    bool is_equal(casadi_int y_nrow, casadi_int y_ncol,
                  const casadi_int* y_colind, const casadi_int* y_row) const;

    /// Check if pattern is repeated
    bool is_stacked(const Sparsity& y, casadi_int n) const;

    /// Enlarge the matrix along the first dimension (i.e. insert rows)
    Sparsity _enlargeRows(casadi_int nrow, const std::vector<casadi_int>& rr, bool ind1) const;

    /// Enlarge the matrix along the second dimension (i.e. insert columns)
    Sparsity _enlargeColumns(casadi_int ncol, const std::vector<casadi_int>& cc, bool ind1) const;

    /// Make a patten dense
    Sparsity makeDense(std::vector<casadi_int>& mapping) const;

    /// Erase rows and/or columns - does bounds checking
    Sparsity _erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                      bool ind1, std::vector<casadi_int>& mapping) const;

    /// Erase elements
    Sparsity _erase(const std::vector<casadi_int>& rr, bool ind1,
                      std::vector<casadi_int>& mapping) const;

    /// Append another sparsity patten vertically (vectors only)
    Sparsity _appendVector(const SparsityInternal& sp) const;

    /// Append another sparsity patten horizontally
    Sparsity _appendColumns(const SparsityInternal& sp) const;

    /** \brief Get a submatrix
    * Does bounds checking
    * rr and rr are not required to be monotonous
    */
    Sparsity sub(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                 std::vector<casadi_int>& mapping, bool ind1) const;

    /** \brief Get a set of elements
    * Does bounds checking
    * rr is not required to be monotonous
    */
    Sparsity sub(const std::vector<casadi_int>& rr, const SparsityInternal& sp,
                 std::vector<casadi_int>& mapping, bool ind1) const;

    /// Get the index of an existing non-zero element
    casadi_int get_nz(casadi_int rr, casadi_int cc) const;

    /// Get a set of non-zero element - does bounds checking
    std::vector<casadi_int> get_nz(const std::vector<casadi_int>& rr,
                                    const std::vector<casadi_int>& cc) const;

    /// Get the nonzero index for a set of elements (see description in public class)
    void get_nz(std::vector<casadi_int>& indices) const;

    /// Does the rows appear sequentially on each col
    bool rowsSequential(bool strictly) const;

    /** \brief Remove duplicate entries
     *
     * The same indices will be removed from the mapping vector,
     * which must have the same length as the number of nonzeros
     */
    Sparsity _removeDuplicates(std::vector<casadi_int>& mapping) const;

    /// Get element index for each nonzero
    void find(std::vector<casadi_int>& loc, bool ind1) const;

    /// Hash the sparsity pattern
    std::size_t hash() const;

    /// Readable name of the internal class
    std::string class_name() const override {return "SparsityInternal";}

    /// Print description
    void disp(std::ostream& stream, bool more) const override;

    /** \brief Perform a unidirectional coloring
     *
     * A greedy distance-2 coloring algorithm
     * (Algorithm 3.1 in A. H. GEBREMEDHIN, F. MANNE, A. POTHEN)
     */
    Sparsity uni_coloring(const Sparsity& AT, casadi_int cutoff) const;

    /** \brief A greedy distance-2 coloring algorithm
     * See description in public class.
     */
    Sparsity star_coloring(casadi_int ordering, casadi_int cutoff) const;

    /** \brief An improved distance-2 coloring algorithm
     * See description in public class.
     */
    Sparsity star_coloring2(casadi_int ordering, casadi_int cutoff) const;

    /// Order the columns by decreasing degree
    std::vector<casadi_int> largest_first() const;

    /// Permute rows and/or columns
    Sparsity pmult(const std::vector<casadi_int>& p, bool permute_rows=true, bool permute_cols=true,
                   bool invert_permutation=false) const;

    /** \brief Print a textual representation of sparsity */
    void spy(std::ostream &stream) const;

    /// Generate a script for Matlab or Octave which visualizes the sparsity using the spy command
    void spy_matlab(const std::string& mfile) const;

    /** \brief Export sparsity in Matlab format */
    void export_code(const std::string& lang, std::ostream &stream,
       const Dict& options) const;

    /// Propagate sparsity through a linear solve
    void spsolve(bvec_t* X, const bvec_t* B, bool tr) const;
};

} // namespace casadi
/// \endcond

#endif // CASADI_SPARSITY_INTERNAL_HPP
