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


#include "sparsity_internal.hpp"
#include "matrix.hpp"
#include "casadi_misc.hpp"
#include "sparse_storage_impl.hpp"
#include <climits>

#define CASADI_THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in Sparsity::" FNAME " at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

using namespace std;

namespace casadi {
  // Instantiate templates
  template class SparseStorage<Sparsity>;

  /// \cond INTERNAL
  // Singletons
  class EmptySparsity : public Sparsity {
  public:
    EmptySparsity() {
      const casadi_int colind[1] = {0};
      own(new SparsityInternal(0, 0, colind, nullptr));
    }
  };

  class ScalarSparsity : public Sparsity {
  public:
    ScalarSparsity() {
      const casadi_int colind[2] = {0, 1};
      const casadi_int row[1] = {0};
      own(new SparsityInternal(1, 1, colind, row));
    }
  };

  class ScalarSparseSparsity : public Sparsity {
  public:
    ScalarSparseSparsity() {
      const casadi_int colind[2] = {0, 0};
      const casadi_int row[1] = {0};
      own(new SparsityInternal(1, 1, colind, row));
    }
  };
  /// \endcond

  Sparsity::Sparsity(casadi_int dummy) {
    casadi_assert_dev(dummy==0);
  }

  Sparsity Sparsity::create(SparsityInternal *node) {
    Sparsity ret;
    ret.own(node);
    return ret;
  }

  Sparsity::Sparsity(casadi_int nrow, casadi_int ncol) {
    casadi_assert_dev(nrow>=0);
    casadi_assert_dev(ncol>=0);
    std::vector<casadi_int> row, colind(ncol+1, 0);
    assign_cached(nrow, ncol, colind, row);
  }

  Sparsity::Sparsity(const std::pair<casadi_int, casadi_int>& rc) {
    casadi_assert_dev(rc.first>=0);
    casadi_assert_dev(rc.second>=0);
    std::vector<casadi_int> row, colind(rc.second+1, 0);
    assign_cached(rc.first, rc.second, colind, row);
  }

  Sparsity::Sparsity(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                     const std::vector<casadi_int>& row, bool order_rows) {
    casadi_assert_dev(nrow>=0);
    casadi_assert_dev(ncol>=0);
    assign_cached(nrow, ncol, colind, row, order_rows);
  }

  Sparsity::Sparsity(casadi_int nrow, casadi_int ncol,
      const casadi_int* colind, const casadi_int* row, bool order_rows) {
    casadi_assert_dev(nrow>=0);
    casadi_assert_dev(ncol>=0);
    if (colind==nullptr || colind[ncol]==nrow*ncol) {
      *this = dense(nrow, ncol);
    } else {
      vector<casadi_int> colindv(colind, colind+ncol+1);
      vector<casadi_int> rowv(row, row+colind[ncol]);
      assign_cached(nrow, ncol, colindv, rowv, order_rows);
    }
  }

  const SparsityInternal* Sparsity::operator->() const {
    return static_cast<const SparsityInternal*>(SharedObject::operator->());
  }

  const SparsityInternal& Sparsity::operator*() const {
    return *static_cast<const SparsityInternal*>(get());
  }

  bool Sparsity::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const SparsityInternal*>(ptr)!=nullptr;
  }

  casadi_int Sparsity::size1() const {
    return (*this)->size1();
  }

  casadi_int Sparsity::size2() const {
    return (*this)->size2();
  }

  casadi_int Sparsity::numel() const {
    return (*this)->numel();
  }

  double Sparsity::density() const {
    double r = 100;
    r *= static_cast<double>(nnz());
    r /= static_cast<double>(size1());
    r /= static_cast<double>(size2());
    return r;
  }

  bool Sparsity::is_empty(bool both) const {
    return (*this)->is_empty(both);
  }

  casadi_int Sparsity::nnz() const {
    return (*this)->nnz();
  }

  std::pair<casadi_int, casadi_int> Sparsity::size() const {
    return (*this)->size();
  }

  casadi_int Sparsity::size(casadi_int axis) const {
    switch (axis) {
    case 1: return size1();
    case 2: return size2();
    }
    casadi_error("Axis must be 1 or 2.");
  }

  const casadi_int* Sparsity::row() const {
    return (*this)->row();
  }

  const casadi_int* Sparsity::colind() const {
    return (*this)->colind();
  }

  casadi_int Sparsity::row(casadi_int el) const {
    if (el<0 || el>=nnz()) {
      throw std::out_of_range("Sparsity::row: Index " + str(el)
        + " out of range [0," + str(nnz()) + ")");
    }
    return row()[el];
  }

  casadi_int Sparsity::colind(casadi_int cc) const {
    if (cc<0 || cc>size2()) {
      throw std::out_of_range("Sparsity::colind: Index "
        + str(cc) + " out of range [0," + str(size2()) + "]");
    }
    return colind()[cc];
  }

  void Sparsity::resize(casadi_int nrow, casadi_int ncol) {
    if (size1()!=nrow || size2() != ncol) {
      *this = (*this)->_resize(nrow, ncol);
    }
  }

  casadi_int Sparsity::add_nz(casadi_int rr, casadi_int cc) {
    // If negative index, count from the back
    if (rr<0) rr += size1();
    if (cc<0) cc += size2();

    // Check consistency
    casadi_assert(rr>=0 && rr<size1(), "Row index out of bounds");
    casadi_assert(cc>=0 && cc<size2(), "Column index out of bounds");

    // Quick return if matrix is dense
    if (is_dense()) return rr+cc*size1();

    // Get sparsity pattern
    casadi_int size1=this->size1(), size2=this->size2(), nnz=this->nnz();
    const casadi_int *colind = this->colind(), *row = this->row();

    // Quick return if we are adding an element to the end
    if (colind[cc]==nnz || (colind[cc+1]==nnz && row[nnz-1]<rr)) {
      std::vector<casadi_int> rowv(nnz+1);
      copy(row, row+nnz, rowv.begin());
      rowv[nnz] = rr;
      std::vector<casadi_int> colindv(colind, colind+size2+1);
      for (casadi_int c=cc; c<size2; ++c) colindv[c+1]++;
      assign_cached(size1, size2, colindv, rowv);
      return rowv.size()-1;
    }

    // go to the place where the element should be
    casadi_int ind;
    for (ind=colind[cc]; ind<colind[cc+1]; ++ind) { // better: loop from the back to the front
      if (row[ind] == rr) {
        return ind; // element exists
      } else if (row[ind] > rr) {
        break;                // break at the place where the element should be added
      }
    }

    // insert the element
    std::vector<casadi_int> rowv = get_row(), colindv = get_colind();
    rowv.insert(rowv.begin()+ind, rr);
    for (casadi_int c=cc+1; c<size2+1; ++c) colindv[c]++;

    // Return the location of the new element
    assign_cached(size1, size2, colindv, rowv);
    return ind;
  }

  bool Sparsity::has_nz(casadi_int rr, casadi_int cc) const {
    return get_nz(rr, cc)!=-1;
  }


  casadi_int Sparsity::get_nz(casadi_int rr, casadi_int cc) const {
    return (*this)->get_nz(rr, cc);
  }

  Sparsity Sparsity::reshape(const Sparsity& x, const Sparsity& sp) {
    casadi_assert_dev(x.is_reshape(sp));
    return sp;
  }

  Sparsity Sparsity::reshape(const Sparsity& x, casadi_int nrow, casadi_int ncol) {
    return x->_reshape(nrow, ncol);
  }

  std::vector<casadi_int> Sparsity::get_nz(const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc) const {
    return (*this)->get_nz(rr, cc);
  }

  bool Sparsity::is_scalar(bool scalar_and_dense) const {
    return (*this)->is_scalar(scalar_and_dense);
  }

  bool Sparsity::is_dense() const {
    return (*this)->is_dense();
  }

  bool Sparsity::is_diag() const {
    return (*this)->is_diag();
  }

  bool Sparsity::is_row() const {
    return (*this)->is_row();
  }

  bool Sparsity::is_column() const {
    return (*this)->is_column();
  }

  bool Sparsity::is_vector() const {
    return (*this)->is_vector();
  }

  bool Sparsity::is_square() const {
    return (*this)->is_square();
  }

  bool Sparsity::is_symmetric() const {
    return (*this)->is_symmetric();
  }

  bool Sparsity::is_tril() const {
    return (*this)->is_tril();
  }

  bool Sparsity::is_triu() const {
    return (*this)->is_triu();
  }

  Sparsity Sparsity::sub(const std::vector<casadi_int>& rr, const Sparsity& sp,
                         std::vector<casadi_int>& mapping, bool ind1) const {
    return (*this)->sub(rr, *sp, mapping, ind1);
  }

  Sparsity Sparsity::sub(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                         std::vector<casadi_int>& mapping, bool ind1) const {
    return (*this)->sub(rr, cc, mapping, ind1);
  }

  std::vector<casadi_int> Sparsity::erase(const std::vector<casadi_int>& rr,
                                  const std::vector<casadi_int>& cc, bool ind1) {
    vector<casadi_int> mapping;
    *this = (*this)->_erase(rr, cc, ind1, mapping);
    return mapping;
  }

  std::vector<casadi_int> Sparsity::erase(const std::vector<casadi_int>& rr, bool ind1) {
    vector<casadi_int> mapping;
    *this = (*this)->_erase(rr, ind1, mapping);
    return mapping;
  }

  casadi_int Sparsity::nnz_lower(bool strictly) const {
    return (*this)->nnz_lower(strictly);
  }

  casadi_int Sparsity::nnz_upper(bool strictly) const {
    return (*this)->nnz_upper(strictly);
  }

  casadi_int Sparsity::nnz_diag() const {
    return (*this)->nnz_diag();
  }

  std::vector<casadi_int> Sparsity::get_colind() const {
    return (*this)->get_colind();
  }

  std::vector<casadi_int> Sparsity::get_col() const {
    return (*this)->get_col();
  }

  std::vector<casadi_int> Sparsity::get_row() const {
    return (*this)->get_row();
  }

  void Sparsity::get_ccs(std::vector<casadi_int>& colind, std::vector<casadi_int>& row) const {
    colind = get_colind();
    row = get_row();
  }

  void Sparsity::get_crs(std::vector<casadi_int>& rowind, std::vector<casadi_int>& col) const {
    T().get_ccs(rowind, col);
  }

  void Sparsity::get_triplet(std::vector<casadi_int>& row, std::vector<casadi_int>& col) const {
    row = get_row();
    col = get_col();
  }

  Sparsity Sparsity::transpose(std::vector<casadi_int>& mapping, bool invert_mapping) const {
    return (*this)->transpose(mapping, invert_mapping);
  }

  Sparsity Sparsity::T() const {
    return (*this)->T();
  }

  Sparsity Sparsity::combine(const Sparsity& y, bool f0x_is_zero,
                                    bool function0_is_zero,
                                    std::vector<unsigned char>& mapping) const {
    return (*this)->combine(y, f0x_is_zero, function0_is_zero, mapping);
  }

  Sparsity Sparsity::combine(const Sparsity& y, bool f0x_is_zero,
                                    bool function0_is_zero) const {
    return (*this)->combine(y, f0x_is_zero, function0_is_zero);
  }

  Sparsity Sparsity::unite(const Sparsity& y, std::vector<unsigned char>& mapping) const {
    return (*this)->combine(y, false, false, mapping);
  }

  Sparsity Sparsity::unite(const Sparsity& y) const {
    return (*this)->combine(y, false, false);
  }

  Sparsity Sparsity::intersect(const Sparsity& y,
                                         std::vector<unsigned char>& mapping) const {
    return (*this)->combine(y, true, true, mapping);
  }

  Sparsity Sparsity::intersect(const Sparsity& y) const {
    return (*this)->combine(y, true, true);
  }

  Sparsity Sparsity::mtimes(const Sparsity& x, const Sparsity& y) {
    // Check matching dimensions
    casadi_assert(x.size2()==y.size1(),
      "Matrix product with incompatible dimensions. Lhs is "
      + x.dim() + " and rhs is " + y.dim() + ".");

    return x->_mtimes(y);
  }

  bool Sparsity::is_stacked(const Sparsity& y, casadi_int n) const {
    return (*this)->is_stacked(y, n);
  }

  bool Sparsity::is_equal(const Sparsity& y) const {
    return (*this)->is_equal(y);
  }

  bool Sparsity::is_equal(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                         const std::vector<casadi_int>& row) const {
    return (*this)->is_equal(nrow, ncol, colind, row);
  }

  bool Sparsity::is_equal(casadi_int nrow, casadi_int ncol,
      const casadi_int* colind, const casadi_int* row) const {
    return (*this)->is_equal(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::operator+(const Sparsity& b) const {
    return unite(b);
  }

  Sparsity Sparsity::operator*(const Sparsity& b) const {
    std::vector< unsigned char > mapping;
    return intersect(b, mapping);
  }

  Sparsity Sparsity::pattern_inverse() const {
    return (*this)->pattern_inverse();
  }

  void Sparsity::append(const Sparsity& sp) {
    if (sp.size1()==0 && sp.size2()==0) {
      // Appending pattern is empty
      return;
    } else if (size1()==0 && size2()==0) {
      // This is empty
      *this = sp;
    } else {
      casadi_assert(size2()==sp.size2(),
                            "Sparsity::append: Dimension mismatch. "
                            "You attempt to append a shape " + sp.dim()
                            + " to a shape " + dim()
                            + ". The number of columns must match.");
      if (sp.size1()==0) {
        // No rows to add
        return;
      } else if (size1()==0) {
        // No rows before
        *this = sp;
      } else if (is_column()) {
        // Append to vector (inefficient)
        *this = (*this)->_appendVector(*sp);
      } else {
        // Append to matrix (inefficient)
        *this = vertcat({*this, sp});
      }
    }
  }

  void Sparsity::appendColumns(const Sparsity& sp) {
    if (sp.size1()==0 && sp.size2()==0) {
      // Appending pattern is empty
      return;
    } else if (size1()==0 && size2()==0) {
      // This is empty
      *this = sp;
    } else {
      casadi_assert(size1()==sp.size1(),
                            "Sparsity::appendColumns: Dimension mismatch. You attempt to "
                            "append a shape " + sp.dim() + " to a shape "
                            + dim() + ". The number of rows must match.");
      if (sp.size2()==0) {
        // No columns to add
        return;
      } else if (size2()==0) {
        // No columns before
        *this = sp;
      } else {
        // Append to matrix (expensive)
        *this = (*this)->_appendColumns(*sp);
      }
    }
  }

  Sparsity::CachingMap& Sparsity::getCache() {
    static CachingMap ret;
    return ret;
  }

  const Sparsity& Sparsity::getScalar() {
    static ScalarSparsity ret;
    return ret;
  }

  const Sparsity& Sparsity::getScalarSparse() {
    static ScalarSparseSparsity ret;
    return ret;
  }

  const Sparsity& Sparsity::getEmpty() {
    static EmptySparsity ret;
    return ret;
  }

  void Sparsity::enlarge(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& rr,
                         const std::vector<casadi_int>& cc, bool ind1) {
    enlargeColumns(ncol, cc, ind1);
    enlargeRows(nrow, rr, ind1);
  }

  void Sparsity::enlargeColumns(casadi_int ncol, const std::vector<casadi_int>& cc, bool ind1) {
    casadi_assert_dev(cc.size() == size2());
    if (cc.empty()) {
      *this = Sparsity(size1(), ncol);
    } else {
      *this = (*this)->_enlargeColumns(ncol, cc, ind1);
    }
  }

  void Sparsity::enlargeRows(casadi_int nrow, const std::vector<casadi_int>& rr, bool ind1) {
    casadi_assert_dev(rr.size() == size1());
    if (rr.empty()) {
      *this = Sparsity(nrow, size2());
    } else {
      *this = (*this)->_enlargeRows(nrow, rr, ind1);
    }
  }

  Sparsity Sparsity::diag(casadi_int nrow, casadi_int ncol) {
    // Smallest dimension
    casadi_int n = min(nrow, ncol);

    // Column offset
    vector<casadi_int> colind(ncol+1, n);
    for (casadi_int cc=0; cc<n; ++cc) colind[cc] = cc;

    // Row
    vector<casadi_int> row = range(n);

    // Create pattern from vectors
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::makeDense(std::vector<casadi_int>& mapping) const {
    return (*this)->makeDense(mapping);
  }

  std::string Sparsity::dim(bool with_nz) const {
    return (*this)->dim(with_nz);
  }

  std::string Sparsity::postfix_dim() const {
    if (is_dense()) {
      if (is_scalar()) {
        return "";
      } else if (is_empty(true)) {
        return "[]";
      } else if (is_column()) {
        return "[" + str(size1()) + "]";
      } else {
        return "[" + dim(false) + "]";
      }
    } else {
      return "[" + dim(true) + "]";
    }
  }

  std::string Sparsity::repr_el(casadi_int k) const {
    return (*this)->repr_el(k);
  }

  Sparsity Sparsity::get_diag(std::vector<casadi_int>& mapping) const {
    return (*this)->get_diag(mapping);
  }

  std::vector<casadi_int> Sparsity::etree(bool ata) const {
    vector<casadi_int> parent(size2()), w(size1() + size2());
    SparsityInternal::etree(*this, get_ptr(parent), get_ptr(w), ata);
    return parent;
  }

  Sparsity Sparsity::ldl(std::vector<casadi_int>& p, bool amd) const {
    casadi_assert(is_symmetric(),
                 "LDL factorization requires a symmetric matrix");
    // Recursive call if AMD
    if (amd) {
      // Get AMD reordering
      p = this->amd();
      // Permute sparsity pattern
      std::vector<casadi_int> tmp;
      Sparsity Aperm = sub(p, p, tmp);
      // Call recursively
      return Aperm.ldl(tmp, false);
    }
    // Dimension
    casadi_int n=size1();
    // Natural ordering
    p = range(n);
    // Work vector
    std::vector<casadi_int> w(3*n);
    // Elimination tree
    std::vector<casadi_int> parent(n);
    // Calculate colind in L (strictly lower entries only)
    std::vector<casadi_int> L_colind(1+n);
    SparsityInternal::ldl_colind(*this, get_ptr(parent), get_ptr(L_colind), get_ptr(w));
    // Get rows in L (strictly lower entries only)
    std::vector<casadi_int> L_row(L_colind.back());
    SparsityInternal::ldl_row(*this, get_ptr(parent), get_ptr(L_colind), get_ptr(L_row),
                    get_ptr(w));
    // Sparsity of L^T
    return Sparsity(n, n, L_colind, L_row, true).T();
  }

  void Sparsity::
  qr_sparse(Sparsity& V, Sparsity& R, std::vector<casadi_int>& prinv,
            std::vector<casadi_int>& pc, bool amd) const {
    // Dimensions
    casadi_int size1=this->size1(), size2=this->size2();

    // Recursive call if AMD
    if (amd) {
      // Get AMD reordering
      pc = mtimes(T(), *this).amd();
      // Permute sparsity pattern
      std::vector<casadi_int> tmp;
      Sparsity Aperm = sub(range(size1), pc, tmp);
      // Call recursively
      return Aperm.qr_sparse(V, R, prinv, tmp, false);
    }

    // No column permutation
    pc = range(size2);

    // Allocate memory
    vector<casadi_int> leftmost(size1);
    vector<casadi_int> parent(size2);
    prinv.resize(size1 + size2);
    vector<casadi_int> iw(size1 + 7*size2 + 1);

    // Initialize QP solve
    casadi_int nrow_ext, v_nnz, r_nnz;
    SparsityInternal::qr_init(*this, T(),
                              get_ptr(leftmost), get_ptr(parent), get_ptr(prinv),
                              &nrow_ext, &v_nnz, &r_nnz, get_ptr(iw));

    // Calculate sparsities
    vector<casadi_int> sp_v(2 + size2 + 1 + v_nnz);
    vector<casadi_int> sp_r(2 + size2 + 1 + r_nnz);
    SparsityInternal::qr_sparsities(*this, nrow_ext, get_ptr(sp_v), get_ptr(sp_r),
                                    get_ptr(leftmost), get_ptr(parent), get_ptr(prinv),
                                    get_ptr(iw));
    prinv.resize(nrow_ext);
    V = compressed(sp_v, true);
    R = compressed(sp_r, true);
  }

  casadi_int Sparsity::dfs(casadi_int j, casadi_int top, std::vector<casadi_int>& xi,
                            std::vector<casadi_int>& pstack,
                            const std::vector<casadi_int>& pinv,
                            std::vector<bool>& marked) const {
    return (*this)->dfs(j, top, xi, pstack, pinv, marked);
  }

  casadi_int Sparsity::scc(std::vector<casadi_int>& p, std::vector<casadi_int>& r) const {
    return (*this)->scc(p, r);
  }

  std::vector<casadi_int> Sparsity::amd() const {
    return (*this)->amd();
  }

  casadi_int Sparsity::btf(std::vector<casadi_int>& rowperm, std::vector<casadi_int>& colperm,
                            std::vector<casadi_int>& rowblock, std::vector<casadi_int>& colblock,
                            std::vector<casadi_int>& coarse_rowblock,
                            std::vector<casadi_int>& coarse_colblock) const {
    try {
      return (*this)->btf(rowperm, colperm, rowblock, colblock,
                          coarse_rowblock, coarse_colblock);
    } catch (exception &e) {
      CASADI_THROW_ERROR("btf", e.what());
    }
  }

  void Sparsity::spsolve(bvec_t* X, const bvec_t* B, bool tr) const {
    (*this)->spsolve(X, B, tr);
  }

  bool Sparsity::rowsSequential(bool strictly) const {
    return (*this)->rowsSequential(strictly);
  }

  void Sparsity::removeDuplicates(std::vector<casadi_int>& mapping) {
    *this = (*this)->_removeDuplicates(mapping);
  }

  std::vector<casadi_int> Sparsity::find(bool ind1) const {
    std::vector<casadi_int> loc;
    find(loc, ind1);
    return loc;
  }

  void Sparsity::find(std::vector<casadi_int>& loc, bool ind1) const {
    (*this)->find(loc, ind1);
  }

  void Sparsity::get_nz(std::vector<casadi_int>& indices) const {
    (*this)->get_nz(indices);
  }

  Sparsity Sparsity::uni_coloring(const Sparsity& AT, casadi_int cutoff) const {
    if (AT.is_null()) {
      return (*this)->uni_coloring(T(), cutoff);
    } else {
      return (*this)->uni_coloring(AT, cutoff);
    }
  }

  Sparsity Sparsity::star_coloring(casadi_int ordering, casadi_int cutoff) const {
    return (*this)->star_coloring(ordering, cutoff);
  }

  Sparsity Sparsity::star_coloring2(casadi_int ordering, casadi_int cutoff) const {
    return (*this)->star_coloring2(ordering, cutoff);
  }

  std::vector<casadi_int> Sparsity::largest_first() const {
    return (*this)->largest_first();
  }

  Sparsity Sparsity::pmult(const std::vector<casadi_int>& p, bool permute_rows,
                            bool permute_columns, bool invert_permutation) const {
    return (*this)->pmult(p, permute_rows, permute_columns, invert_permutation);
  }

  void Sparsity::spy_matlab(const std::string& mfile) const {
    (*this)->spy_matlab(mfile);
  }

  void Sparsity::export_code(const std::string& lang, std::ostream &stream,
      const Dict& options) const {
    (*this)->export_code(lang, stream, options);
  }

  void Sparsity::spy(std::ostream &stream) const {
    (*this)->spy(stream);
  }

  bool Sparsity::is_transpose(const Sparsity& y) const {
    return (*this)->is_transpose(*y);
  }

  bool Sparsity::is_reshape(const Sparsity& y) const {
    return (*this)->is_reshape(*y);
  }

  std::size_t Sparsity::hash() const {
    return (*this)->hash();
  }

  void Sparsity::assign_cached(casadi_int nrow, casadi_int ncol,
                                const std::vector<casadi_int>& colind,
                                const std::vector<casadi_int>& row, bool order_rows) {
    casadi_assert_dev(colind.size()==ncol+1);
    casadi_assert_dev(row.size()==colind.back());
    assign_cached(nrow, ncol, get_ptr(colind), get_ptr(row), order_rows);
  }

  void Sparsity::assign_cached(casadi_int nrow, casadi_int ncol,
      const casadi_int* colind, const casadi_int* row, bool order_rows) {
    // Scalars and empty patterns are handled separately
    if (ncol==0 && nrow==0) {
      // If empty
      *this = getEmpty();
      return;
    } else if (ncol==1 && nrow==1) {
      if (colind[ncol]==0) {
        // If sparse scalar
        *this = getScalarSparse();
        return;
      } else {
        // If dense scalar
        *this = getScalar();
        return;
      }
    }

    // Make sure colind starts with zero
    casadi_assert(colind[0]==0,
                  "Compressed Column Storage is not sane. "
                  "First element of colind must be zero.");

    // Make sure colind is montone
    for (casadi_int c=0; c<ncol; c++) {
      casadi_assert(colind[c+1]>=colind[c],
                    "Compressed Column Storage is not sane. "
                    "colind must be monotone.");
    }

    // Check if rows correct and ordered without duplicates
    bool rows_ordered = true;
    for (casadi_int c=0; c<ncol; ++c) {
      casadi_int last_r = -1;
      for (casadi_int k=colind[c]; k<colind[c+1]; ++k) {
        casadi_int r = row[k];
        // Make sure values are in within the [0,ncol) range
        casadi_assert(r>=0 && r<nrow,
                      "Compressed Column Storage is not sane.\n"
                      "row[ " + str(k) + " == " + str(r) + " not in range "
                      "[0, " + str(nrow) + ")");
        // Check if ordered
        if (r<=last_r) rows_ordered = false;
        last_r = r;
      }
    }

    // If unordered, need to order
    if (!rows_ordered) {
      casadi_assert(order_rows,
                    "Compressed Column Storage is not sane.\n"
                    "Row indices not strictly monotonically increasing for each "
                    "column and reordering is not enabled.");
      // Number of nonzeros
      casadi_int nnz = colind[ncol];
      // Get all columns
      vector<casadi_int> col(nnz);
      for (casadi_int c=0; c<ncol; ++c) {
        for (casadi_int k=colind[c]; k<colind[c+1]; ++k) col[k] = c;
      }
      // Make sane via triplet format
      *this = triplet(nrow, ncol, vector<casadi_int>(row, row+nnz), col);
      return;
    }

    // Hash the pattern
    std::size_t h = hash_sparsity(nrow, ncol, colind, row);

    // Get a reference to the cache
    CachingMap& cache = getCache();

    // Record the current number of buckets (for garbage collection below)
    casadi_int bucket_count_before = cache.bucket_count();

    // WORKAROUND, functions do not appear to work when bucket_count==0
    if (bucket_count_before>0) {

      // Find the range of patterns equal to the key (normally only zero or one)
      pair<CachingMap::iterator, CachingMap::iterator> eq = cache.equal_range(h);

      // Loop over maching patterns
      for (CachingMap::iterator i=eq.first; i!=eq.second; ++i) {

        // Get a weak reference to the cached sparsity pattern
        WeakRef& wref = i->second;

        // Check if the pattern still exists
        if (wref.alive()) {

          // Get an owning reference to the cached pattern
          Sparsity ref = shared_cast<Sparsity>(wref.shared());

          // Check if the pattern matches
          if (ref.is_equal(nrow, ncol, colind, row)) {

            // Found match!
            own(ref.get());
            return;

          } else { // There is a hash rowision (unlikely, but possible)
            // Leave the pattern alone, continue to the next matching pattern
            continue;
          }
        } else {

          // Check if one of the other cache entries indeed has a matching sparsity
          CachingMap::iterator j=i;
          j++; // Start at the next matching key
          for (; j!=eq.second; ++j) {
            if (j->second.alive()) {

              // Recover cached sparsity
              Sparsity ref = shared_cast<Sparsity>(j->second.shared());

              // Match found if sparsity matches
              if (ref.is_equal(nrow, ncol, colind, row)) {
                own(ref.get());
                return;
              }
            }
          }

          // The cached entry has been deleted, create a new one
          own(new SparsityInternal(nrow, ncol, colind, row));

          // Cache this pattern
          wref = *this;

          // Return
          return;
        }
      }
    }

    // No matching sparsity pattern could be found, create a new one
    own(new SparsityInternal(nrow, ncol, colind, row));

    // Cache this pattern
    cache.insert(std::pair<std::size_t, WeakRef>(h, *this));

    // Garbage collection (currently only supported for unordered_multimap)
    casadi_int bucket_count_after = cache.bucket_count();

    // We we increased the number of buckets, take time to garbage-collect deleted references
    if (bucket_count_before!=bucket_count_after) {
      CachingMap::const_iterator i=cache.begin();
      while (i!=cache.end()) {
        if (!i->second.alive()) {
          i = cache.erase(i);
        } else {
          i++;
        }
      }
    }
  }

  Sparsity Sparsity::tril(const Sparsity& x, bool includeDiagonal) {
    return x->_tril(includeDiagonal);
  }

  Sparsity Sparsity::triu(const Sparsity& x, bool includeDiagonal) {
    return x->_triu(includeDiagonal);
  }

  std::vector<casadi_int> Sparsity::get_lower() const {
    return (*this)->get_lower();
  }

  std::vector<casadi_int> Sparsity::get_upper() const {
    return (*this)->get_upper();
  }


  std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& colind,
                            const std::vector<casadi_int>& row) {
    return hash_sparsity(nrow, ncol, get_ptr(colind), get_ptr(row));
  }

  std::size_t hash_sparsity(casadi_int nrow, casadi_int ncol,
      const casadi_int* colind, const casadi_int* row) {
    // Condense the sparsity pattern to a single, deterministric number
    std::size_t ret=0;
    hash_combine(ret, nrow);
    hash_combine(ret, ncol);
    hash_combine(ret, colind, ncol+1);
    hash_combine(ret, row, colind[ncol]);
    return ret;
  }

  Sparsity Sparsity::dense(casadi_int nrow, casadi_int ncol) {
    casadi_assert_dev(nrow>=0);
    casadi_assert_dev(ncol>=0);
    // Column offset
    std::vector<casadi_int> colind(ncol+1);
    for (casadi_int cc=0; cc<ncol+1; ++cc) colind[cc] = cc*nrow;

    // Row
    std::vector<casadi_int> row(ncol*nrow);
    for (casadi_int cc=0; cc<ncol; ++cc)
      for (casadi_int rr=0; rr<nrow; ++rr)
        row[rr+cc*nrow] = rr;

    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::upper(casadi_int n) {
    casadi_assert(n>=0, "Sparsity::upper expects a positive integer as argument");
    casadi_int nrow=n, ncol=n;
    std::vector<casadi_int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);

    // Loop over columns
    colind.push_back(0);
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Loop over rows for the upper triangular half
      for (casadi_int rr=0; rr<=cc; ++rr) {
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::lower(casadi_int n) {
    casadi_assert(n>=0, "Sparsity::lower expects a positive integer as argument");
    casadi_int nrow=n, ncol=n;
    std::vector<casadi_int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);

    // Loop over columns
    colind.push_back(0);
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Loop over rows for the lower triangular half
      for (casadi_int rr=cc; rr<nrow; ++rr) {
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::band(casadi_int n, casadi_int p) {
    casadi_assert(n>=0, "Sparsity::band expects a positive integer as argument");
    casadi_assert((p<0? -p : p)<n,
                          "Sparsity::band: position of band schould be smaller then size argument");

    casadi_int nc = n-(p<0? -p : p);

    std::vector< casadi_int >          row(nc);

    casadi_int offset = max(p, casadi_int(0));
    for (casadi_int i=0;i<nc;i++) {
      row[i]=i+offset;
    }

    std::vector< casadi_int >          colind(n+1);

    offset = min(p, casadi_int(0));
    for (casadi_int i=0;i<n+1;i++) {
      colind[i]=max(min(i+offset, nc), casadi_int(0));
    }

    return Sparsity(n, n, colind, row);

  }

  Sparsity Sparsity::banded(casadi_int n, casadi_int p) {
    // This is not an efficient implementation
    Sparsity ret = Sparsity(n, n);
    for (casadi_int i=-p;i<=p;++i) {
      ret = ret + Sparsity::band(n, i);
    }
    return ret;
  }

  Sparsity Sparsity::unit(casadi_int n, casadi_int el) {
    std::vector<casadi_int> row(1, el), colind(2);
    colind[0] = 0;
    colind[1] = 1;
    return Sparsity(n, 1, colind, row);
  }

  Sparsity Sparsity::rowcol(const std::vector<casadi_int>& row, const std::vector<casadi_int>& col,
                            casadi_int nrow, casadi_int ncol) {
    std::vector<casadi_int> all_rows, all_cols;
    all_rows.reserve(row.size()*col.size());
    all_cols.reserve(row.size()*col.size());
    for (std::vector<casadi_int>::const_iterator c_it=col.begin(); c_it!=col.end(); ++c_it) {
      casadi_assert(*c_it>=0 && *c_it<ncol, "Sparsity::rowcol: Column index out of bounds");
      for (std::vector<casadi_int>::const_iterator r_it=row.begin(); r_it!=row.end(); ++r_it) {
        casadi_assert(*r_it>=0 && *r_it<nrow, "Sparsity::rowcol: Row index out of bounds");
        all_rows.push_back(*r_it);
        all_cols.push_back(*c_it);
      }
    }
    return Sparsity::triplet(nrow, ncol, all_rows, all_cols);
  }

  Sparsity Sparsity::triplet(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& row,
                             const std::vector<casadi_int>& col, std::vector<casadi_int>& mapping,
                             bool invert_mapping) {
    // Assert dimensions
    casadi_assert_dev(nrow>=0);
    casadi_assert_dev(ncol>=0);
    casadi_assert(col.size()==row.size(), "inconsistent lengths");

    // Create the return sparsity pattern and access vectors
    std::vector<casadi_int> r_colind(ncol+1, 0);
    std::vector<casadi_int> r_row;
    r_row.reserve(row.size());

    // Consistency check and check if elements are already perfectly ordered with no duplicates
    casadi_int last_col=-1, last_row=-1;
    bool perfectly_ordered=true;
    for (casadi_int k=0; k<col.size(); ++k) {
      // Consistency check
      casadi_assert(col[k]>=0 && col[k]<ncol, "Column index out of bounds");
      casadi_assert(row[k]>=0 && row[k]<nrow, "Row index out of bounds");

      // Check if ordering is already perfect
      perfectly_ordered = perfectly_ordered && (col[k]<last_col ||
                                                (col[k]==last_col && row[k]<=last_row));
      last_col = col[k];
      last_row = row[k];
    }

    // Quick return if perfectly ordered
    if (perfectly_ordered) {
      // Save rows
      r_row.resize(row.size());
      copy(row.begin(), row.end(), r_row.begin());

      // Find offset index
      casadi_int el=0;
      for (casadi_int i=0; i<ncol; ++i) {
        while (el<col.size() && col[el]==i) el++;
        r_colind[i+1] = el;
      }

      // Identity mapping
      mapping.resize(row.size());
      for (casadi_int k=0; k<row.size(); ++k) mapping[k] = k;

      // Quick return
      return Sparsity(nrow, ncol, r_colind, r_row);
    }

    // Reuse data
    std::vector<casadi_int>& mapping1 = invert_mapping ? r_row : mapping;
    std::vector<casadi_int>& mapping2 = invert_mapping ? mapping : r_row;

    // Make sure that enough memory is allocated to use as a work vector
    mapping1.reserve(std::max(nrow+1, static_cast<casadi_int>(col.size())));

    // Number of elements in each row
    std::vector<casadi_int>& rowcount = mapping1; // reuse memory
    rowcount.resize(nrow+1);
    fill(rowcount.begin(), rowcount.end(), 0);
    for (std::vector<casadi_int>::const_iterator it=row.begin(); it!=row.end(); ++it) {
      rowcount[*it+1]++;
    }

    // Cumsum to get index offset for each row
    for (casadi_int i=0; i<nrow; ++i) {
      rowcount[i+1] += rowcount[i];
    }

    // New row for each old row
    mapping2.resize(row.size());
    for (casadi_int k=0; k<row.size(); ++k) {
      mapping2[rowcount[row[k]]++] = k;
    }

    // Number of elements in each col
    // reuse memory, r_colind is already the right size
    std::vector<casadi_int>& colcount = r_colind;
                                           // and is filled with zeros
    for (std::vector<casadi_int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it) {
      colcount[col[*it]+1]++;
    }

    // Cumsum to get index offset for each col
    for (casadi_int i=0; i<ncol; ++i) {
      colcount[i+1] += colcount[i];
    }

    // New col for each old col
    mapping1.resize(col.size());
    for (std::vector<casadi_int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it) {
      mapping1[colcount[col[*it]]++] = *it;
    }

    // Current element in the return matrix
    casadi_int r_el = 0;
    r_row.resize(col.size());

    // Current nonzero
    std::vector<casadi_int>::const_iterator it=mapping1.begin();

    // Loop over columns
    r_colind[0] = 0;
    for (casadi_int i=0; i<ncol; ++i) {

      // Previous row (to detect duplicates)
      casadi_int j_prev = -1;

      // Loop over nonzero elements of the col
      while (it!=mapping1.end() && col[*it]==i) {

        // Get the element
        casadi_int el = *it;
        it++;

        // Get the row
        casadi_int j = row[el];

        // If not a duplicate, save to return matrix
        if (j!=j_prev)
          r_row[r_el++] = j;

        if (invert_mapping) {
          // Save to the inverse mapping
          mapping2[el] = r_el-1;
        } else {
          // If not a duplicate, save to the mapping vector
          if (j!=j_prev)
            mapping1[r_el-1] = el;
        }

        // Save row
        j_prev = j;
      }

      // Update col offset
      r_colind[i+1] = r_el;
    }

    // Resize the row vector
    r_row.resize(r_el);

    // Resize mapping matrix
    if (!invert_mapping) {
      mapping1.resize(r_el);
    }

    return Sparsity(nrow, ncol, r_colind, r_row);
  }

  Sparsity Sparsity::triplet(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& row,
                             const std::vector<casadi_int>& col) {
    std::vector<casadi_int> mapping;
    return Sparsity::triplet(nrow, ncol, row, col, mapping, false);
  }

  Sparsity Sparsity::nonzeros(casadi_int nrow, casadi_int ncol,
      const std::vector<casadi_int>& nz, bool ind1) {
    casadi_assert(nrow>0, "nrow must be >0.");
    std::vector<casadi_int> row(nz.size());
    std::vector<casadi_int> col(nz.size());
    for (casadi_int i=0;i<nz.size();++i) {
      casadi_int k = nz[i];
      k-= ind1;
      row[i] = k % nrow;
      col[i] = k / nrow;
    }
    return triplet(nrow, ncol, row, col);
  }

  bool Sparsity::is_singular() const {
    casadi_assert(is_square(),
      "is_singular: only defined for square matrices, but got " + dim());
    return sprank(*this)!=size2();
  }

  std::vector<casadi_int> Sparsity::compress() const {
    return (*this)->sp();
  }

  Sparsity::operator const std::vector<casadi_int>&() const {
    return (*this)->sp();
  }

  Sparsity::operator SparsityStruct() const {
    const casadi_int* sp = *this;
    casadi_int nrow = sp[0], ncol = sp[1];
    const casadi_int* colind = sp+2, *row = sp+2+ncol+1;
    return SparsityStruct{nrow, ncol, colind, row};
  }

  Sparsity Sparsity::compressed(const std::vector<casadi_int>& v, bool order_rows) {
    // Check consistency
    casadi_assert_dev(v.size() >= 2);
    casadi_int nrow = v[0];
    casadi_int ncol = v[1];
    casadi_assert_dev(v.size() >= 2 + ncol+1);
    casadi_int nnz = v[2 + ncol];
    bool dense = v.size() == 2 + ncol+1 && nrow*ncol==nnz;
    bool sparse = v.size() == 2 + ncol+1 + nnz;
    casadi_assert_dev(dense || sparse);

    // Call array version
    return compressed(&v.front(), order_rows);
  }

  Sparsity Sparsity::compressed(const casadi_int* v, bool order_rows) {
    casadi_assert_dev(v!=nullptr);

    // Get sparsity pattern
    casadi_int nrow = v[0];
    casadi_int ncol = v[1];
    const casadi_int *colind = v+2;
    if (colind[0]==1) {
      // Dense matrix
      return Sparsity::dense(nrow, ncol);
    }
    casadi_int nnz = colind[ncol];
    if (nrow*ncol == nnz) {
      // Dense matrix
      return Sparsity::dense(nrow, ncol);
    } else {
      // Sparse matrix
      const casadi_int *row = v + 2 + ncol+1;
      return Sparsity(nrow, ncol,
                      vector<casadi_int>(colind, colind+ncol+1),
                      vector<casadi_int>(row, row+nnz), order_rows);
    }
  }

  casadi_int Sparsity::bw_upper() const {
    return (*this)->bw_upper();
  }

  casadi_int Sparsity::bw_lower() const {
    return (*this)->bw_lower();
  }

  Sparsity Sparsity::horzcat(const std::vector<Sparsity> & sp) {
    // Quick return if possible
    if (sp.empty()) return Sparsity(0, 0);
    if (sp.size()==1) return sp.front();

    // Count total nnz
    casadi_int nnz_total = 0;
    for (casadi_int i=0; i<sp.size(); ++i) nnz_total += sp[i].nnz();

    // Construct from vectors (triplet format)
    vector<casadi_int> ret_row, ret_col;
    ret_row.reserve(nnz_total);
    ret_col.reserve(nnz_total);
    casadi_int ret_ncol = 0;
    casadi_int ret_nrow = 0;
    for (casadi_int i=0; i<sp.size() && ret_nrow==0; ++i)
      ret_nrow = sp[i].size1();

    // Append all patterns
    for (vector<Sparsity>::const_iterator i=sp.begin(); i!=sp.end(); ++i) {
      // Get sparsity pattern
      casadi_int sp_nrow = i->size1();
      casadi_int sp_ncol = i->size2();
      const casadi_int* sp_colind = i->colind();
      const casadi_int* sp_row = i->row();
      casadi_assert(sp_nrow==ret_nrow || sp_nrow==0,
                            "Sparsity::horzcat: Mismatching number of rows");

      // Add entries to pattern
      for (casadi_int cc=0; cc<sp_ncol; ++cc) {
        for (casadi_int k=sp_colind[cc]; k<sp_colind[cc+1]; ++k) {
          ret_row.push_back(sp_row[k]);
          ret_col.push_back(cc + ret_ncol);
        }
      }

      // Update offset
      ret_ncol += sp_ncol;
    }
    return Sparsity::triplet(ret_nrow, ret_ncol, ret_row, ret_col);
  }

  Sparsity Sparsity::kron(const Sparsity& a, const Sparsity& b) {
    casadi_int a_ncol = a.size2();
    casadi_int b_ncol = b.size2();
    casadi_int a_nrow = a.size1();
    casadi_int b_nrow = b.size1();
    if (a.is_dense() && b.is_dense()) return Sparsity::dense(a_nrow*b_nrow, a_ncol*b_ncol);

    const casadi_int* a_colind = a.colind();
    const casadi_int* a_row = a.row();
    const casadi_int* b_colind = b.colind();
    const casadi_int* b_row = b.row();

    std::vector<casadi_int> r_colind(a_ncol*b_ncol+1, 0);
    std::vector<casadi_int> r_row(a.nnz()*b.nnz());

    casadi_int* r_colind_ptr = get_ptr(r_colind);
    casadi_int* r_row_ptr = get_ptr(r_row);

    casadi_int i=0;
    casadi_int j=0;
    // Loop over the columns
    for (casadi_int a_cc=0; a_cc<a_ncol; ++a_cc) {
      casadi_int a_start = a_colind[a_cc];
      casadi_int a_stop  = a_colind[a_cc+1];
      // Loop over the columns
      for (casadi_int b_cc=0; b_cc<b_ncol; ++b_cc) {
        casadi_int b_start = b_colind[b_cc];
        casadi_int b_stop  = b_colind[b_cc+1];
        // Loop over existing nonzeros
        for (casadi_int a_el=a_start; a_el<a_stop; ++a_el) {
          casadi_int a_r = a_row[a_el];
          // Loop over existing nonzeros
          for (casadi_int b_el=b_start; b_el<b_stop; ++b_el) {
            casadi_int b_r = b_row[b_el];
            r_row_ptr[i++] = a_r*b_nrow+b_r;
          }
        }
        j+=1;
        r_colind_ptr[j] = r_colind_ptr[j-1] + (b_stop-b_start)*(a_stop-a_start);
      }
    }
    return Sparsity(a_nrow*b_nrow, a_ncol*b_ncol, r_colind, r_row);
  }

  Sparsity Sparsity::vertcat(const std::vector<Sparsity> & sp) {
    // Quick return if possible
    if (sp.empty()) return Sparsity(0, 0);
    if (sp.size()==1) return sp.front();

    // Count total nnz
    casadi_int nnz_total = 0;
    for (casadi_int i=0; i<sp.size(); ++i) nnz_total += sp[i].nnz();

    // Construct from vectors (triplet format)
    vector<casadi_int> ret_row, ret_col;
    ret_row.reserve(nnz_total);
    ret_col.reserve(nnz_total);
    casadi_int ret_nrow = 0;
    casadi_int ret_ncol = 0;
    for (casadi_int i=0; i<sp.size() && ret_ncol==0; ++i)
      ret_ncol = sp[i].size2();

    // Append all patterns
    for (vector<Sparsity>::const_iterator i=sp.begin(); i!=sp.end(); ++i) {
      // Get sparsity pattern
      casadi_int sp_nrow = i->size1();
      casadi_int sp_ncol = i->size2();
      const casadi_int* sp_colind = i->colind();
      const casadi_int* sp_row = i->row();
      casadi_assert(sp_ncol==ret_ncol || sp_ncol==0,
                            "Sparsity::vertcat: Mismatching number of columns");

      // Add entries to pattern
      for (casadi_int cc=0; cc<sp_ncol; ++cc) {
        for (casadi_int k=sp_colind[cc]; k<sp_colind[cc+1]; ++k) {
          ret_row.push_back(sp_row[k] + ret_nrow);
          ret_col.push_back(cc);
        }
      }

      // Update offset
      ret_nrow += sp_nrow;
    }
    return Sparsity::triplet(ret_nrow, ret_ncol, ret_row, ret_col);
  }

  Sparsity Sparsity::diagcat(const std::vector< Sparsity > &v) {
    casadi_int n = 0;
    casadi_int m = 0;

    std::vector<casadi_int> colind(1, 0);
    std::vector<casadi_int> row;

    casadi_int nz = 0;
    for (casadi_int i=0;i<v.size();++i) {
      const casadi_int* colind_ = v[i].colind();
      casadi_int ncol = v[i].size2();
      const casadi_int* row_ = v[i].row();
      casadi_int sz = v[i].nnz();
      for (casadi_int k=1; k<ncol+1; ++k) {
        colind.push_back(colind_[k]+nz);
      }
      for (casadi_int k=0; k<sz; ++k) {
        row.push_back(row_[k]+m);
      }
      n+= v[i].size2();
      m+= v[i].size1();
      nz+= v[i].nnz();
    }

    return Sparsity(m, n, colind, row);
  }

  std::vector<Sparsity> Sparsity::horzsplit(const Sparsity& x,
      const std::vector<casadi_int>& offset) {
    // Consistency check
    casadi_assert_dev(offset.size()>=1);
    casadi_assert_dev(offset.front()==0);
    casadi_assert(offset.back()==x.size2(),
                          "horzsplit(Sparsity, std::vector<casadi_int>): Last elements of offset "
                          "(" + str(offset.back()) + ") must equal the number of columns "
                          "(" + str(x.size2()) + ")");
    casadi_assert_dev(is_monotone(offset));

    // Number of outputs
    casadi_int n = offset.size()-1;

    // Get the sparsity of the input
    const casadi_int* colind_x = x.colind();
    const casadi_int* row_x = x.row();

    // Allocate result
    std::vector<Sparsity> ret;
    ret.reserve(n);

    // Sparsity pattern as CCS vectors
    vector<casadi_int> colind, row;
    casadi_int ncol, nrow = x.size1();

    // Get the sparsity patterns of the outputs
    for (casadi_int i=0; i<n; ++i) {
      casadi_int first_col = offset[i];
      casadi_int last_col = offset[i+1];
      ncol = last_col - first_col;

      // Construct the sparsity pattern
      colind.resize(ncol+1);
      copy(colind_x+first_col, colind_x+last_col+1, colind.begin());
      for (vector<casadi_int>::iterator it=colind.begin()+1; it!=colind.end(); ++it)
        *it -= colind[0];
      colind[0] = 0;
      row.resize(colind.back());
      copy(row_x+colind_x[first_col], row_x+colind_x[last_col], row.begin());

      // Append to the list
      ret.push_back(Sparsity(nrow, ncol, colind, row));
    }

    // Return (RVO)
    return ret;
  }

  std::vector<Sparsity> Sparsity::vertsplit(const Sparsity& x,
      const std::vector<casadi_int>& offset) {
    std::vector<Sparsity> ret = horzsplit(x.T(), offset);
    for (std::vector<Sparsity>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = it->T();
    }
    return ret;
  }

  Sparsity Sparsity::blockcat(const std::vector< std::vector< Sparsity > > &v) {
    std::vector< Sparsity > ret;
    for (casadi_int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  std::vector<Sparsity> Sparsity::diagsplit(const Sparsity& x,
                                            const std::vector<casadi_int>& offset1,
                                            const std::vector<casadi_int>& offset2) {
    // Consistency check
    casadi_assert_dev(offset1.size()>=1);
    casadi_assert_dev(offset1.front()==0);
    casadi_assert(offset1.back()==x.size1(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset1 "
                          "(" + str(offset1.back()) + ") must equal the number of rows "
                          "(" + str(x.size1()) + ")");
    casadi_assert(offset2.back()==x.size2(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset2 "
                          "(" + str(offset2.back()) + ") must equal the number of rows "
                          "(" + str(x.size2()) + ")");
    casadi_assert_dev(is_monotone(offset1));
    casadi_assert_dev(is_monotone(offset2));
    casadi_assert_dev(offset1.size()==offset2.size());

    // Number of outputs
    casadi_int n = offset1.size()-1;

    // Return value
    std::vector<Sparsity> ret;

    // Caveat: this is a very silly implementation
    IM x2 = IM::zeros(x);

    for (casadi_int i=0; i<n; ++i) {
      ret.push_back(x2(Slice(offset1[i], offset1[i+1]),
                       Slice(offset2[i], offset2[i+1])).sparsity());
    }

    return ret;
  }

  Sparsity Sparsity::sum2(const Sparsity &x) {
    return mtimes(x, Sparsity::dense(x.size2(), 1));

  }
  Sparsity Sparsity::sum1(const Sparsity &x) {
    return mtimes(Sparsity::dense(1, x.size1()), x);
  }

  casadi_int Sparsity::sprank(const Sparsity& x) {
    std::vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    x.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

  Sparsity::operator const casadi_int*() const {
    return &(*this)->sp().front();
  }

  casadi_int Sparsity::norm_0_mul(const Sparsity& x, const Sparsity& A) {
    // Implementation borrowed from Scipy's sparsetools/csr.h
    casadi_assert(A.size1()==x.size2(), "Dimension error. Got " + x.dim()
                          + " times " + A.dim() + ".");

    casadi_int n_row = A.size2();
    casadi_int n_col = x.size1();

    // Allocate work vectors
    std::vector<bool> Bwork(n_col);
    std::vector<casadi_int> Iwork(n_row+1+n_col);

    const casadi_int* Aj = A.row();
    const casadi_int* Ap = A.colind();
    const casadi_int* Bj = x.row();
    const casadi_int* Bp = x.colind();
    casadi_int *Cp = get_ptr(Iwork);
    casadi_int *mask = Cp+n_row+1;

    // Pass 1
    // method that uses O(n) temp storage
    std::fill(mask, mask+n_col, -1);

    Cp[0] = 0;
    casadi_int nnz = 0;
    for (casadi_int i = 0; i < n_row; i++) {
      casadi_int row_nnz = 0;
      for (casadi_int jj = Ap[i]; jj < Ap[i+1]; jj++) {
        casadi_int j = Aj[jj];
        for (casadi_int kk = Bp[j]; kk < Bp[j+1]; kk++) {
          casadi_int k = Bj[kk];
          if (mask[k] != i) {
            mask[k] = i;
            row_nnz++;
          }
        }
      }
      casadi_int next_nnz = nnz + row_nnz;
      nnz = next_nnz;
      Cp[i+1] = nnz;
    }

    // Pass 2
    casadi_int *next = get_ptr(Iwork) + n_row+1;
    std::fill(next, next+n_col, -1);
    std::vector<bool> & sums = Bwork;
    std::fill(sums.begin(), sums.end(), false);
    nnz = 0;
    Cp[0] = 0;
    for (casadi_int i = 0; i < n_row; i++) {
      casadi_int head   = -2;
      casadi_int length =  0;
      casadi_int jj_start = Ap[i];
      casadi_int jj_end   = Ap[i+1];
      for (casadi_int jj = jj_start; jj < jj_end; jj++) {
        casadi_int j = Aj[jj];
        casadi_int kk_start = Bp[j];
        casadi_int kk_end   = Bp[j+1];
        for (casadi_int kk = kk_start; kk < kk_end; kk++) {
          casadi_int k = Bj[kk];
          sums[k] = true;
          if (next[k] == -1) {
            next[k] = head;
            head  = k;
            length++;
          }
        }
      }
      for (casadi_int jj = 0; jj < length; jj++) {
        if (sums[head]) {
          nnz++;
        }
        casadi_int temp = head;
        head = next[head];
        next[temp] = -1; //clear arrays
        sums[temp] =  0;
      }
      Cp[i+1] = nnz;
    }
    return nnz;
  }

  void Sparsity::mul_sparsityF(const bvec_t* x, const Sparsity& x_sp,
                               const bvec_t* y, const Sparsity& y_sp,
                               bvec_t* z, const Sparsity& z_sp,
                               bvec_t* w) {
    // Assert dimensions
    casadi_assert(z_sp.size1()==x_sp.size1() && x_sp.size2()==y_sp.size1()
                          && y_sp.size2()==z_sp.size2(),
                          "Dimension error. Got x=" + x_sp.dim() + ", y=" + y_sp.dim()
                          + " and z=" + z_sp.dim() + ".");

    // Direct access to the arrays
    const casadi_int* y_colind = y_sp.colind();
    const casadi_int* y_row = y_sp.row();
    const casadi_int* x_colind = x_sp.colind();
    const casadi_int* x_row = x_sp.row();
    const casadi_int* z_colind = z_sp.colind();
    const casadi_int* z_row = z_sp.row();

    // Loop over the columns of y and z
    casadi_int ncol = z_sp.size2();
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Get the dense column of z
      for (casadi_int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        w[z_row[kk]] = z[kk];
      }

      // Loop over the nonzeros of y
      for (casadi_int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
        casadi_int rr = y_row[kk];

        // Loop over corresponding columns of x
        bvec_t yy = y[kk];
        for (casadi_int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
          w[x_row[kk1]] |= x[kk1] | yy;
        }
      }

      // Get the sparse column of z
      for (casadi_int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        z[kk] = w[z_row[kk]];
      }
    }
  }

  void Sparsity::mul_sparsityR(bvec_t* x, const Sparsity& x_sp,
                               bvec_t* y, const Sparsity& y_sp,
                               bvec_t* z, const Sparsity& z_sp,
                               bvec_t* w) {
    // Assert dimensions
    casadi_assert(z_sp.size1()==x_sp.size1() && x_sp.size2()==y_sp.size1()
                          && y_sp.size2()==z_sp.size2(),
                          "Dimension error. Got x=" + x_sp.dim() + ", y=" + y_sp.dim()
                          + " and z=" + z_sp.dim() + ".");

    // Direct access to the arrays
    const casadi_int* y_colind = y_sp.colind();
    const casadi_int* y_row = y_sp.row();
    const casadi_int* x_colind = x_sp.colind();
    const casadi_int* x_row = x_sp.row();
    const casadi_int* z_colind = z_sp.colind();
    const casadi_int* z_row = z_sp.row();

    // Loop over the columns of y and z
    casadi_int ncol = z_sp.size2();
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Get the dense column of z
      for (casadi_int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        w[z_row[kk]] = z[kk];
      }

      // Loop over the nonzeros of y
      for (casadi_int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
        casadi_int rr = y_row[kk];

        // Loop over corresponding columns of x
        bvec_t yy = 0;
        for (casadi_int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
          yy |= w[x_row[kk1]];
          x[kk1] |= w[x_row[kk1]];
        }
        y[kk] |= yy;
      }

      // Get the sparse column of z
      for (casadi_int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        z[kk] = w[z_row[kk]];
      }
    }
  }

  Dict Sparsity::info() const {
    if (is_null()) return Dict();
    return {{"nrow", size1()}, {"ncol", size2()}, {"colind", get_colind()}, {"row", get_row()}};
  }

  Sparsity Sparsity::from_info(const Dict& info) {
    auto it = info.find("nrow");
    if (it==info.end()) return Sparsity();
    casadi_int nrow = info.at("nrow");
    casadi_int ncol = info.at("ncol");
    std::vector<casadi_int> row, colind;
    if (info.at("row").is_int_vector()) {
      row = info.at("row");
    } else {
      row.push_back(info.at("row"));
    }
    if (info.at("colind").is_int_vector()) {
      colind = info.at("colind");
    } else {
      colind.push_back(info.at("colind"));
    }
    return Sparsity(nrow, ncol, colind, row);
  }

  std::set<std::string> Sparsity::file_formats = {"mtx"};

  std::string Sparsity::file_format(const std::string& filename, const std::string& format_hint) {
    if (format_hint=="") {
      std::string extension = filename.substr(filename.rfind(".")+1);
      auto it = file_formats.find(extension);
      casadi_assert(it!=file_formats.end(),
        "Extension '" + extension + "' not recognised. "
        "Valid options: " + str(file_formats) + ".");
      return extension;
    } else {
      auto it = file_formats.find(format_hint);
      casadi_assert(it!=file_formats.end(),
        "File format hint '" + format_hint + "' not recognised. "
        "Valid options: " + str(file_formats) + ".");
      return format_hint;
    }

  }
  void Sparsity::to_file(const std::string& filename, const std::string& format_hint) const {
    std::string format = file_format(filename, format_hint);
    std::ofstream out(filename);
    if (format=="mtx") {
      out << std::scientific << std::setprecision(15);
      out << "%%MatrixMarket matrix coordinate pattern general" << std::endl;
      out << size1() << " " << size2() << " " << nnz() << std::endl;
      std::vector<casadi_int> row = get_row();
      std::vector<casadi_int> col = get_col();

      for (casadi_int k=0;k<row.size();++k) {
        out << row[k]+1 << " " << col[k]+1 << std::endl;
      }
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }

  Sparsity Sparsity::from_file(const std::string& filename, const std::string& format_hint) {
    std::string format = file_format(filename, format_hint);
    std::ifstream in(filename);
    if (format=="mtx") {
      std::string line;
      std::vector<casadi_int> row, col;
      casadi_int size1, size2, nnz;
      int line_num = 0;
      while (std::getline(in, line)) {
        if (line_num==0) {
          casadi_assert(line=="%%MatrixMarket matrix coordinate pattern general", "Wrong header");
          line_num = 1;
        } else if (line_num==1) {
          std::stringstream stream(line);
          stream >> size1;
          stream >> size2;
          stream >> nnz;
          row.reserve(nnz);
          col.reserve(nnz);
          line_num = 2;
        } else {
          std::stringstream stream(line);
          casadi_int r, c;
          stream >> r;
          stream >> c;
          row.push_back(r-1);
          col.push_back(c-1);
        }
      }
      return triplet(size1, size2, row, col);
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }

  Sparsity Sparsity::kkt(const Sparsity& H, const Sparsity& J,
                         bool with_x_diag, bool with_lam_g_diag) {
    // Consistency check
    casadi_assert(H.is_square(), "H must be square");
    casadi_assert(H.size1() == J.size2(), "Dimension mismatch");

    // Add diagonal to H recursively
    if (with_x_diag) return kkt(H + diag(H.size()), J, false, with_lam_g_diag);

    // Lower right entry
    int ng = J.size1();
    Sparsity B = with_lam_g_diag ? diag(ng, ng) : Sparsity(ng, ng);

    // Concatenate
    return blockcat({{H, J.T()}, {J, B}});
  }

  void Sparsity::serialize(std::ostream &stream) const {
    casadi_int size1=this->size1(), size2=this->size2(), nnz=this->nnz();
    const casadi_int *colind = this->colind(), *row = this->row();
    stream << "sp";
    stream << size1 << "x" << size2;
    stream << "n" << nnz;
    for (int i=0; i<nnz; ++i)
      stream << ":" << row[i];
    for (int i=0; i<size2+1; ++i)
      stream << ":" << colind[i];
    stream << "s";
  }

  Sparsity Sparsity::deserialize(std::istream &stream) {
    char ch;
    stream >> ch;
    stream >> ch;
    casadi_int nrow, ncol, nnz;
    stream >> nrow; stream >> ch;
    stream >> ncol; stream >> ch;
    stream >> nnz;
    std::vector<casadi_int> row(nnz), colind(ncol+1);
    for (int i=0; i<nnz; ++i) {
      stream >> ch;
      stream >> row[i];
    }
    for (int i=0; i<ncol+1; ++i) {
      stream >> ch;
      stream >> colind[i];
    }
    stream >> ch;
    return Sparsity(nrow, ncol, colind, row);
  }

  std::string Sparsity::serialize() const {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }

  Sparsity Sparsity::deserialize(const std::string& s) {
    std::stringstream ss;
    ss << s;
    return deserialize(ss);
  }
} // namespace casadi
