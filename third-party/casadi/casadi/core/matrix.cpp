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

#define CASADI_MATRIX_CPP
#include "sx_function.hpp"
#include "sx_node.hpp"
#include "linsol.hpp"
#include "expm.hpp"
#include <chrono>

using namespace std;

namespace casadi {
  template<typename Scalar>
  void Matrix<Scalar>::set_precision(casadi_int precision) { stream_precision_ = precision; }

  template<typename Scalar>
  void Matrix<Scalar>::set_width(casadi_int width) { stream_width_ = width; }

  template<typename Scalar>
  void Matrix<Scalar>::set_scientific(bool scientific) { stream_scientific_ = scientific; }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::stream_precision_ = 6;
  template<typename Scalar>
  casadi_int Matrix<Scalar>::stream_width_ = 0;
  template<typename Scalar>
  bool Matrix<Scalar>::stream_scientific_ = false;

  template<typename Scalar>
  std::default_random_engine Matrix<Scalar>::rng_(
    // Seed with current time
    std::chrono::system_clock::now().time_since_epoch().count());

  template<typename Scalar>
  void Matrix<Scalar>::rng(casadi_int seed) {
    rng_.seed(seed);
  }

  template<typename Scalar>
  bool Matrix<Scalar>::__nonzero__() const {
    if (numel()!=1) {
      casadi_error("Only scalar Matrix could have a truth value, but you "
                   "provided a shape" + dim());
    }
    return nonzeros().at(0)!=0;
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Slice& rr, const Slice& cc) const {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2())) {
      casadi_int k = sparsity().get_nz(rr.scalar(size1()), cc.scalar(size2()));
      if (k>=0) {
        m = nonzeros().at(k);
      } else {
        m = Matrix<Scalar>(1, 1);
      }
      return;
    }

    // Fall back on IM-IM
    get(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Slice& rr, const Matrix<casadi_int>& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Slice& cc) const {
    // Fall back on IM-IM
    get(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1,
                             const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const {
    // Scalar
    if (rr.is_scalar(true) && cc.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1), to_slice(cc, ind1));
    }

    // Make sure dense vectors
    casadi_assert(rr.is_dense() && rr.is_vector(),
                          "Marix::get: First index must be a dense vector");
    casadi_assert(cc.is_dense() && cc.is_vector(),
                          "Marix::get: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), cc.nonzeros(), mapping, ind1);

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(sp);
    for (casadi_int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Slice& rr) const {
    // Scalar
    if (rr.is_scalar(numel())) {
      casadi_int r = rr.scalar(numel());
      casadi_int k = sparsity().get_nz(r % size1(), r / size1());
      if (k>=0) {
        m = nonzeros().at(k);
      } else {
        m = Matrix<Scalar>(1, 1);
      }
      return;
    }

    // Fall back on IM
    get(m, ind1, rr.all(numel(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr) const {
    // Scalar
    if (rr.is_scalar(true)) {
      return get(m, ind1, to_slice(rr, ind1));
    }

    // If the indexed matrix is dense, use nonzero indexing
    if (is_dense()) {
      return get_nz(m, ind1, rr);
    }

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), rr.sparsity(), mapping, ind1);

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && rr.is_row()) || (is_row() && rr.is_column());

    // Copy nonzeros
    m = Matrix<Scalar>::zeros(tr ? sp.T() : sp);
    for (casadi_int k=0; k<mapping.size(); ++k) m->at(k) = nonzeros().at(mapping[k]);
  }

  template<typename Scalar>
  void Matrix<Scalar>::get(Matrix<Scalar>& m, bool ind1, const Sparsity& sp) const {
    casadi_assert(size()==sp.size(),
                          "Shape mismatch. This matrix has shape "
                          + str(size()) + ", but supplied sparsity index has shape "
                          + str(sp.size()) + ".");
    m = project(*this, sp);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                             const Slice& rr, const Slice& cc) {
    // Both are scalar
    if (rr.is_scalar(size1()) && cc.is_scalar(size2()) && m.is_dense()) {
      casadi_int oldsize = sparsity_.nnz();
      casadi_int ind = sparsity_.add_nz(rr.scalar(size1()), cc.scalar(size2()));
      if (oldsize == sparsity_.nnz()) {
        nonzeros_.at(ind) = m.scalar();
      } else {
        nonzeros_.insert(nonzeros_.begin()+ind, m.scalar());
      }
      return;
    }

    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                             const Slice& rr, const Matrix<casadi_int>& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr, cc.all(size2(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1,
                                const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) {
    // Scalar
    if (rr.is_scalar(true) && cc.is_scalar(true) && m.is_dense()) {
      return set(m, ind1, to_slice(rr, ind1), to_slice(cc, ind1));
    }

    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return set(m, ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return set(m, ind1, rr, cc.T());
    }

    // Make sure rr and cc are dense vectors
    casadi_assert(rr.is_dense() && rr.is_column(),
                          "Matrix::set: First index not dense vector");
    casadi_assert(cc.is_dense() && cc.is_column(),
                          "Matrix::set: Second index not dense vector");

    // Assert dimensions of assigning matrix
    if (rr.size1() != m.size1() || cc.size1() != m.size2()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        return set(repmat(m, rr.size1(), cc.size1()), ind1, rr, cc);
      } else if (rr.size1() == m.size2() && cc.size1() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr, cc);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(rr.size1()) + "-by-"
                     + str(cc.size1()) + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions
    casadi_int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    casadi_assert_in_range(rr.nonzeros(), -sz1+ind1, sz1+ind1);
    casadi_assert_in_range(cc.nonzeros(), -sz2+ind1, sz2+ind1);

    // If we are assigning with something sparse, first remove existing entries
    if (!m.is_dense()) {
      erase(rr.nonzeros(), cc.nonzeros(), ind1);
    }

    // Collect all assignments
    IM el = IM::zeros(m.sparsity());
    for (casadi_int j=0; j<el.size2(); ++j) { // Loop over columns of m
      casadi_int this_j = cc->at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (casadi_int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        casadi_int i = m.row(k);
        casadi_int this_i = rr->at(i) - ind1; // Corresponding row in this
        if (this_i<0) this_i += sz1;
        el->at(k) = this_i + this_j*sz1;
      }
    }
    return set(m, false, el);
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Slice& rr) {
    // Scalar
    if (rr.is_scalar(numel()) && m.is_dense()) {
      casadi_int r = rr.scalar(numel());
      casadi_int oldsize = sparsity_.nnz();
      casadi_int ind = sparsity_.add_nz(r % size1(), r / size1());
      if (oldsize == sparsity_.nnz()) {
        nonzeros_.at(ind) = m.scalar();
      } else {
        nonzeros_.insert(nonzeros_.begin()+ind, m.scalar());
      }
      return;
    }

    // Fall back on IM
    set(m, ind1, rr.all(numel(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& rr) {
    // Scalar
    if (rr.is_scalar(true) && m.is_dense()) {
      return set(m, ind1, to_slice(rr, ind1));
    }

    // Assert dimensions of assigning matrix
    if (rr.sparsity() != m.sparsity()) {
      if (rr.size() == m.size()) {
        // Remove submatrix to be replaced
        erase(rr.nonzeros(), ind1);

        // Find the intersection between rr's and m's sparsity patterns
        Sparsity sp = rr.sparsity() * m.sparsity();

        // Project both matrices to this sparsity
        return set(project(m, sp), ind1, Matrix<casadi_int>::project(rr, sp));
      } else if (m.is_scalar()) {
        // m scalar means "set all"
        if (m.is_dense()) {
          return set(Matrix<Scalar>(rr.sparsity(), m), ind1, rr);
        } else {
          return set(Matrix<Scalar>(rr.size()), ind1, rr);
        }
      } else if (rr.size1() == m.size2() && rr.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(rr.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions of this
    casadi_int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    casadi_assert_in_range(rr.nonzeros(), -nel+ind1, nel+ind1);

    // Dense mode
    if (is_dense() && m.is_dense()) {
      return set_nz(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<casadi_int> new_row =
      sparsity().get_row(), new_col=sparsity().get_col(), nz(rr.nonzeros());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    for (std::vector<casadi_int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (ind1) (*i)--;
      if (*i<0) *i += nel;
      new_row.push_back(*i % sz1);
      new_col.push_back(*i / sz1);
    }
    Sparsity sp = Sparsity::triplet(sz1, sz2, new_row, new_col);

    // If needed, update pattern
    if (sp != sparsity()) *this = project(*this, sp);

    // Find the nonzeros corresponding to rr
    sparsity().get_nz(nz);

    // Carry out the assignments
    for (casadi_int i=0; i<nz.size(); ++i) {
      nonzeros().at(nz[i]) = m->at(i);
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set(const Matrix<Scalar>& m, bool ind1, const Sparsity& sp) {
    casadi_assert(size()==sp.size(),
                          "set(Sparsity sp): shape mismatch. This matrix has shape "
                          + str(size()) + ", but supplied sparsity index has shape "
                          + str(sp.size()) + ".");
    std::vector<casadi_int> ii = sp.find();
    if (m.is_scalar()) {
      (*this)(ii) = densify(m);
    } else {
      (*this)(ii) = densify(m(ii));
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::get_nz(Matrix<Scalar>& m, bool ind1, const Slice& kk) const {
    // Scalar
    if (kk.is_scalar(nnz())) {
      m = nonzeros().at(kk.scalar(nnz()));
      return;
    }

    // Fall back on IM
    get_nz(m, ind1, kk.all(nnz(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::get_nz(Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& kk) const {
    // Scalar
    if (kk.is_scalar(true)) {
      return get_nz(m, ind1, to_slice(kk, ind1));
    }

    // Get nonzeros of kk
    const std::vector<casadi_int>& k = kk.nonzeros();
    casadi_int sz = nnz();

    // Check bounds
    casadi_assert_in_range(k, -sz+ind1, sz+ind1);

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && kk.is_row()) || (is_row() && kk.is_column());

    // Copy nonzeros
    m = zeros(tr ? kk.sparsity().T() : kk.sparsity());
    for (casadi_int el=0; el<k.size(); ++el) {
      casadi_assert(!(ind1 && k[el]<=0), "Matlab is 1-based, but requested index "
                                                + str(k[el]) + ". Note that negative slices are"
                                                " disabled in the Matlab interface. "
                                                "Possibly you may want to use 'end'.");
      casadi_int k_el = k[el]-ind1;
      m->at(el) = nonzeros().at(k_el>=0 ? k_el : k_el+sz);
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_nz(const Matrix<Scalar>& m, bool ind1, const Slice& kk) {
    // Scalar
    if (kk.is_scalar(nnz())) {
      nonzeros().at(kk.scalar(nnz())) = m.scalar();
      return;
    }

    // Fallback on IM
    set_nz(m, ind1, kk.all(nnz(), ind1));
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_nz(const Matrix<Scalar>& m, bool ind1, const Matrix<casadi_int>& kk) {
    // Scalar
    if (kk.is_scalar(true)) {
      return set_nz(m, ind1, to_slice(kk, ind1));
    }

    // Assert dimensions of assigning matrix
    if (kk.sparsity() != m.sparsity()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        if (!m.is_dense()) return; // Nothing to set
        return set_nz(Matrix<Scalar>(kk.sparsity(), m), ind1, kk);
      } else if (kk.size() == m.size()) {
        // Project sparsity if needed
        return set_nz(project(m, kk.sparsity()), ind1, kk);
      } else if (kk.size1() == m.size2() && kk.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set_nz(m.T(), ind1, kk);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(kk.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Get nonzeros
    const std::vector<casadi_int>& k = kk.nonzeros();
    casadi_int sz = nnz();

    // Check bounds
    casadi_assert_in_range(k, -sz+ind1, sz+ind1);

    // Set nonzeros, ignoring negative indices
    for (casadi_int el=0; el<k.size(); ++el) {
      casadi_assert(!(ind1 && k[el]<=0),
        "Matlab is 1-based, but requested index " + str(k[el])
        +  ". Note that negative slices are disabled in the Matlab interface. "
           "Possibly you may want to use 'end'.");
      casadi_int k_el = k[el]-ind1;
      nonzeros().at(k_el>=0 ? k_el : k_el+sz) = m->at(el);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::densify(const Matrix<Scalar>& x) {
    return densify(x, 0);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::densify(const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& val) {
    // Check argument
    casadi_assert_dev(val.is_scalar());

    // Quick return if possible
    if (x.is_dense()) return x;

    // Get sparsity pattern
    casadi_int nrow = x.size1();
    casadi_int ncol = x.size2();
    const casadi_int* colind = x.colind();
    const casadi_int* row = x.row();
    auto it = x.nonzeros().cbegin();

    // New data vector
    std::vector<Scalar> d(nrow*ncol, val.scalar());

    // Copy nonzeros
    for (casadi_int cc=0; cc<ncol; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        d[cc*nrow + row[el]] = *it++;
      }
    }

    // Construct return matrix
    return Matrix<Scalar>(Sparsity::dense(x.size()), d);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::cumsum(const Matrix<Scalar> &x, casadi_int axis) {
    if (axis==-1) axis = x.is_row();
    Matrix<Scalar> ret = x;
    if (axis==0) {
      for (casadi_int i=1;i<x.size1();++i)
        ret(i, Slice()) += ret(i-1, Slice());
    } else {
      for (casadi_int i=1;i<x.size2();++i)
        ret(Slice(), i) += ret(Slice(), i-1);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::einstein(
      const Matrix<Scalar>& A, const Matrix<Scalar>& B, const Matrix<Scalar>& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    std::vector<casadi_int> iter_dims;
    std::vector<casadi_int> strides_a;
    std::vector<casadi_int> strides_b;
    std::vector<casadi_int> strides_c;
    casadi_int n_iter = einstein_process(A, B, C, dim_a, dim_b, dim_c, a, b, c,
          iter_dims, strides_a, strides_b, strides_c);

    const std::vector<Scalar>& Av = A.nonzeros();
    const std::vector<Scalar>& Bv = B.nonzeros();

    Matrix<Scalar> ret = C;
    std::vector<Scalar>& Cv = ret.nonzeros();

    einstein_eval(n_iter, iter_dims, strides_a, strides_b, strides_c,
      get_ptr(Av), get_ptr(Bv), get_ptr(Cv));
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::einstein(const Matrix<Scalar>& A, const Matrix<Scalar>& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    return Matrix<Scalar>::einstein(A, B, Matrix<Scalar>::zeros(product(dim_c), 1),
      dim_a, dim_b, dim_c, a, b, c);
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix() : sparsity_(Sparsity(0, 0)) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Matrix<Scalar>& m) : sparsity_(m.sparsity_), nonzeros_(m.nonzeros_) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::vector<Scalar>& x) :
      sparsity_(Sparsity::dense(x.size(), 1)), nonzeros_(x) {
  }

  template<typename Scalar>
  Matrix<Scalar>& Matrix<Scalar>::operator=(const Matrix<Scalar>& m) {
    sparsity_ = m.sparsity_;
    nonzeros_ = m.nonzeros_;
    return *this;
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::type_name() { return matrixName<Scalar>(); }

  template<typename Scalar>
  void Matrix<Scalar>::print_scalar(std::ostream &stream) const {
    casadi_assert(numel()==1, "Not a scalar");

    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();

    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }

    if (nnz()==0) {
      stream << "00";
    } else {
      stream << scalar();
    }
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_vector(std::ostream &stream, bool truncate) const {
    casadi_assert(is_column(), "Not a vector");

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Access data structures
    const casadi_int* row = this->row();
    casadi_int nnz = this->nnz();
    casadi_int size1 = this->size1();

    // No need to truncate if less than 1000 entries
    const casadi_int max_numel = 1000;
    if (truncate && size1<=max_numel) truncate=false;

    // Nonzero
    casadi_int el=0;

    // Loop over rows
    stream << "[";
    for (casadi_int rr=0; rr<size1; ++rr) {
      // String representation
      std::string s = el<nnz && rr==row[el] ? nz.at(el++) : "00";

      // Truncate?
      if (truncate && rr>=3 && rr<size1-3) {
        // Do not print
        if (rr==3) stream << ", ...";
      } else {
        // Print
        if (rr!=0) stream << ", ";
        stream << s;
      }
    }
    stream << "]" << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_dense(std::ostream &stream, bool truncate) const {
    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << "@" << (i+1) << "=" << inter[i] << ", ";
    inter.clear();

    // Access data structures
    casadi_int size1 = this->size1();
    casadi_int size2 = this->size2();
    const casadi_int* colind = this->colind();
    const casadi_int* row = this->row();

    // No need to truncate if less than 1000 entries
    const casadi_int max_numel = 1000;
    if (truncate && size1*size2<=max_numel) truncate=false;

    // Truncate rows and/or columns
    bool truncate_rows = truncate && size1>=7;
    bool truncate_columns = truncate && size2>=7;

    // Index counter for each column
    std::vector<casadi_int> ind(colind, colind+size2+1);

    // Print as a single line?
    bool oneliner=size1<=1;

    // Loop over rows
    for (casadi_int rr=0; rr<size1; ++rr) {
      // Print row?
      bool print_row = !(truncate_rows && rr>=3 && rr<size1-3);

      // Beginning of row
      if (rr==0) {
        if (!oneliner) stream << std::endl;
        stream << "[[";
      } else if (print_row) {
        stream << " [";
      }

      // Loop over columns
      for (casadi_int cc=0; cc<size2; ++cc) {
        // String representation of element
        std::string s = ind[cc]<colind[cc+1] && row[ind[cc]]==rr
          ? nz.at(ind[cc]++) : "00";

        // Skip whole row?
        if (!print_row) continue;

        // Print column?
        bool print_column = !(truncate_columns && cc>=3 && cc<size2-3);

        // Print element
        if (print_column) {
          if (cc!=0) stream << ", ";
          stream << s;
        } else if (cc==3) {
          stream << ", ...";
        }
      }

      // End of row
      if (rr<size1-1) {
        if (print_row) {
          stream << "], ";
          if (!oneliner) stream << std::endl;
        } else if (rr==3) {
          stream << " ...," << std::endl;
        }
      } else {
        stream << "]]";
      }
    }
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_sparse(std::ostream &stream, bool truncate) const {
    // Access data structures
    casadi_int size1 = this->size1();
    casadi_int size2 = this->size2();
    const casadi_int* colind = this->colind();
    const casadi_int* row = this->row();
    casadi_int nnz = this->nnz();

    // Quick return if all zero sparse
    if (nnz==0) {
      stream << "all zero sparse: " << size1 << "-by-" << size2 << std::flush;
      return;
    }

    // Print header
    stream << "sparse: " << size1 << "-by-" << size2 << ", " << nnz << " nnz";

    // Get components
    std::vector<std::string> nz, inter;
    print_split(nz, inter);

    // Print intermediate expressions
    for (casadi_int i=0; i<inter.size(); ++i)
      stream << std::endl << " @" << (i+1) << "=" << inter[i] << ",";
    inter.clear();

    // No need to truncate if less than 1000 nonzeros
    const casadi_int max_nnz = 1000;
    if (truncate && nnz<=max_nnz) truncate=false;

    // Print nonzeros
    for (casadi_int cc=0; cc<size2; ++cc) {
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        if (truncate && el>=3 && el<nnz-3) {
          if (el==3) stream << std::endl << " ...";
        } else {
          stream << std::endl << " (" << row[el] << ", " << cc << ") -> " << nz.at(el);
          InterruptHandler::check();
        }
      }
    }
    stream << std::flush;
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_split(std::vector<std::string>& nz,
                                    std::vector<std::string>& inter) const {
    nz.resize(nnz());
    inter.resize(0);

    // Temporary
    std::stringstream ss;
    ss.precision(stream_precision_);
    ss.width(stream_width_);
    if (stream_scientific_) {
      ss.setf(std::ios::scientific);
    } else {
      ss.unsetf(std::ios::scientific);
    }

    // Print nonzeros
    for (casadi_int i=0; i<nz.size(); ++i) {
      ss.str(std::string());
      ss << nonzeros().at(i);
      nz[i] = ss.str();
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::print_scalar(std::ostream &stream, const Scalar& e) {
    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();

    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }
    stream << e;
    stream << std::flush;

    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename Scalar>
  void Matrix<Scalar>::disp(std::ostream& stream, bool more) const {
    if (is_empty()) {
      stream << "[]";
    } else if (numel()==1) {
      print_scalar(stream);
    } else if (is_column()) {
      print_vector(stream);
    } else if (std::max(size1(), size2())<=10 ||
        static_cast<double>(nnz())/static_cast<double>(numel())>=0.5) {
      // if "small" or "dense"
      print_dense(stream);
    } else {
      print_sparse(stream);
    }
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::get_str(bool more) const {
    std::stringstream ss;
    disp(ss, more);
    return ss.str();
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(casadi_int nnz) {
    reserve(nnz, size2());
  }

  template<typename Scalar>
  void Matrix<Scalar>::reserve(casadi_int nnz, casadi_int ncol) {
    nonzeros().reserve(nnz);
  }

  template<typename Scalar>
  void Matrix<Scalar>::resize(casadi_int nrow, casadi_int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename Scalar>
  void Matrix<Scalar>::clear() {
    sparsity_ = Sparsity(0, 0);
    nonzeros().clear();
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(double val) :
      sparsity_(
        Sparsity::dense(1, 1)),
        nonzeros_(std::vector<Scalar>(1, static_cast<Scalar>(val))) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::vector< std::vector<double> >& d) {
    // Get dimensions
    casadi_int nrow=d.size();
    casadi_int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (casadi_int rr=0; rr<nrow; ++rr) {
      casadi_assert(ncol==d[rr].size(),
        "Shape mismatch.\n"
        "Attempting to construct a matrix from a nested list.\n"
        "I got convinced that the desired size is (" + str(nrow) + " x " + str(ncol)
        + " ), but now I encounter a vector of size (" + str(d[rr].size()) +  " )");
    }

    // Form matrix
    sparsity_ = Sparsity::dense(nrow, ncol);
    nonzeros().resize(nrow*ncol);
    typename std::vector<Scalar>::iterator it=nonzeros_.begin();
    for (casadi_int cc=0; cc<ncol; ++cc) {
      for (casadi_int rr=0; rr<nrow; ++rr) {
        *it++ = static_cast<Scalar>(d[rr][cc]);
      }
    }
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp) : sparsity_(sp), nonzeros_(sp.nnz(), 1) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(casadi_int nrow, casadi_int ncol) : sparsity_(nrow, ncol) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const std::pair<casadi_int, casadi_int>& rc) : sparsity_(rc) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Scalar& val, bool dummy) :
      sparsity_(sp), nonzeros_(sp.nnz(), val) {
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const std::vector<Scalar>& d, bool dummy) :
      sparsity_(sp), nonzeros_(d) {
    casadi_assert(sp.nnz()==d.size(), "Size mismatch.\n"
                          "You supplied a sparsity of " + sp.dim()
                          + ", but the supplied vector is of length " + str(d.size()));
  }

  template<typename Scalar>
  Matrix<Scalar>::Matrix(const Sparsity& sp, const Matrix<Scalar>& d) {
    if (d.is_scalar()) {
      *this = Matrix<Scalar>(sp, d.scalar(), false);
    } else if (d.is_column() || d.size1()==1) {
      casadi_assert_dev(sp.nnz()==d.numel());
      if (d.is_dense()) {
        *this = Matrix<Scalar>(sp, d.nonzeros(), false);
      } else {
        *this = Matrix<Scalar>(sp, densify(d).nonzeros(), false);
      }
    } else {
      casadi_error("Matrix(Sparsity, Matrix): Only allowed for scalars and vectors");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::unary(casadi_int op, const Matrix<Scalar> &x) {
    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<x.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_data[el], x_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.is_dense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, 0, 0, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::operator-() const {
    return unary(OP_NEG, *this);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::operator+() const {
    return *this;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mrdivide(const Matrix<Scalar>& b,
                                              const Matrix<Scalar>& a) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a.T(), b.T()).T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mldivide(const Matrix<Scalar>& a,
                                              const Matrix<Scalar>& b) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a, b);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::printme(const Matrix<Scalar>& y) const {
    return binary(OP_PRINTME, *this, y);
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc, bool ind1) {
    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sparsity_.erase(rr, cc, ind1);

    // Update non-zero entries
    for (casadi_int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::erase(const std::vector<casadi_int>& rr, bool ind1) {
    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sparsity_.erase(rr, ind1);

    // Update non-zero entries
    for (casadi_int k=0; k<mapping.size(); ++k)
      nonzeros()[k] = nonzeros()[mapping[k]];

    // Truncate nonzero vector
    nonzeros().resize(mapping.size());
  }

  template<typename Scalar>
  void Matrix<Scalar>::remove(const std::vector<casadi_int>& rr,
      const std::vector<casadi_int>& cc) {
    casadi_assert_bounded(rr, size1());
    casadi_assert_bounded(cc, size2());

    // Remove by performing a complementary slice
    std::vector<casadi_int> rrc = complement(rr, size1());
    std::vector<casadi_int> ccc = complement(cc, size2());

    Matrix<Scalar> ret = (*this)(rrc, ccc);

    operator=(ret);

  }

  template<typename Scalar>
  void Matrix<Scalar>::enlarge(casadi_int nrow, casadi_int ncol, const std::vector<casadi_int>& rr,
                                 const std::vector<casadi_int>& cc, bool ind1) {
    sparsity_.enlarge(nrow, ncol, rr, cc, ind1);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtimes(const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return x*y;
    } else {
      Matrix<Scalar> z = Matrix<Scalar>::zeros(Sparsity::mtimes(x.sparsity(), y.sparsity()));
      return mac(x, y, z);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mac(const Matrix<Scalar> &x,
                                         const Matrix<Scalar> &y,
                                         const Matrix<Scalar> &z) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + x*y;
    }

    // Check matching dimensions
    casadi_assert(x.size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          + x.dim() + " and rhs is " + y.dim() + ".");

    casadi_assert(y.size2()==z.size2(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          + mtimes(x, y).dim() + " and rhs is " + z.dim() + ".");

    casadi_assert(x.size1()==z.size1(),
                          "Matrix addition with incompatible dimensions. Lhs is "
                          + mtimes(x, y).dim() + " and rhs is " + z.dim() + ".");

    // Check if we can simplify the product
    if (x.is_eye()) {
      return y + z;
    } else if (y.is_eye()) {
      return x + z;
    } else if (x.is_zero() || y.is_zero()) {
      return z;
    } else {
      // Carry out the matrix product
      Matrix<Scalar> ret = z;
      std::vector<Scalar> work(x.size1());
      casadi_mtimes(x.ptr(), x.sparsity(), y.ptr(), y.sparsity(),
                    ret.ptr(), ret.sparsity(), get_ptr(work), false);
      return ret;
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  _bilin(const Matrix<Scalar>& A, const Matrix<Scalar>& x,
         const Matrix<Scalar>& y) {
    return casadi_bilin(A.ptr(), A.sparsity(), x.ptr(), y.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  _rank1(const Matrix<Scalar>& A, const Matrix<Scalar>& alpha,
         const Matrix<Scalar>& x, const Matrix<Scalar>& y) {
    Matrix<Scalar> ret = A;
    casadi_rank1(ret.ptr(), ret.sparsity(), *alpha.ptr(), x.ptr(), y.ptr());
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::T() const {
    // quick return if empty or scalar
    if ((size1()==0 && size2()==0) || is_scalar()) return *this;

    // Create the new sparsity pattern and the mapping
    std::vector<casadi_int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<Scalar> ret = zeros(s);

    // Copy the content
    for (casadi_int i=0; i<mapping.size(); ++i)
      ret->at(i) = nonzeros().at(mapping[i]);

    return ret;
  }

  template<typename Scalar>
  const Scalar Matrix<Scalar>::scalar() const {
    // Make sure that the matrix is 1-by-1
    casadi_assert(is_scalar(), "Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if (nnz()==1)
      return nonzeros()[0];
    else
      return casadi_limits<Scalar>::zero;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::binary(casadi_int op,
                                            const Matrix<Scalar> &x,
                                            const Matrix<Scalar> &y) {
    if (x.is_scalar()) {
      return scalar_matrix(op, x, y);
    } else if (y.is_scalar()) {
      return matrix_scalar(op, x, y);
    } else {
      return matrix_matrix(op, x, y);
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  scalar_matrix(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    if ( (operation_checker<FX0Checker>(op) && y.nnz()==0) ||
         (operation_checker<F0XChecker>(op) && x.nnz()==0))
            return Matrix<Scalar>::zeros(Sparsity(y.size()));

    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(y.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const Scalar& x_val = x_data.empty() ? casadi_limits<Scalar>::zero : x->front();
    const std::vector<Scalar>& y_data = y.nonzeros();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<y.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_val, y_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!y.is_dense() && !operation_checker<FX0Checker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, x_val, casadi_limits<Scalar>::zero, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  matrix_scalar(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {

    if ( (operation_checker<FX0Checker>(op) && y.nnz()==0) ||
         (operation_checker<F0XChecker>(op) && x.nnz()==0))
            return Matrix<Scalar>::zeros(Sparsity(x.size()));

    // Return value
    Matrix<Scalar> ret = Matrix<Scalar>::zeros(x.sparsity());

    // Nonzeros
    std::vector<Scalar>& ret_data = ret.nonzeros();
    const std::vector<Scalar>& x_data = x.nonzeros();
    const std::vector<Scalar>& y_data = y.nonzeros();
    const Scalar& y_val = y_data.empty() ? casadi_limits<Scalar>::zero : y->front();

    // Do the operation on all non-zero elements
    for (casadi_int el=0; el<x.nnz(); ++el) {
      casadi_math<Scalar>::fun(op, x_data[el], y_val, ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.is_dense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, casadi_limits<Scalar>::zero, y_val, fcn_0);
      if (!casadi_limits<Scalar>::is_zero(fcn_0)) { // Remove this if?
        ret = densify(ret, fcn_0);
      }
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  matrix_matrix(casadi_int op, const Matrix<Scalar> &x, const Matrix<Scalar> &y) {
    // Check, correct dimensions
    if (x.size() != y.size()) {
      // x and y are horizontal multiples of each other?
      if (!x.is_empty() && !y.is_empty()) {
        if (x.size1() == y.size1() && x.size2() % y.size2() == 0) {
          return matrix_matrix(op, x, repmat(y, 1, x.size2() / y.size2()));
        } else if (y.size1() == x.size1() && y.size2() % x.size2() == 0) {
          return matrix_matrix(op, repmat(x, 1, y.size2() / x.size2()), y);
        }
      }
      // Dimension mismatch
      casadi_error("Dimension mismatch for " + casadi_math<Scalar>::print(op, "x", "y") +
                   ", x is " + x.dim() + ", while y is " + y.dim());
    }

    // Get the sparsity pattern of the result
    // (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.combine(y_sp, operation_checker<F0XChecker>(op),
                                 operation_checker<FX0Checker>(op));

    // Return value
    Matrix<Scalar> r = zeros(r_sp);

    // Perform the operations elementwise
    if (x_sp==y_sp) {
      // Matching sparsities
      casadi_math<Scalar>::fun(op, x.ptr(), y.ptr(), r.ptr(), r_sp.nnz());
    } else if (y_sp==r_sp) {
      // Project first argument
      Matrix<Scalar> x_mod = x(r_sp);
      casadi_math<Scalar>::fun(op, x_mod.ptr(), y.ptr(), r.ptr(), r_sp.nnz());
    } else if (x_sp==r_sp) {
      // Project second argument
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, x.ptr(), y_mod.ptr(), r.ptr(), r_sp.nnz());
    } else {
      // Project both arguments
      Matrix<Scalar> x_mod = x(r_sp);
      Matrix<Scalar> y_mod = y(r_sp);
      casadi_math<Scalar>::fun(op, x_mod.ptr(), y_mod.ptr(), r.ptr(), r_sp.nnz());
    }

    // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
    if (!r.is_dense() && !operation_checker<F00Checker>(op)) {
      // Get the value for the structural zeros
      Scalar fcn_0;
      casadi_math<Scalar>::fun(op, casadi_limits<Scalar>::zero,
                               casadi_limits<Scalar>::zero, fcn_0);
      r = densify(r, fcn_0);
    }

    return r;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d) {
    return triplet(row, col, d, *std::max_element(row.begin(), row.end()),
                   *std::max_element(col.begin(), col.end()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d,
                                             const std::pair<casadi_int, casadi_int>& rc) {
    return triplet(row, col, d, rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triplet(const std::vector<casadi_int>& row,
                                             const std::vector<casadi_int>& col,
                                             const Matrix<Scalar>& d,
                                             casadi_int nrow, casadi_int ncol) {
    casadi_assert(col.size()==row.size() && col.size()==d.nnz(),
                          "Argument error in Matrix<Scalar>::triplet(row, col, d): "
                          "supplied lists must all be of equal length, but got: "
                          + str(row.size()) + ", " + str(col.size()) + " and " + str(d.nnz()));
    std::vector<casadi_int> mapping;
    Sparsity sp = Sparsity::triplet(nrow, ncol, row, col, mapping, false);
    return Matrix<Scalar>(sp, d.nz(mapping));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eye(casadi_int n) {
    return Matrix<Scalar>::ones(Sparsity::diag(n));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const Sparsity& sp) {
    casadi_assert(std::numeric_limits<Scalar>::has_infinity,
                          "Datatype cannot represent infinity");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::infinity(), false);
  }


  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(const std::pair<casadi_int, casadi_int>& rc) {
    return inf(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inf(casadi_int nrow, casadi_int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const Sparsity& sp) {
    casadi_assert(std::numeric_limits<Scalar>::has_quiet_NaN,
                          "Datatype cannot represent not-a-number");
    return Matrix<Scalar>(sp, std::numeric_limits<Scalar>::quiet_NaN(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(const std::pair<casadi_int, casadi_int>& rc) {
    return nan(rc.first, rc.second);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nan(casadi_int nrow, casadi_int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_regular() const {
    return casadi::is_regular(nonzeros_);
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_smooth() const {
    return true;
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::element_hash() const {
    casadi_error("'element_hash' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_leaf() const {
    casadi_error("'is_leaf' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_commutative() const {
    casadi_error("'is_commutative' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_symbolic() const {
    return false;
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::op() const {
    casadi_error("'op' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_op(casadi_int k) const {
    casadi_error("'is_op' not defined for " + type_name());
  }

  template<typename Scalar>
  void Matrix<Scalar>::export_code(const std::string& lang,
       std::ostream &stream, const Dict& options) const {
    casadi_error("'export_code' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_valid_input() const {
    return false;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::has_duplicates() const {
    casadi_error("'has_duplicates' not defined for " + type_name());
  }

  template<typename Scalar>
  void Matrix<Scalar>::reset_input() const {
    casadi_error("'reset_input' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_integer() const {
    // Look for non-integers
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_integer(e)) return false;

    // Integer if reached this point
    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_constant() const {
    // Look for non-constants
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_constant(e)) return false;

    // Constant if we reach this point
    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_one() const {
    if (!is_dense()) return false;

    // Look for non-ones
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_minus_one() const {
    if (!is_dense()) return false;

    // Look for non-minus-ones
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_minus_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_zero() const {

    // Look for non-zeros
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_zero(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_eye() const {

    // Make sure that the matrix is diagonal
    if (!sparsity().is_diag()) return false;

    // Make sure that all entries are one
    for (auto&& e : nonzeros()) if (!casadi_limits<Scalar>::is_one(e)) return false;

    return true;
  }

  template<typename Scalar>
  bool Matrix<Scalar>::is_equal(const Matrix<Scalar> &x, const Matrix<Scalar> &y,
      casadi_int depth) {
    // Assert matching dimensions
    casadi_assert(x.size() == y.size(), "Dimension mismatch");

    // Project to union of patterns and call recursively if different sparsity
    if (x.sparsity() != y.sparsity()) {
      Sparsity sp = x.sparsity() + y.sparsity();
      return is_equal(project(x, sp), project(y, sp), depth);
    }

    // Check individual elements
    auto y_it = y.nonzeros().begin();
    for (auto&& e : x.nonzeros()) {
      if (!casadi_limits<Scalar>::is_equal(e, *y_it++, depth)) return false;
    }

    // True if reched this point
    return true;
  }

  // To avoid overloaded function name conflicts
  template<typename Scalar>
  inline Matrix<Scalar> mmin_nonstatic(const Matrix<Scalar> &x) {
    if (x.is_empty()) return Matrix<Scalar>();
    return casadi_mmin(x.ptr(), x.nnz(), x.is_dense());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mmin(const Matrix<Scalar> &x) {
    return mmin_nonstatic(x);
  }

  // To avoid overloaded function name conflicts
  template<typename Scalar>
  inline Matrix<Scalar> mmax_nonstatic(const Matrix<Scalar> &x) {
    if (x.is_empty()) return Matrix<Scalar>();
    return casadi_mmax(x.ptr(), x.nnz(), x.is_dense());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mmax(const Matrix<Scalar> &x) {
    return mmax_nonstatic(x);
  }

  template<typename Scalar>
  bool Matrix<Scalar>::has_zeros() const {
    // Check if the structural nonzero is known to be zero
    for (auto&& e : nonzeros()) if (casadi_limits<Scalar>::is_zero(e)) return true;

    // No known zeros amongst the structurally nonzero entries
    return false;
  }

  template<typename Scalar>
  std::string Matrix<Scalar>::name() const {
    casadi_error("'name' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dep(casadi_int ch) const {
    casadi_error("'dep' not defined for " + type_name());
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::n_dep() const {
    casadi_error("'n_dep' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::project(const Matrix<Scalar>& x,
                                         const Sparsity& sp, bool intersect) {
    if (intersect) {
      return project(x, sp.intersect(x.sparsity()), false);
    } else {
      casadi_assert(sp.size()==x.size(), "Dimension mismatch");
      Matrix<Scalar> ret = Matrix<Scalar>::zeros(sp);
      std::vector<Scalar> w(x.size1());
      casadi_project(x.ptr(), x.sparsity(), ret.ptr(), sp, get_ptr(w));
      return ret;
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::set_max_depth(casadi_int eq_depth) {
    casadi_error("'set_max_depth' not defined for " + type_name());
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::get_max_depth() {
    casadi_error("'get_max_depth' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::det(const Matrix<Scalar>& x) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "matrix must be square");

    // Trivial return if scalar
    if (x.is_scalar()) return x;

    // Trivial case 2 x 2
    if (n==2) return x(0, 0) * x(1, 1) - x(0, 1) * x(1, 0);

    // Return expression
    Matrix<Scalar> ret = 0;

    // Find out which is the best direction to expand along

    // Build up an IM with ones on the non-zeros
    Matrix<casadi_int> sp = IM::ones(x.sparsity());

    // Have a count of the nonzeros for each row
    Matrix<casadi_int> row_count = Matrix<casadi_int>::sum2(sp);

    // A blank row? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    // Have a count of the nonzeros for each col
    Matrix<casadi_int> col_count = Matrix<casadi_int>::sum1(sp).T();

    // A blank col? determinant is structurally zero
    if (!row_count.is_dense()) return 0;

    casadi_int min_row = std::distance(row_count.nonzeros().begin(),
                                std::min_element(row_count.nonzeros().begin(),
                                                 row_count.nonzeros().end()));
    casadi_int min_col = std::distance(col_count.nonzeros().begin(),
                                std::min_element(col_count.nonzeros().begin(),
                                                 col_count.nonzeros().end()));

    if (min_row <= min_col) {
      // Expand along row j
      casadi_int j = row_count.sparsity().row(min_row);

      Matrix<Scalar> row = x(j, Slice(0, n));

      std::vector< casadi_int > col_i = row.sparsity().get_col();

      for (casadi_int k=0; k<row.nnz(); ++k) {
        // Sum up the cofactors
        ret += row->at(k)*cofactor(x, col_i.at(k), j);
      }
      return ret;
    } else {
      // Expand along col i
      casadi_int i = col_count.sparsity().row(min_col);

      Matrix<Scalar> col = x(Slice(0, n), i);

      const casadi_int* row_i = col.row();

      for (casadi_int k=0; k<col.nnz(); ++k) {
        // Sum up the cofactors
        ret += col->at(k)*cofactor(x, i, row_i[k]);
      }
      return ret;
    }

  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sum2(const Matrix<Scalar>& x) {
    return mtimes(x, Matrix<Scalar>::ones(x.size2(), 1));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sum1(const Matrix<Scalar>& x) {
    return mtimes(Matrix<Scalar>::ones(1, x.size1()), x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::minor(const Matrix<Scalar>& x,
                                              casadi_int i, casadi_int j) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "minor: matrix must be square");

    // Trivial return if scalar
    if (n==1) return 1;

    // Remove col i and row j
    Matrix<Scalar> M = Matrix<Scalar>(n-1, n-1);

    std::vector<casadi_int> col = x.sparsity().get_col();
    const casadi_int* row = x.sparsity().row();

    for (casadi_int k=0; k<x.nnz(); ++k) {
      casadi_int i1 = col[k];
      casadi_int j1 = row[k];

      if (i1 == i || j1 == j) continue;

      casadi_int i2 = (i1<i)?i1:i1-1;
      casadi_int j2 = (j1<j)?j1:j1-1;

      M(j2, i2) = x(j1, i1);
    }
    return det(M);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::cofactor(const Matrix<Scalar>& A, casadi_int i, casadi_int j) {

    // Calculate the i, j minor
    Matrix<Scalar> minor_ij = minor(A, i, j);
    // Calculate the cofactor
    casadi_int sign_i = 1-2*((i+j) % 2);

    return sign_i * minor_ij;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::adj(const Matrix<Scalar>& x) {
    casadi_int n = x.size2();
    casadi_assert(n == x.size1(), "adj: matrix must be square");

    // Temporary placeholder
    Matrix<Scalar> temp;

    // Cofactor matrix
    Matrix<Scalar> C = Matrix<Scalar>(n, n);
    for (casadi_int i=0; i<n; ++i)
      for (casadi_int j=0; j<n; ++j) {
        temp = cofactor(x, i, j);
        if (!temp.is_zero()) C(j, i) = temp;
      }

    return C.T();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::inv_minor(const Matrix<Scalar>& x) {
    // laplace formula
    return adj(x)/det(x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x,
      casadi_int nrow, casadi_int ncol) {
    Sparsity sp = Sparsity::reshape(x.sparsity(), nrow, ncol);
    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::reshape(const Matrix<Scalar>& x, const Sparsity& sp) {
    // quick return if already the right shape
    if (sp==x.sparsity()) return x;

    // make sure that the patterns match
    casadi_assert_dev(sp.is_reshape(x.sparsity()));

    return Matrix<Scalar>(sp, x.nonzeros(), false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::trace(const Matrix<Scalar>& x) {
    casadi_assert(x.is_square(), "trace: must be square");
    Scalar res=0;
    const Scalar* d=x.ptr();
    casadi_int size2 = x.size2();
    const casadi_int *colind=x.colind(), *row=x.row();
    for (casadi_int c=0; c<size2; c++) {
      for (casadi_int k=colind[c]; k!=colind[c+1]; ++k) {
        if (row[k]==c) {
          res += d[k];
        }
      }
    }
    return res;
  }

  template<typename Scalar>
  Matrix<Scalar>
  Matrix<Scalar>::blockcat(const std::vector< std::vector<Matrix<Scalar> > > &v) {
    std::vector< Matrix<Scalar> > ret;
    for (casadi_int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::horzcat(const std::vector<Matrix<Scalar> > &v) {
    // Concatenate sparsity patterns
    std::vector<Sparsity> sp(v.size());
    for (casadi_int i=0; i<v.size(); ++i) sp[i] = v[i].sparsity();
    Matrix<Scalar> ret = zeros(Sparsity::horzcat(sp));

    // Copy nonzeros
    auto i=ret->begin();
    for (auto&& j : v) {
      std::copy(j->begin(), j->end(), i);
      i += j.nnz();
    }
    return ret;
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> >
  Matrix<Scalar>::horzsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset) {
    // Split up the sparsity pattern
    std::vector<Sparsity> sp = Sparsity::horzsplit(x.sparsity(), offset);

    // Return object
    std::vector<Matrix<Scalar> > ret;
    ret.reserve(sp.size());

    // Copy data
    auto i=x.nonzeros().begin();
    for (auto&& j : sp) {
      auto i_next = i + j.nnz();
      ret.push_back(Matrix<Scalar>(j, std::vector<Scalar>(i, i_next), false));
      i = i_next;
    }

    // Return the assembled matrix
    casadi_assert_dev(i==x.nonzeros().end());
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::vertcat(const std::vector<Matrix<Scalar> > &v) {
    std::vector<Matrix<Scalar> > vT(v.size());
    for (casadi_int i=0; i<v.size(); ++i) vT[i] = v[i].T();
    return horzcat(vT).T();
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::vertsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset) {
    std::vector< Matrix<Scalar> > ret = horzsplit(x.T(), offset);
    for (auto&& e : ret) e = e.T();
    return ret;
  }

  template<typename Scalar>
  std::vector< Matrix<Scalar> >
  Matrix<Scalar>::diagsplit(const Matrix<Scalar>& x, const std::vector<casadi_int>& offset1,
                              const std::vector<casadi_int>& offset2) {
    // Consistency check
    casadi_assert_dev(offset1.size()>=1);
    casadi_assert_dev(offset1.front()==0);
    casadi_assert_dev(offset1.back()==x.size1());
    casadi_assert_dev(is_monotone(offset1));

    // Consistency check
    casadi_assert_dev(offset2.size()>=1);
    casadi_assert_dev(offset2.front()==0);
    casadi_assert_dev(offset2.back()==x.size2());
    casadi_assert_dev(is_monotone(offset2));

    // Number of outputs
    casadi_int n = offset1.size()-1;

    // Return value
    std::vector< Matrix<Scalar> > ret;

    // Caveat: this is a very silly implementation
    for (casadi_int i=0; i<n; ++i) {
      ret.push_back(x(Slice(offset1[i], offset1[i+1]), Slice(offset2[i], offset2[i+1])));
    }

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::dot(const Matrix<Scalar> &x,
                                     const Matrix<Scalar> &y) {
    casadi_assert(x.size()==y.size(), "dot: Dimension mismatch");
    if (x.sparsity()!=y.sparsity()) {
      Sparsity sp = x.sparsity() * y.sparsity();
      return dot(project(x, sp), project(y, sp));
    }
    return casadi_dot(x.nnz(), x.ptr(), y.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::all(const Matrix<Scalar>& x) {
    if (!x.is_dense()) return false;
    Scalar ret=1;
    for (casadi_int i=0; i<x.nnz(); ++i) {
      ret = ret && x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::any(const Matrix<Scalar>& x) {
    if (!x.is_dense()) return false;
    Scalar ret=0;
    for (casadi_int i=0; i<x.nnz(); ++i) {
      ret = ret || x->at(i)==1;
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_1(const Matrix<Scalar>& x) {
    return casadi_norm_1(x.nnz(), x.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_2(const Matrix<Scalar>& x) {
    if (x.is_vector()) {
      return norm_fro(x);
    } else {
      casadi_error("2-norms currently only supported for vectors. "
                   "Did you intend to calculate a Frobenius norms (norm_fro)?");
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_fro(const Matrix<Scalar>& x) {
    return casadi_norm_2(x.nnz(), x.ptr());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_inf(const Matrix<Scalar>& x) {
    // Get largest element by absolute value
    Matrix<Scalar> s = 0;
    for (auto i=x.nonzeros().begin(); i!=x.nonzeros().end(); ++i) {
      s = fmax(s, fabs(Matrix<Scalar>(*i)));
    }
    return s;
  }

  template<typename Scalar>
  void Matrix<Scalar>::
  qr_sparse(const Matrix<Scalar>& A,
    Matrix<Scalar>& V, Matrix<Scalar> &R, Matrix<Scalar>& beta,
    std::vector<casadi_int>& prinv, std::vector<casadi_int>& pc, bool amd) {
    // Calculate the pattern
    Sparsity spV, spR;
    A.sparsity().qr_sparse(spV, spR, prinv, pc, amd);
    // Calculate the nonzeros
    casadi_int nrow_ext = spV.size1(), ncol = spV.size2();
    V = nan(spV);
    R = nan(spR);
    beta = nan(ncol, 1);
    vector<Scalar> w(nrow_ext);
    casadi_qr(A.sparsity(), A.ptr(), get_ptr(w), spV, V.ptr(),
              spR, R.ptr(), beta.ptr(),
              get_ptr(prinv), get_ptr(pc));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  qr_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& v,
           const Matrix<Scalar>& r, const Matrix<Scalar>& beta,
           const std::vector<casadi_int>& prinv, const std::vector<casadi_int>& pc,
           bool tr) {
    // Get dimensions, check consistency
    casadi_int ncol = v.size2();
    casadi_int nrow = b.size1(), nrhs = b.size2();
    casadi_assert(r.size()==v.size(), "'r', 'v' dimension mismatch");
    casadi_assert(beta.is_vector() && beta.numel()==ncol, "'beta' has wrong dimension");
    casadi_assert(prinv.size()==r.size1(), "'pinv' has wrong dimension");
    // Work vector
    std::vector<Scalar> w(nrow+ncol);
    // Return value
    Matrix<Scalar> x = densify(b);
    casadi_qr_solve(x.ptr(), nrhs, tr, v.sparsity(), v.ptr(), r.sparsity(), r.ptr(),
                    beta.ptr(), get_ptr(prinv), get_ptr(pc), get_ptr(w));
    return x;
  }

  template<typename Scalar>
  void Matrix<Scalar>::qr(const Matrix<Scalar>& A,
                            Matrix<Scalar>& Q, Matrix<Scalar> &R) {
    // The following algorithm is taken from J. Demmel:
    // Applied Numerical Linear Algebra (algorithm 3.1.)
    casadi_assert(A.size1()>=A.size2(), "qr: fewer rows than columns");

    // compute Q and R column by column
    Q = R = Matrix<Scalar>();
    for (casadi_int i=0; i<A.size2(); ++i) {
      // Initialize qi to be the i-th column of *this
      Matrix<Scalar> ai = A(Slice(), i);
      Matrix<Scalar> qi = ai;
      // The i-th column of R
      Matrix<Scalar> ri = Matrix<Scalar>(A.size2(), 1);

      // subtract the projection of qi in the previous directions from ai
      for (casadi_int j=0; j<i; ++j) {

        // Get the j-th column of Q
        Matrix<Scalar> qj = Q(Slice(), j);

        ri(j, 0) = mtimes(qi.T(), qj); // Modified Gram-Schmidt
        // ri[j] = dot(qj, ai); // Classical Gram-Schmidt

        // Remove projection in direction j
        if (ri.has_nz(j, 0))
          qi -= ri(j, 0) * qj;
      }

      // Normalize qi
      ri(i, 0) = norm_2(qi);
      qi /= ri(i, 0);

      // Update R and Q
      Q = Matrix<Scalar>::horzcat({Q, qi});
      R = Matrix<Scalar>::horzcat({R, ri});
    }
  }

  template<typename Scalar>
  void Matrix<Scalar>::ldl(const Matrix<Scalar>& A, Matrix<Scalar> &D,
    Matrix<Scalar>& LT, std::vector<casadi_int>& p, bool amd) {
    // Symbolic factorization
    Sparsity Lt_sp = A.sparsity().ldl(p, amd);

    // Get dimension
    casadi_int n=A.size1();

    // Calculate entries in L and D
    vector<Scalar> D_nz(n), L_nz(Lt_sp.nnz()), w(n);
    casadi_ldl(A.sparsity(), get_ptr(A.nonzeros()), Lt_sp,
              get_ptr(L_nz), get_ptr(D_nz), get_ptr(p), get_ptr(w));

    // Assemble L and D
    LT = Matrix<Scalar>(Lt_sp, L_nz);
    D = D_nz;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  ldl_solve(const Matrix<Scalar>& b, const Matrix<Scalar>& D, const Matrix<Scalar>& LT,
            const std::vector<casadi_int>& p) {
    // Get dimensions, check consistency
    casadi_int n = b.size1(), nrhs = b.size2();
    casadi_assert(p.size()==n, "'p' has wrong dimension");
    casadi_assert(LT.size1()==n && LT.size2()==n, "'LT' has wrong dimension");
    casadi_assert(D.is_vector() && D.numel()==n, "'D' has wrong dimension");
    // Solve for all right-hand-sides
    Matrix<Scalar> x = densify(b);
    std::vector<Scalar> w(n);
    casadi_ldl_solve(x.ptr(), nrhs, LT.sparsity(), LT.ptr(), D.ptr(), get_ptr(p), get_ptr(w));
    return x;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::nullspace(const Matrix<Scalar>& A) {
    Matrix<Scalar> X = A;
    casadi_int n = X.size1();
    casadi_int m = X.size2();
    casadi_assert(m>=n, "nullspace(): expecting a flat matrix (more columns than rows), "
                          "but got " + str(X.dim()) + ".");

    Matrix<Scalar> seed = DM::eye(m)(Slice(0, m), Slice(n, m));

    std::vector< Matrix<Scalar> > us;
    std::vector< Matrix<Scalar> > betas;

    Matrix<Scalar> beta;

    for (casadi_int i=0;i<n;++i) {
      Matrix<Scalar> x = X(i, Slice(i, m));
      Matrix<Scalar> u = Matrix<Scalar>(x);
      Matrix<Scalar> sigma = sqrt(sum2(x*x));
      const Matrix<Scalar>& x0 = x(0, 0);
      u(0, 0) = 1;

      Matrix<Scalar> b = -copysign(sigma, x0);

      u(Slice(0), Slice(1, m-i))*= 1/(x0-b);
      beta = 1-x0/b;

      X(Slice(i, n), Slice(i, m)) -=
        beta*mtimes(mtimes(X(Slice(i, n), Slice(i, m)), u.T()), u);
      us.push_back(u);
      betas.push_back(beta);
    }

    for (casadi_int i=n-1;i>=0;--i) {
      seed(Slice(i, m), Slice(0, m-n)) -=
        betas[i]*mtimes(us[i].T(), mtimes(us[i], seed(Slice(i, m), Slice(0, m-n))));
    }

    return seed;

  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::chol(const Matrix<Scalar>& A) {
    // Perform an LDL transformation
    Matrix<Scalar> D, LT;
    std::vector<casadi_int> p;
    ldl(A, D, LT, p, false);
    // Add unit diagonal
    LT += Matrix<Scalar>::eye(D.size1());
    // Get the cholesky factor: R*R' = L*D*L' = (sqrt(D)*L')'*(sqrt(D)*L')
    return mtimes(diag(sqrt(D)), LT);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    // check dimensions
    casadi_assert(a.size1() == b.size1(), "solve Ax=b: dimension mismatch: b has "
                          + str(b.size1()) + " rows while A has " + str(a.size1()) + ".");
    casadi_assert(a.size1() == a.size2(), "solve: A not square but " + str(a.dim()));

    if (a.is_tril()) {
      // forward substitution if lower triangular
      Matrix<Scalar> x = b;
      const casadi_int*  Arow = a.row();
      const casadi_int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (casadi_int i=0; i<a.size2(); ++i) { // loop over columns forwards
        for (casadi_int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (casadi_int kk=Acolind[i+1]-1; kk>=Acolind[i] && Arow[kk]>i; --kk) {
            casadi_int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (a.is_triu()) {
      // backward substitution if upper triangular
      Matrix<Scalar> x = b;
      const casadi_int*  Arow = a.row();
      const casadi_int*  Acolind = a.colind();
      const std::vector<Scalar> & Adata = a.nonzeros();
      for (casadi_int i=a.size2()-1; i>=0; --i) { // loop over columns backwards
        for (casadi_int k=0; k<b.size2(); ++k) { // for every right hand side
          if (!x.has_nz(i, k)) continue;
          x(i, k) /= a(i, i);
          for (casadi_int kk=Acolind[i]; kk<Acolind[i+1] && Arow[kk]<i; ++kk) {
            casadi_int j = Arow[kk];
            x(j, k) -= Adata[kk]*x(i, k);
          }
        }
      }
      return x;
    } else if (a.has_zeros()) {

      // If there are structurally nonzero entries that are known to be zero,
      // remove these and rerun the algorithm
      return solve(sparsify(a), b);

    } else {

      // Make a BLT transformation of A
      std::vector<casadi_int> rowperm, colperm, rowblock, colblock;
      std::vector<casadi_int> coarse_rowblock, coarse_colblock;
      a.sparsity().btf(rowperm, colperm, rowblock, colblock,
                       coarse_rowblock, coarse_colblock);

      // Permute the right hand side
      Matrix<Scalar> bperm = b(rowperm, Slice());

      // Permute the linear system
      Matrix<Scalar> Aperm = a(rowperm, colperm);

      // Solution
      Matrix<Scalar> xperm;

      // Solve permuted system
      if (Aperm.is_tril()) {

        // Forward substitution if lower triangular
        xperm = solve(Aperm, bperm);

      } else if (a.size2()<=3) {

        // Form inverse by minor expansion and multiply if very small (up to 3-by-3)
        xperm = mtimes(inv_minor(Aperm), bperm);

      } else {

        // Make a QR factorization
        Matrix<Scalar> Q, R;
        qr(Aperm, Q, R);

        // Solve the factorized system (note that solve will now be fast since it is triangular)
        xperm = solve(R, mtimes(Q.T(), bperm));
      }

      // get the inverted column permutation
      std::vector<casadi_int> inv_colperm(colperm.size());
      for (casadi_int k=0; k<colperm.size(); ++k)
        inv_colperm[colperm[k]] = k;

      // Permute back the solution and return
      Matrix<Scalar> x = xperm(inv_colperm, Slice());
      return x;
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  solve(const Matrix<Scalar>& a, const Matrix<Scalar>& b,
           const std::string& lsolver, const Dict& dict) {
    casadi_error("'solve' with plugin not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  inv(const Matrix<Scalar>& a) {
    return solve(a, Matrix<Scalar>::eye(a.size1()));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  inv(const Matrix<Scalar>& a,
           const std::string& lsolver, const Dict& dict) {
    casadi_error("'inv' with plugin not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pinv(const Matrix<Scalar>& A) {
    if (A.size2()>=A.size1()) {
      return solve(mtimes(A, A.T()), A).T();
    } else {
      return solve(mtimes(A.T(), A), A.T());
    }
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  pinv(const Matrix<Scalar>& A, const std::string& lsolver, const Dict& dict) {
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  expm_const(const Matrix<Scalar>& A, const Matrix<Scalar>& t) {
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  expm(const Matrix<Scalar>& A) {
    casadi_error("'solve' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::kron(const Matrix<Scalar>& a, const Matrix<Scalar>& b) {
    Sparsity sp_ret = Sparsity::kron(a.sparsity(), b.sparsity());

    casadi_int a_ncol = a.size2();
    casadi_int b_ncol = b.size2();
    const casadi_int* a_colind = a.colind();
    //const casadi_int* a_row = a.row();
    const casadi_int* b_colind = b.colind();
    //const casadi_int* b_row = b.row();

    std::vector<Scalar> ret(a.nnz()*b.nnz());
    Scalar* ret_ptr = get_ptr(ret);

    casadi_int k=0;

    const Scalar* v_a = get_ptr(a);
    const Scalar* v_b = get_ptr(b);

    // Loop over the columns
    for (casadi_int a_cc=0; a_cc<a_ncol; ++a_cc) {
      // Loop over the columns
      for (casadi_int b_cc=0; b_cc<b_ncol; ++b_cc) {
        // Loop over existing nonzeros
        for (casadi_int a_el=a_colind[a_cc]; a_el<a_colind[a_cc+1]; ++a_el) {
          Scalar a_v = v_a[a_el];
          // Loop over existing nonzeros
          for (casadi_int b_el=b_colind[b_cc]; b_el<b_colind[b_cc+1]; ++b_el) {
            Scalar b_v = v_b[b_el];
            ret_ptr[k++] = a_v*b_v;
          }
        }
      }
    }
    return Matrix<Scalar>(sp_ret, ret, false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diag(const Matrix<Scalar>& A) {
    // Nonzero mapping
    std::vector<casadi_int> mapping;
    // Get the sparsity
    Sparsity sp = A.sparsity().get_diag(mapping);

    Matrix<Scalar> ret = zeros(sp);

    for (casadi_int k=0; k<mapping.size(); k++) ret.nz(k) = A.nz(mapping[k]);
    return ret;
  }

  /** \brief   Construct a matrix with given block on the diagonal */
  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::diagcat(const std::vector< Matrix<Scalar> > &A) {
    std::vector<Scalar> data;

    std::vector<Sparsity> sp;
    for (casadi_int i=0;i<A.size();++i) {
      data.insert(data.end(), A[i].nonzeros().begin(), A[i].nonzeros().end());
      sp.push_back(A[i].sparsity());
    }

    return Matrix<Scalar>(Sparsity::diagcat(sp), data, false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::unite(const Matrix<Scalar>& A, const Matrix<Scalar>& B) {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = A.sparsity().unite(B.sparsity(), mapping);

    // Create return matrix
    Matrix<Scalar> ret = zeros(sp);

    // Copy sparsity
    casadi_int elA=0, elB=0;
    for (casadi_int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        ret.nonzeros()[k] = A.nonzeros()[elA++];
      } else if (mapping[k]==2) {
        ret.nonzeros()[k] = B.nonzeros()[elB++];
      } else {
        casadi_error("Pattern intersection not empty");
      }
    }

    casadi_assert_dev(A.nnz()==elA);
    casadi_assert_dev(B.nnz()==elB);

    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::polyval(const Matrix<Scalar>& p, const Matrix<Scalar>& x) {
    casadi_assert(p.is_dense(), "polynomial coefficients vector must be dense");
    casadi_assert(p.is_vector() && p.nnz()>0, "polynomial coefficients must be a vector");
    Matrix<Scalar> ret = x;
    for (auto&& e : ret.nonzeros()) {
      e = casadi_polyval(p.ptr(), p.numel()-1, e);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::norm_inf_mul(const Matrix<Scalar>& x,
                                                  const Matrix<Scalar>& y) {
    casadi_assert(y.size1()==x.size2(), "Dimension error. Got " + x.dim()
                          + " times " + y.dim() + ".");

    // Allocate work vectors
    std::vector<Scalar> dwork(x.size1());
    std::vector<casadi_int> iwork(x.size1()+1+y.size2());

    // Call C runtime
    return casadi_norm_inf_mul(x.ptr(), x.sparsity(), y.ptr(), y.sparsity(),
                               get_ptr(dwork), get_ptr(iwork));
  }

  template<typename Scalar>
  void Matrix<Scalar>::expand(const Matrix<Scalar>& ex,
                                Matrix<Scalar> &weights, Matrix<Scalar>& terms) {
    casadi_error("'expand' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_const(const Matrix<Scalar>& ex,
                                              const Matrix<Scalar>& tval,
                                              const Matrix<Scalar>& val) {
    casadi_error("'pw_const' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::pw_lin(const Matrix<Scalar>& ex,
                                            const Matrix<Scalar>& tval,
                                            const Matrix<Scalar>& val) {
    casadi_error("'pw_lin' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::if_else(const Matrix<Scalar> &cond,
                                             const Matrix<Scalar> &if_true,
                                             const Matrix<Scalar> &if_false,
                                             bool short_circuit) {
    return if_else_zero(cond, if_true) + if_else_zero(!cond, if_false);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::conditional(const Matrix<Scalar>& ind,
                                                 const std::vector<Matrix<Scalar> >& x,
                                                 const Matrix<Scalar>& x_default,
                                                 bool short_circuit) {
    casadi_assert(!short_circuit,
      "Short-circuiting 'conditional' not supported for " + type_name());
    casadi_assert(ind.is_scalar(true),
      "conditional: first argument must be scalar. Got " + ind.dim()+ " instead.");

    Matrix<Scalar> ret = x_default;
    for (casadi_int k=0; k<x.size(); ++k) {
      ret = if_else(ind==k, x[k], ret, short_circuit);
    }
    return ret;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::heaviside(const Matrix<Scalar>& x) {
    return (1+sign(x))/2;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::rectangle(const Matrix<Scalar>& x) {
    return 0.5*(sign(x+0.5)-sign(x-0.5));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::triangle(const Matrix<Scalar>& x) {
    return rectangle(x/2)*(1-fabs(x));
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::ramp(const Matrix<Scalar>& x) {
    return x*heaviside(x);
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  gauss_quadrature(const Matrix<Scalar> &f,
                   const Matrix<Scalar> &x, const Matrix<Scalar> &a,
                   const Matrix<Scalar> &b, casadi_int order) {
    return gauss_quadrature(f, x, a, b, order, Matrix<Scalar>());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::gauss_quadrature(const Matrix<Scalar>& f,
                                                      const Matrix<Scalar>& x,
                                                      const Matrix<Scalar>& a,
                                                      const Matrix<Scalar>& b, casadi_int order,
                                                      const Matrix<Scalar>& w) {
    casadi_error("'gauss_quadrature' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::simplify(const Matrix<Scalar> &x) {
    return x;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::substitute(const Matrix<Scalar>& ex,
                                                const Matrix<Scalar>& v,
                                                const Matrix<Scalar>& vdef) {
    casadi_error("'substitute' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> >
  Matrix<Scalar>::substitute(const std::vector<Matrix<Scalar> >& ex,
                               const std::vector<Matrix<Scalar> >& v,
                               const std::vector<Matrix<Scalar> >& vdef) {
    casadi_error("'substitute' not defined for " + type_name());
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::substitute_inplace(const std::vector<Matrix<Scalar> >& v,
                                           std::vector<Matrix<Scalar> >& vdef,
                                           std::vector<Matrix<Scalar> >& ex,
                                           bool reverse) {
    casadi_error("'substitute_inplace' not defined for " + type_name());
  }

  template<typename Scalar>
  bool Matrix<Scalar>::depends_on(const Matrix<Scalar> &x, const Matrix<Scalar> &arg) {
    casadi_error("'depends_on' not defined for " + type_name());
    return false;
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::
  jacobian(const Matrix<Scalar> &f, const Matrix<Scalar> &x, const Dict& opts) {
    casadi_error("'jacobian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x) {
    casadi_error("'hessian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::hessian(const Matrix<Scalar> &f,
                                             const Matrix<Scalar> &x,
                                             Matrix<Scalar> &g) {
    casadi_error("'hessian' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  std::vector<std::vector<Matrix<Scalar> > >
  Matrix<Scalar>::
  forward(const std::vector<Matrix<Scalar> > &ex,
          const std::vector<Matrix<Scalar> > &arg,
          const std::vector<std::vector<Matrix<Scalar> > > &v,
          const Dict& opts) {
    casadi_error("'forward' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<std::vector<Matrix<Scalar> > >
  Matrix<Scalar>::
  reverse(const std::vector<Matrix<Scalar> > &ex,
          const std::vector<Matrix<Scalar> > &arg,
          const std::vector<std::vector<Matrix<Scalar> > > &v,
          const Dict& opts) {
    casadi_error("'reverse' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<bool>
  Matrix<Scalar>::which_depends(const Matrix<Scalar> &expr, const Matrix<Scalar> &var,
      casadi_int order, bool tr) {
    casadi_error("'which_depends' not defined for " + type_name());
    return std::vector<bool>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::taylor(const Matrix<Scalar>& f,
                                            const Matrix<Scalar>& x,
                                            const Matrix<Scalar>& a, casadi_int order) {
    casadi_error("'taylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, casadi_int order) {
    casadi_error("'mtaylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::mtaylor(const Matrix<Scalar>& f,
                                             const Matrix<Scalar>& x,
                                             const Matrix<Scalar>& a, casadi_int order,
                                             const std::vector<casadi_int>&order_contributions) {
    casadi_error("'mtaylor' not defined for " + type_name());
    return Matrix<Scalar>();
  }

  template<typename Scalar>
  casadi_int Matrix<Scalar>::n_nodes(const Matrix<Scalar>& x) {
    casadi_error("'n_nodes' not defined for " + type_name());
    return 0;
  }

  template<typename Scalar>
  std::string
  Matrix<Scalar>::print_operator(const Matrix<Scalar>& x,
                                   const std::vector<std::string>& args) {
    casadi_error("'print_operator' not defined for " + type_name());
    return std::string();
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::symvar(const Matrix<Scalar>& x) {
    casadi_error("'symvar' not defined for " + type_name());
    return std::vector<Matrix<Scalar> >();
  }

  template<typename Scalar>
  void Matrix<Scalar>::shared(std::vector<Matrix<Scalar> >& ex,
                                       std::vector<Matrix<Scalar> >& v,
                                       std::vector<Matrix<Scalar> >& vdef,
                                       const std::string& v_prefix,
                                       const std::string& v_suffix) {
    casadi_error("'shared' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_coeff(const Matrix<Scalar>& f,
                                                const Matrix<Scalar>&x) {
    casadi_error("'poly_coeff' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::poly_roots(const Matrix<Scalar>& p) {
    casadi_error("'poly_roots' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::eig_symbolic(const Matrix<Scalar>& m) {
    casadi_error("'eig_symbolic' not defined for " + type_name());
  }

  template<typename Scalar>
  DM Matrix<Scalar>::evalf(const Matrix<Scalar>& m) {
    Function f("f", std::vector<SX>{}, std::vector<SX>{m});
    return f(std::vector<DM>{})[0];
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::sparsify(const Matrix<Scalar>& x, double tol) {
    // Quick return if there are no entries to be removed
    bool remove_nothing = true;
    for (auto it=x.nonzeros().begin(); it!=x.nonzeros().end() && remove_nothing; ++it) {
      remove_nothing = !casadi_limits<Scalar>::is_almost_zero(*it, tol);
    }
    if (remove_nothing) return x;

    // Get the current sparsity pattern
    casadi_int size1 = x.size1();
    casadi_int size2 = x.size2();
    const casadi_int* colind = x.colind();
    const casadi_int* row = x.row();

    // Construct the new sparsity pattern
    std::vector<casadi_int> new_colind(1, 0), new_row;
    std::vector<Scalar> new_data;

    // Loop over the columns
    for (casadi_int cc=0; cc<size2; ++cc) {
      // Loop over existing nonzeros
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        // If it is not known to be a zero
        if (!casadi_limits<Scalar>::is_almost_zero(x->at(el), tol)) {
          // Save the nonzero in its new location
          new_data.push_back(x->at(el));

          // Add to pattern
          new_row.push_back(row[el]);
        }
      }
      // Save the new column offset
      new_colind.push_back(new_row.size());
    }

    // Construct the sparsity pattern
    Sparsity sp(size1, size2, new_colind, new_row);

    // Construct matrix and return
    return Matrix<Scalar>(sp, new_data);
  }


  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::get_input(const Function& f) {
    casadi_error("'get_input' not defined for " + type_name());
  }

  template<typename Scalar>
  std::vector<Matrix<Scalar> > Matrix<Scalar>::get_free(const Function& f) {
    casadi_error("'get_free' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar>::operator double() const {
    casadi_assert_dev(is_scalar());
    return static_cast<double>(scalar());
  }

  template<typename Scalar>
  Matrix<Scalar>::operator casadi_int() const {
    casadi_assert_dev(is_scalar());
    return static_cast<casadi_int>(scalar());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::_sym(const std::string& name, const Sparsity& sp) {
    casadi_error("'sym' not defined for " + type_name());
  }

  template<typename Scalar>
  Matrix<Scalar> Matrix<Scalar>::rand(const Sparsity& sp) { // NOLINT(runtime/threadsafe_fn)

    casadi_error("'rand' not defined for " + type_name());
  }

  // Template specializations
  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  solve(const Matrix<double>& a, const Matrix<double>& b,
        const std::string& lsolver, const Dict& dict);

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  pinv(const Matrix<double>& A, const std::string& lsolver, const Dict& dict);

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  rand(const Sparsity& sp); // NOLINT(runtime/threadsafe_fn)

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  expm_const(const Matrix<double>& A, const Matrix<double>& t);

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  expm(const Matrix<double>& A);

  bool CASADI_EXPORT is_slice(const IM& x, bool ind1) {
    return x.is_scalar() || (x.is_column() && x.is_dense() && is_slice(x.nonzeros(), ind1));
  }

  Slice CASADI_EXPORT to_slice(const IM& x, bool ind1) {
    return x.is_scalar() ? Slice(x.scalar(), ind1) : to_slice(x.nonzeros(), ind1);
  }

  template<>
  Matrix<double> Matrix<double>::
  solve(const Matrix<double>& A, const Matrix<double>& b,
        const string& lsolver, const Dict& dict) {
    Linsol mysolver("tmp", lsolver, A.sparsity(), dict);
    return mysolver.solve(A, b, false);
  }

  template<>
  Matrix<double> Matrix<double>::
  inv(const Matrix<double>& A,
        const string& lsolver, const Dict& dict) {
    return solve(A, DM::eye(A.size1()), lsolver, dict);
  }

  template<>
  Matrix<double> Matrix<double>::
  pinv(const Matrix<double>& A, const string& lsolver,
       const Dict& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mtimes(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mtimes(A, A.T()), A, lsolver, dict).T();
    }
  }

  template<>
  CASADI_EXPORT Matrix<double> Matrix<double>::
  rand(const Sparsity& sp) { // NOLINT(runtime/threadsafe_fn)
    // C++11 random number generator
    std::uniform_real_distribution<double> distribution(0., 1.);
    // Nonzeros
    std::vector<double> nz(sp.nnz());
    for (double& e : nz) e = distribution(rng_);
    // Construct return object
    return Matrix<double>(sp, nz, false);
  }

  template<>
  Matrix<double> Matrix<double>::
  expm_const(const Matrix<double>& A, const Matrix<double>& t) {
    return expm(A*t);
  }

  template<>
  Matrix<double> Matrix<double>::
  expm(const Matrix<double>& A) {
    Function ret = expmsol("mysolver", "slicot", A.sparsity());
    return ret(std::vector<DM>{A, 1})[0];
  }

  template<>
  bool Matrix<SXElem>::__nonzero__() const {
    casadi_assert(numel()==1,
      "Only scalar Matrix could have a truth value, but you "
      "provided a shape" + dim());
    return nonzeros().at(0).__nonzero__();
  }

  template<>
  void SX::set_max_depth(casadi_int eq_depth) {
    SXNode::eq_depth_ = eq_depth;
  }

  template<>
  casadi_int SX::get_max_depth() {
    return SXNode::eq_depth_;
  }

  template<>
  SX SX::_sym(const string& name, const Sparsity& sp) {
    // Create a dense n-by-m matrix
    vector<SXElem> retv;

    // Check if individial names have been provided
    if (name[0]=='[') {

      // Make a copy of the string and modify it as to remove the special characters
      string modname = name;
      for (string::iterator it=modname.begin(); it!=modname.end(); ++it) {
        switch (*it) {
        case '(': case ')': case '[': case ']': case '{': case '}': case ',': case ';': *it = ' ';
        }
      }

      istringstream iss(modname);
      string varname;

      // Loop over elements
      while (!iss.fail()) {
        // Read the name
        iss >> varname;

        // Append to the return vector
        if (!iss.fail())
          retv.push_back(SXElem::sym(varname));
      }
    } else if (sp.is_scalar(true)) {
      retv.push_back(SXElem::sym(name));
    } else {
      // Scalar
      stringstream ss;
      for (casadi_int k=0; k<sp.nnz(); ++k) {
        ss.str("");
        ss << name << "_" << k;
        retv.push_back(SXElem::sym(ss.str()));
      }
    }

    // Determine dimensions automatically if empty
    if (sp.is_scalar(true)) {
      return SX(retv);
    } else {
      return SX(sp, retv, false);
    }
  }

  template<>
  bool SX::is_regular() const {
    // First pass: ignore symbolics
    for (casadi_int i=0; i<nnz(); ++i) {
      const SXElem& x = nonzeros().at(i);
      if (x.is_constant()) {
        if (x.is_nan() || x.is_inf() || x.is_minus_inf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (casadi_int i=0; i<nnz(); ++i) {
      if (!nonzeros().at(i).is_regular()) return false;
    }
    return true;
  }

  template<>
  bool SX::is_smooth() const {
    // Make a function
    Function temp("temp", {SX()}, {*this});

    // Run the function on the temporary variable
    SXFunction* t = temp.get<SXFunction>();
    return t->is_smooth();
  }

  template<>
  casadi_int SX::element_hash() const {
    return scalar().__hash__();
  }

  template<>
  bool SX::is_leaf() const {
    return scalar().is_leaf();
  }

  template<>
  bool SX::is_commutative() const {
    return scalar().is_commutative();
  }

  template<>
  bool SX::is_symbolic() const {
    if (is_dense()) {
      return is_valid_input();
    } else {
      return false;
    }
  }

  template<>
  casadi_int SX::op() const {
    return scalar().op();
  }

  template<>
  bool SX::is_op(casadi_int op) const {
    return scalar().is_op(op);
  }

  template<>
  bool SX::is_valid_input() const {
    for (casadi_int k=0; k<nnz(); ++k) // loop over non-zero elements
      if (!nonzeros().at(k)->is_symbolic()) // if an element is not symbolic
        return false;

    return true;
  }

  template<> bool SX::has_duplicates() const {
    bool has_duplicates = false;
    for (auto&& i : nonzeros_) {
      bool is_duplicate = i.get_temp()!=0;
      if (is_duplicate) {
        casadi_warning("Duplicate expression: " + str(i));
      }
      has_duplicates = has_duplicates || is_duplicate;
      i.set_temp(1);
    }
    return has_duplicates;
  }

  template<> void SX::reset_input() const {
    for (auto&& i : nonzeros_) {
      i.set_temp(0);
    }
  }

  template<>
  string SX::name() const {
    return scalar().name();
  }

  template<>
  SX SX::dep(casadi_int ch) const {
    return scalar().dep(ch);
  }

  template<>
  casadi_int SX::n_dep() const {
    return scalar().n_dep();
  }

  template<>
  void SX::expand(const SX& ex2, SX& ww, SX& tt) {
    casadi_assert_dev(ex2.is_scalar());
    SXElem ex = ex2.scalar();

    // Terms, weights and indices of the nodes that are already expanded
    vector<vector<SXNode*> > terms;
    vector<vector<double> > weights;
    std::map<SXNode*, casadi_int> indices;

    // Stack of nodes that are not yet expanded
    stack<SXNode*> to_be_expanded;
    to_be_expanded.push(ex.get());

    while (!to_be_expanded.empty()) { // as long as there are nodes to be expanded

      // Check if the last element on the stack is already expanded
      if (indices.find(to_be_expanded.top()) != indices.end()) {
        // Remove from stack
        to_be_expanded.pop();
        continue;
      }

      // Weights and terms
      vector<double> w; // weights
      vector<SXNode*> f; // terms

      if (to_be_expanded.top()->is_constant()) { // constant nodes are seen as multiples of one
        w.push_back(to_be_expanded.top()->to_double());
        f.push_back(casadi_limits<SXElem>::one.get());
      } else if (to_be_expanded.top()->is_symbolic()) {
        // symbolic nodes have weight one and itself as factor
        w.push_back(1);
        f.push_back(to_be_expanded.top());
      } else { // unary or binary node

        casadi_assert_dev(to_be_expanded.top()->n_dep()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        SXNode* node = to_be_expanded.top();
        // If we have a binary node that we can factorize
        if (node->op() == OP_ADD || node->op() == OP_SUB ||
           (node->op() == OP_MUL  && (node->dep(0)->is_constant() ||
                                         node->dep(1)->is_constant()))) {
          // Make sure that both children are factorized, if not - add to stack
          if (indices.find(node->dep(0).get()) == indices.end()) {
            to_be_expanded.push(node->dep(0).get());
            continue;
          }
          if (indices.find(node->dep(1).get()) == indices.end()) {
            to_be_expanded.push(node->dep(1).get());
            continue;
          }

          // Get indices of children
          casadi_int ind1 = indices[node->dep(0).get()];
          casadi_int ind2 = indices[node->dep(1).get()];

          // If multiplication
          if (node->op() == OP_MUL) {
            double fac;
            // Multiplication where the first factor is a constant
            if (node->dep(0)->is_constant()) {
              fac = node->dep(0)->to_double();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = node->dep(1)->to_double();
              f = terms[ind1];
              w = weights[ind1];
            }
            for (casadi_int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if (node->op() == OP_ADD) {          // Addition: join both sums
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];    w.insert(w.end(), weights[ind2].begin(), weights[ind2].end());
            } else {      // Subtraction: join both sums with negative weights for second term
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];
              w.reserve(f.size());
              for (casadi_int i=0; i<weights[ind2].size(); ++i) w.push_back(-weights[ind2][i]);
            }
            // Eliminate multiple elements
            vector<double> w_new; w_new.reserve(w.size());   // weights
            vector<SXNode*> f_new;  f_new.reserve(f.size());   // terms
            std::map<SXNode*, casadi_int> f_ind; // index in f_new

            for (casadi_int i=0; i<w.size(); i++) {
              // Try to locate the node
              auto it = f_ind.find(f[i]);
              if (it == f_ind.end()) { // if the term wasn't found
                w_new.push_back(w[i]);
                f_new.push_back(f[i]);
                f_ind[f[i]] = f_new.size()-1;
              } else { // if the term already exists
                w_new[it->second] += w[i]; // just add the weight
              }
            }
            w = w_new;
            f = f_new;
          }
        } else { // if we have a binary node that we cannot factorize
          // By default,
          w.push_back(1);
          f.push_back(node);

        }
      }

      // Save factorization of the node
      weights.push_back(w);
      terms.push_back(f);
      indices[to_be_expanded.top()] = terms.size()-1;

      // Remove node from stack
      to_be_expanded.pop();
    }

    // Save expansion to output
    casadi_int thisind = indices[ex.get()];
    ww = SX(weights[thisind]);

    vector<SXElem> termsv(terms[thisind].size());
    for (casadi_int i=0; i<termsv.size(); ++i)
      termsv[i] = SXElem::create(terms[thisind][i]);
    tt = SX(termsv);
  }

  template<>
  SX SX::pw_const(const SX& t, const SX& tval, const SX& val) {
    // number of intervals
    casadi_int n = val.numel();

    casadi_assert(t.is_scalar(), "t must be a scalar");
    casadi_assert(tval.numel() == n-1, "dimensions do not match");

    SX ret = val->at(0);
    for (casadi_int i=0; i<n-1; ++i) {
      ret += (val(i+1)-val(i)) * (t>=tval(i));
    }

    return ret;
  }

  template<>
  SX SX::pw_lin(const SX& t, const SX& tval, const SX& val) {
    // Number of points
    casadi_int N = tval.numel();
    casadi_assert(N>=2, "pw_lin: N>=2");
    casadi_assert(val.numel() == N, "dimensions do not match");

    // Gradient for each line segment
    SX g = SX(1, N-1);
    for (casadi_int i=0; i<N-1; ++i)
      g(i) = (val(i+1)- val(i))/(tval(i+1)-tval(i));

    // Line segments
    SX lseg = SX(1, N-1);
    for (casadi_int i=0; i<N-1; ++i)
      lseg(i) = val(i) + g(i)*(t-tval(i));

    // Return piecewise linear function
    return pw_const(t, tval(range(1, N-1)), lseg);
  }

  template<>
  SX SX::gauss_quadrature(const SX& f, const SX& x, const SX& a, const SX& b, casadi_int order,
                          const SX& w) {
    casadi_assert(order == 5, "gauss_quadrature: order must be 5");
    casadi_assert(w.is_empty(), "gauss_quadrature: empty weights");

    // Change variables to [-1, 1]
    if (!is_equal(a.scalar(), -1) || !is_equal(b.scalar(), 1)) {
      SX q1 = (b-a)/2;
      SX q2 = (b+a)/2;

      Function fcn("gauss_quadrature", {x}, {f});

      return q1*gauss_quadrature(fcn(q1*x+q2).at(0), x, -1, 1);
    }

    // Gauss points
    vector<double> xi;
    xi.push_back(-std::sqrt(5 + 2*std::sqrt(10.0/7))/3);
    xi.push_back(-std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
    xi.push_back(0);
    xi.push_back(std::sqrt(5 - 2*std::sqrt(10.0/7))/3);
    xi.push_back(std::sqrt(5 + 2*std::sqrt(10.0/7))/3);

    // Gauss weights
    vector<double> wi;
    wi.push_back((322-13*std::sqrt(70.0))/900.0);
    wi.push_back((322+13*std::sqrt(70.0))/900.0);
    wi.push_back(128/225.0);
    wi.push_back((322+13*std::sqrt(70.0))/900.0);
    wi.push_back((322-13*std::sqrt(70.0))/900.0);

    // Evaluate at the Gauss points
    Function fcn("gauss_quadrature", {x}, {f});
    vector<SXElem> f_val(5);
    for (casadi_int i=0; i<5; ++i)
      f_val[i] = fcn(SX(xi[i])).at(0).scalar();

    // Weighted sum
    SXElem sum;
    for (casadi_int i=0; i<5; ++i)
      sum += wi[i]*f_val[i];

    return sum;
  }

  template<>
  SX SX::simplify(const SX& x) {
    SX r = x;
    for (casadi_int el=0; el<r.nnz(); ++el) {
      // Start by expanding the node to a weighted sum
      SX terms, weights;
      expand(r.nz(el), weights, terms);

      // Make a scalar product to get the simplified expression
      r.nz(el) = mtimes(terms.T(), weights);
    }
    return r;
  }

  template<>
  SX SX::substitute(const SX& ex, const SX& v, const SX& vdef) {
    return substitute(vector<SX>{ex}, vector<SX>{v}, vector<SX>{vdef}).front();
  }

  template<>
  vector<SX>
  SX::substitute(const vector<SX>& ex, const vector<SX>& v, const vector<SX>& vdef) {

    // Assert consistent dimensions
    if (v.size()!=vdef.size()) {
      casadi_warning("subtitute: number of symbols to replace ( " + str(v.size()) + ") "
                     "must match number of expressions (" + str(vdef.size()) + ") "
                     "to replace them with.");
    }

    // Quick return if all equal
    bool all_equal = true;
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !is_equal(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Check sparsities
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].sparsity()!=vdef[k].sparsity()) {
        // Expand vdef to sparsity of v if vdef is scalar
        if (vdef[k].is_scalar() && vdef[k].nnz()==1) {
          vector<SX> vdef_mod = vdef;
          vdef_mod[k] = SX(v[k].sparsity(), vdef[k]->at(0), false);
          return substitute(ex, v, vdef_mod);
        } else {
          casadi_error("Sparsities of v and vdef must match. Got v: "
                       + v[k].dim() + " and vdef: " + vdef[k].dim() + ".");
        }
      }
    }


    // Otherwise, evaluate symbolically
    Function F("tmp", v, ex);
    return F(vdef);
  }

  template<>
  void SX::substitute_inplace(const vector<SX >& v, vector<SX >& vdef,
                             vector<SX >& ex, bool reverse) {
    // Assert correctness
    casadi_assert_dev(v.size()==vdef.size());
    for (casadi_int i=0; i<v.size(); ++i) {
      casadi_assert(v[i].is_symbolic(), "the variable is not symbolic");
      casadi_assert(v[i].sparsity() == vdef[i].sparsity(), "the sparsity patterns of the "
                            "expression and its defining bexpression do not match");
    }

    // Quick return if empty or single expression
    if (v.empty()) return;

    // Function inputs
    vector<SX> f_in;
    if (!reverse) f_in.insert(f_in.end(), v.begin(), v.end());

    // Function outputs
    vector<SX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    Function f("tmp", f_in, f_out);

    // Get references to the internal data structures
    SXFunction *ff = f.get<SXFunction>();
    const vector<ScalarAtomic>& algorithm = ff->algorithm_;
    vector<SXElem> work(f.sz_w());

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=ff->operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = ff->constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = ff->free_vars_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        // reverse is false, substitute out
        work[it->i0] = vdef.at(it->i1)->at(it->i2);
        break;
      case OP_OUTPUT:
        if (it->i0 < v.size()) {
          vdef.at(it->i0)->at(it->i2) = work[it->i1];
          if (reverse) {
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i0)->at(it->i2);
          }
        } else {
          // Auxillary output
          ex.at(it->i0 - v.size())->at(it->i2) = work[it->i1];
        }
        break;
      case OP_CONST:      work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }

          // Avoid creating duplicates
          const casadi_int depth = 2; // NOTE: a higher depth could possibly give more savings
          work[it->i0].assignIfDuplicate(*b_it++, depth);
        }
      }
    }
  }

  template<>
  bool SX::depends_on(const SX &x, const SX &arg) {
    if (x.nnz()==0) return false;

    // Construct a temporary algorithm
    Function temp("temp", {arg}, {x});

    // Perform a single dependency sweep
    vector<bvec_t> t_in(arg.nnz(), 1), t_out(x.nnz());
    temp({get_ptr(t_in)}, {get_ptr(t_out)});

    // Loop over results
    for (casadi_int i=0; i<t_out.size(); ++i) {
      if (t_out[i]) return true;
    }

    return false;
  }

  template<>
  SX SX::jacobian(const SX &f, const SX &x, const Dict& opts) {
    // Propagate verbose option to helper function
    Dict h_opts;
    if (opts.count("verbose")) h_opts["verbose"] = opts.at("verbose");
    Function h("jac_helper", {x}, {f}, h_opts);
    return h.get<SXFunction>()->jac(0, 0, opts);
  }

  template<>
  SX SX::hessian(const SX &ex, const SX &arg) {
    SX g;
    return hessian(ex, arg, g);
  }

  template<>
  SX SX::hessian(const SX &ex, const SX &arg, SX &g) {
    g = gradient(ex, arg);
    return jacobian(g, arg, {{"symmetric", true}});
  }

  template<>
  std::vector<std::vector<SX> >
  SX::forward(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts) {
    // Read options
    bool always_inline = false;
    bool never_inline = false;
    for (auto&& op : opts) {
      if (op.first=="always_inline") {
        always_inline = op.second;
      } else if (op.first=="never_inline") {
        never_inline = op.second;
      } else if (op.first=="verbose") {
        continue;
      }  else {
        casadi_error("No such option: " + string(op.first));
      }
    }
    // Call internal function on a temporary object
    Function temp("forward_temp", arg, ex);
    std::vector<std::vector<SX> > ret;
    temp->call_forward(arg, ex, v, ret, always_inline, never_inline);
    return ret;
  }

  template<>
  std::vector<std::vector<SX> >
  SX::reverse(const std::vector<SX> &ex, const std::vector<SX> &arg,
          const std::vector<std::vector<SX> > &v, const Dict& opts) {
    // Read options
    bool always_inline = false;
    bool never_inline = false;
    for (auto&& op : opts) {
      if (op.first=="always_inline") {
        always_inline = op.second;
      } else if (op.first=="never_inline") {
        never_inline = op.second;
      } else if (op.first=="verbose") {
        continue;
      } else {
        casadi_error("No such option: " + string(op.first));
      }
    }
    // Call internal function on a temporary object
    Function temp("reverse_temp", arg, ex);
    std::vector<std::vector<SX> > ret;
    temp->call_reverse(arg, ex, v, ret, always_inline, never_inline);
    return ret;
  }

  template<>
  std::vector<bool> SX::which_depends(const SX &expr, const SX &var, casadi_int order, bool tr) {
    return _which_depends(expr, var, order, tr);
  }

  template<>
  SX SX::taylor(const SX& f, const SX& x,
                const SX& a, casadi_int order) {
    casadi_assert_dev(x.is_scalar() && a.is_scalar());
    if (f.nnz()!=f.numel())
      throw CasadiException("taylor: not implemented for sparse matrices");
    SX ff = vec(f.T());

    SX result = substitute(ff, x, a);
    double nf=1;
    SX dx = (x-a);
    SX dxa = (x-a);
    for (casadi_int i=1; i<=order; i++) {
      ff = jacobian(ff, x);
      nf*=static_cast<double>(i);
      result+=1/nf * substitute(ff, x, a) * dxa;
      dxa*=dx;
    }
    return reshape(result, f.size2(), f.size1()).T();
  }

  template<>
  SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order) {
    return mtaylor(f, x, a, order, vector<casadi_int>(x.nnz(), 1));
  }

  SX mtaylor_recursive(const SX& ex, const SX& x, const SX& a, casadi_int order,
                       const vector<casadi_int>&order_contributions,
                       const SXElem & current_dx=casadi_limits<SXElem>::one,
                       double current_denom=1, casadi_int current_order=1) {
    SX result = substitute(ex, x, a)*current_dx/current_denom;
    for (casadi_int i=0;i<x.nnz();i++) {
      if (order_contributions[i]<=order) {
        result += mtaylor_recursive(SX::jacobian(ex, x->at(i)),
                                    x, a,
                                    order-order_contributions[i],
                                    order_contributions,
                                    current_dx*(x->at(i)-a->at(i)),
                                    current_denom*static_cast<double>(current_order),
                                    current_order+1);
      }
    }
    return result;
  }

  template<>
  SX SX::mtaylor(const SX& f, const SX& x, const SX& a, casadi_int order,
                 const vector<casadi_int>& order_contributions) {
    casadi_assert(f.nnz()==f.numel() && x.nnz()==x.numel(),
                          "mtaylor: not implemented for sparse matrices");

    casadi_assert(x.nnz()==order_contributions.size(),
                          "mtaylor: number of non-zero elements in x (" + str(x.nnz())
                          + ") must match size of order_contributions ("
                          + str(order_contributions.size()) + ")");

    return reshape(mtaylor_recursive(vec(f), x, a, order,
                                     order_contributions),
                   f.size2(), f.size1()).T();
  }

  template<>
  casadi_int SX::n_nodes(const SX& x) {
    Function f("tmp", {SX()}, {x});
    return f.n_nodes();
  }

  template<>
  string
  SX::print_operator(const SX& X, const vector<string>& args) {
    SXElem x = X.scalar();
    casadi_int ndeps = casadi_math<double>::ndeps(x.op());
    casadi_assert(ndeps==1 || ndeps==2, "Not a unary or binary operator");
    casadi_assert(args.size()==ndeps, "Wrong number of arguments");
    if (ndeps==1) {
      return casadi_math<double>::print(x.op(), args.at(0));
    } else {
      return casadi_math<double>::print(x.op(), args.at(0), args.at(1));
    }
  }

  template<>
  vector<SX> SX::symvar(const SX& x) {
    Function f("tmp", vector<SX>{}, {x});
    return f.free_sx();
  }

  template<>
  void SX::shared(vector<SX >& ex,
                         vector<SX >& v_sx,
                         vector<SX >& vdef_sx,
                         const string& v_prefix,
                         const string& v_suffix) {

    // Sort the expression
    Function f("tmp", vector<SX>(), ex);
    SXFunction *ff = f.get<SXFunction>();

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = ff->algorithm_;
    vector<SXElem> work(f.sz_w());
    vector<SXElem> work2 = work;

    // Iterator to the binary operations
    vector<SXElem>::const_iterator b_it=ff->operations_.begin();

    // Iterator to stack of constants
    vector<SXElem>::const_iterator c_it = ff->constants_.begin();

    // Iterator to free variables
    vector<SXElem>::const_iterator p_it = ff->free_vars_.begin();

    // Count how many times an expression has been used
    vector<casadi_int> usecount(work.size(), 0);

    // Evaluate the algorithm
    vector<SXElem> v, vdef;
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      // Increase usage counters
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
          if (usecount[it->i2]==0) {
            usecount[it->i2]=1;
          } else if (usecount[it->i2]==1) {
            // Get a suitable name
            vdef.push_back(work[it->i2]);
            usecount[it->i2]=-1; // Extracted, do not extract again
          }
        // fall-through
      case OP_OUTPUT:
      default: // Unary operation, binary operation or output
        if (usecount[it->i1]==0) {
          usecount[it->i1]=1;
        } else if (usecount[it->i1]==1) {
          vdef.push_back(work[it->i1]);
          usecount[it->i1]=-1; // Extracted, do not extract again
        }
      }

      // Perform the operation
      switch (it->op) {
      case OP_OUTPUT:
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->i0] = -1; // Never extract since it is a primitive type
        break;
      default:
        work[it->i0] = *b_it++;
        usecount[it->i0] = 0; // Not (yet) extracted
        break;
      }
    }

    // Create intermediate variables
    stringstream v_name;
    for (casadi_int i=0; i<vdef.size(); ++i) {
      v_name.str(string());
      v_name << v_prefix << i << v_suffix;
      v.push_back(SXElem::sym(v_name.str()));
    }

    casadi_assert(vdef.size()<numeric_limits<int>::max(), "Integer overflow");

    // Mark the above expressions
    for (casadi_int i=0; i<vdef.size(); ++i) {
      vdef[i].set_temp(static_cast<int>(i)+1);
    }

    // Save the marked nodes for later cleanup
    vector<SXElem> marked = vdef;

    // Reset iterator
    b_it=ff->operations_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_OUTPUT:     ex.at(it->i0)->at(it->i2) = work[it->i1];      break;
      case OP_CONST:      work2[it->i0] = work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work2[it->i0] = work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }
          work2[it->i0] = *b_it++;

          // Replace with intermediate variables
          casadi_int ind = work2[it->i0].get_temp()-1;
          if (ind>=0) {
            vdef.at(ind) = work[it->i0];
            work[it->i0] = v.at(ind);
          }
        }
      }
    }

    // Unmark the expressions
    for (vector<SXElem>::iterator it=marked.begin(); it!=marked.end(); ++it) {
      it->set_temp(0);
    }

    // Save v, vdef
    v_sx.resize(v.size());
    copy(v.begin(), v.end(), v_sx.begin());
    vdef_sx.resize(vdef.size());
    copy(vdef.begin(), vdef.end(), vdef_sx.begin());
  }

  template<>
  SX SX::poly_coeff(const SX& ex, const SX& x) {
    casadi_assert_dev(ex.is_scalar());
    casadi_assert_dev(x.is_scalar());
    casadi_assert_dev(x.is_symbolic());

    vector<SXElem> r;

    SX j = ex;
    casadi_int mult = 1;
    bool success = false;
    for (casadi_int i=0; i<1000; ++i) {
      r.push_back((substitute(j, x, 0)/static_cast<double>(mult)).scalar());
      j = jacobian(j, x);
      if (j.nnz()==0) {
        success = true;
        break;
      }
      mult*=i+1;
    }

    if (!success) casadi_error("poly: supplied expression does not appear to be polynomial.");

    std::reverse(r.begin(), r.end());

    return r;
  }

  template<>
  SX SX::poly_roots(const SX& p) {
    casadi_assert(p.size2()==1,
                          "poly_root(): supplied parameter must be column vector but got "
                          + p.dim() + ".");
    casadi_assert_dev(p.is_dense());
    if (p.size1()==2) { // a*x + b
      SX a = p(0);
      SX b = p(1);
      return -b/a;
    } else if (p.size1()==3) { // a*x^2 + b*x + c
      SX a = p(0);
      SX b = p(1);
      SX c = p(2);
      SX ds = sqrt(b*b-4*a*c);
      SX bm = -b;
      SX a2 = 2*a;
      SX ret = SX::vertcat({(bm-ds)/a2, (bm+ds)/a2});
      return ret;
    } else if (p.size1()==4) {
      // www.cs.iastate.edu/~cs577/handouts/polyroots.pdf
      SX ai = 1/p(0);

      SX p_ = p(1)*ai;
      SX q  = p(2)*ai;
      SX r  = p(3)*ai;

      SX pp = p_*p_;

      SX a = q - pp/3;
      SX b = r + 2.0/27*pp*p_-p_*q/3;

      SX a3 = a/3;

      SX phi = acos(-b/2/sqrt(-a3*a3*a3));

      SX ret = SX::vertcat({cos(phi/3), cos((phi+2*pi)/3), cos((phi+4*pi)/3)});
      ret*= 2*sqrt(-a3);

      ret-= p_/3;
      return ret;
    } else if (p.size1()==5) {
      SX ai = 1/p(0);
      SX b = p(1)*ai;
      SX c = p(2)*ai;
      SX d = p(3)*ai;
      SX e = p(4)*ai;

      SX bb= b*b;
      SX f = c - (3*bb/8);
      SX g = d + (bb*b / 8) - b*c/2;
      SX h = e - (3*bb*bb/256) + (bb * c/16) - (b*d/4);
      SX poly = SX::vertcat({1, f/2, ((f*f -4*h)/16), -g*g/64});
      SX y = poly_roots(poly);

      SX r0 = y(0);
      SX r1 = y(2);

      SX p = sqrt(r0); // two non-zero-roots
      SX q = sqrt(r1);

      SX r = -g/(8*p*q);

      SX s = b/4;

      SX ret = SX::vertcat({
          p + q + r -s,
            p - q - r -s,
            -p + q - r -s,
            -p - q + r -s});
      return ret;
    } else if (is_equal(p(p.nnz()-1)->at(0), 0)) {
      SX ret = SX::vertcat({poly_roots(p(range(p.nnz()-1))), 0});
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. "
                   "Got order " + str(p.size1()-1) + ".");
    }

  }

  template<>
  SX SX::eig_symbolic(const SX& m) {
    casadi_assert(m.size1()==m.size2(), "eig(): supplied matrix must be square");

    vector<SX> ret;

    /// Bring m in block diagonal form, calculating eigenvalues of each block separately
    vector<casadi_int> offset;
    vector<casadi_int> index;
    casadi_int nb = m.sparsity().scc(offset, index);

    SX m_perm = m(offset, offset);

    SX l = SX::sym("l");

    for (casadi_int k=0; k<nb; ++k) {
      vector<casadi_int> r = range(index.at(k), index.at(k+1));
      // det(lambda*I-m) = 0
      ret.push_back(poly_roots(poly_coeff(det(SX::eye(r.size())*l-m_perm(r, r)), l)));
    }

    return vertcat(ret);
  }

  template<>
  void SX::print_split(vector<string>& nz,
                      vector<string>& inter) const {
    // Find out which noded can be inlined
    std::map<const SXNode*, casadi_int> nodeind;
    for (auto&& i : nonzeros_) i->can_inline(nodeind);

    // Print expression
    nz.resize(0);
    nz.reserve(nnz());
    inter.resize(0);
    for (auto&& i : nonzeros_) nz.push_back(i->print_compact(nodeind, inter));
  }

  template<> vector<SX> SX::get_input(const Function& f) {
    return f.sx_in();
  }

  template<> vector<SX> SX::get_free(const Function& f) {
    return f.free_sx();
  }

  template<> void DM::export_code(const std::string& lang,
       std::ostream &stream, const Dict& options) const {

    casadi_assert(lang=="matlab", "Only matlab language supported for now.");

    // Default values for options
    bool opt_inline = false;
    std::string name = "m";
    casadi_int indent_level = 0;
    bool spoof_zero = false;

    // Read options
    for (auto&& op : options) {
      if (op.first=="inline") {
        opt_inline = op.second;
      } else if (op.first=="name") {
        name = op.second.to_string();
      } else if (op.first=="indent_level") {
        indent_level = op.second;
      } else if (op.first=="spoof_zero") {
        spoof_zero = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Construct indent string
    std::string indent = "";
    for (casadi_int i=0;i<indent_level;++i) {
      indent += "  ";
    }

    casadi_assert(!opt_inline, "Inline not supported for now.");

    // Prepare stream for emitting full precision
    std::ios_base::fmtflags fmtfl = stream.flags();
    stream << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    // Obtain nonzeros of matrix
    std::vector<double> d = nonzeros();

    // Spoof numericals
    if (spoof_zero) {
      for (double& e : d) {
        if (e==0) e=1e-200;
      }
    }

    // Short-circuit for (dense) scalars
    if (is_scalar(true)) {
      stream << indent << name << " = " << d[0] << ";" << std::endl;
      stream.flags(fmtfl);
      return;
    }

    // Are all nonzeros equal?
    bool all_equal = true;
    for (double e : d) {
      if (e!=d[0]) {
        all_equal = false;
        break;
      }
    }

    if (all_equal && d.size()>0) {
      // No need to export all individual nonzeros if they are all equal
      stream << indent << name << "_nz = ones(1, " << d.size() << ")*" << d[0] << ";" << std::endl;
    } else {
      // Export nonzeros
      stream << indent << name << "_nz = [";
      for (casadi_int i=0;i<d.size();++i) {
        stream << d[i] << " ";
        if ((i+1)%20 == 0) stream << "..." << std::endl << indent << "  ";
      }
      stream << "];" << std::endl;
    }

    // Reset stream properties
    stream.flags(fmtfl);

    // Cast nonzeros in correct shape
    if (is_dense()) {
      // Special case for dense (for readibility of exported code)
      stream << indent << name << " = reshape(";
      stream << name << "_nz, ";
      stream << size1() << ", " << size2() << ");" << endl;
    } else {
      // For sparse matrices, export Sparsity and use sparse constructor
      Dict opts;
      opts["as_matrix"] = false;
      opts["indent_level"] = indent_level;
      opts["name"] = name;
      opts["indent_level"] = opt_inline;
      sparsity().export_code(lang, stream, opts);
      stream << indent << name << " = sparse(" << name << "_i, " << name << "_j, ";
      stream << name << "_nz, ";
      stream << size1() << ", " << size2() << ");" << endl;
    }
  }

  template<>
  Dict DM::info() const {
    return {{"sparsity", sparsity().info()}, {"data", nonzeros()}};
  }

  template<>
  DM DM::from_info(const Dict& info) {
    Sparsity sp = Sparsity::from_info(info.at("sparsity"));
    std::vector<double> data = info.at("data");
    return DM(sp, data);
  }

  template<>
  Dict IM::info() const {
    return {{"sparsity", sparsity().info()}, {"data", nonzeros()}};
  }

  template<>
  IM IM::from_info(const Dict& info) {
    Sparsity sp = Sparsity::from_info(info.at("sparsity"));
    std::vector<casadi_int> data = info.at("data");
    return IM(sp, data);
  }

  template<>
  Dict SX::info() const {
    return {{"function", Function("f", std::vector<SX>{}, std::vector<SX>{*this})}};
  }

  template<>
  SX SX::from_info(const Dict& info) {
    casadi_error("Not implemented");
    return SX();
  }

  template<>
  void DM::to_file(const std::string& filename, const std::string& format_hint) const {
    std::string format = format_hint;
    if (format_hint=="") {
      std::string extension = filename.substr(filename.rfind(".")+1);
      if (extension=="mtx") {
        format = "mtx";
      } else {
        casadi_error("Could not detect format from extension '" + extension + "'");
      }
    }
    std::ofstream out(filename);
    if (format=="mtx") {
      out << std::scientific << std::setprecision(15);
      out << "%%MatrixMarket matrix coordinate real general" << std::endl;
      out << size1() << " " << size2() << " " << nnz() << std::endl;
      std::vector<casadi_int> row = sparsity().get_row();
      std::vector<casadi_int> col = sparsity().get_col();

      for (casadi_int k=0;k<row.size();++k) {
        out << row[k]+1 << " " << col[k]+1 << " " << nonzeros_[k] << std::endl;
      }
    } else {
      casadi_error("Unknown format '" + format + "'");
    }
  }

  template<>
  void SX::to_file(const std::string& filename, const std::string& format_hint) const {
    casadi_error("Not implemented");
  }

  template<>
  void IM::to_file(const std::string& filename, const std::string& format_hint) const {
    casadi_error("Not implemented");
  }

  // Instantiate templates
  template class casadi_limits<double>;
  template class casadi_limits<casadi_int>;
  template class Matrix<double>;
  template class Matrix<casadi_int>;
  template class Matrix< SXElem >;
} // namespace casadi
