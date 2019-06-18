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


#include "csparse_interface.hpp"
#include "casadi/core/global_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_CSPARSE_EXPORT
  casadi_register_linsol_csparse(LinsolInternal::Plugin* plugin) {
    plugin->creator = CsparseInterface::creator;
    plugin->name = "csparse";
    plugin->doc = CsparseInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &CsparseInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_CSPARSE_EXPORT casadi_load_linsol_csparse() {
    LinsolInternal::registerPlugin(casadi_register_linsol_csparse);
  }

  CsparseInterface::CsparseInterface(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  CsparseInterface::~CsparseInterface() {
    clear_mem();
  }

  CsparseMemory::~CsparseMemory() {
    if (this->S) cs_sfree(this->S);
    if (this->N) cs_nfree(this->N);
  }

  void CsparseInterface::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);
  }

  int CsparseInterface::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<CsparseMemory*>(mem);

    m->N = nullptr;
    m->S = nullptr;
    m->A.nzmax = this->nnz();  // maximum number of entries
    m->A.m = this->nrow(); // number of rows
    m->A.n = this->ncol(); // number of columns
    m->colind.resize(this->ncol()+1);
    m->row.resize(this->nnz());
    copy_vector(this->colind(), m->colind);
    copy_vector(this->row(), m->row);
    m->A.p = get_ptr(m->colind); // row pointers (size n+1)
    m->A.i = get_ptr(m->row); // row pointers (size n+1)
    m->A.x = nullptr; // numerical values, size nzmax
    m->A.nz = -1; // of entries in triplet matrix, -1 for compressed-column

    // Temporary
    m->temp_.resize(m->A.n);
    return 0;
  }

  int CsparseInterface::sfact(void* mem, const double* A) const {
    auto m = static_cast<CsparseMemory*>(mem);

    // Set the nonzeros of the matrix
    m->A.x = const_cast<double*>(A);

    // ordering and symbolic analysis
    casadi_int order = 0; // ordering?
    if (m->S) cs_sfree(m->S);
    m->S = cs_sqr(order, &m->A, 0);
    return 0;
  }

  int CsparseInterface::nfact(void* mem, const double* A) const {
    auto m = static_cast<CsparseMemory*>(mem);

    // Set the nonzeros of the matrix
    m->A.x = const_cast<double*>(A);

    // Make sure that all entries of the linear system are valid
    for (casadi_int k=0; k<this->nnz(); ++k) {
      casadi_assert(!isnan(A[k]),
        "Nonzero " + str(k) + " is not-a-number");
      casadi_assert(!isinf(A[k]),
        "Nonzero " + str(k) + " is infinite");
    }

    if (verbose_) {
      uout() << "CsparseInterface::prepare: numeric factorization" << endl;
      uout() << "linear system to be factorized = " << endl;
      DM(sp_, vector<double>(A, A+nnz())).print_sparse(uout());
    }

    double tol = 1e-8;

    if (m->N) cs_nfree(m->N);
    m->N = cs_lu(&m->A, m->S, tol) ;                 // numeric LU factorization
    if (m->N==nullptr) {
      DM temp(sp_, vector<double>(A, A+nnz()));
      temp = sparsify(temp);
      if (temp.sparsity().is_singular()) {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed due to matrix"
          " being singular. Matrix contains numerical zeros which are "
            "structurally non-zero. Promoting these zeros to be structural "
            "zeros, the matrix was found to be structurally rank deficient."
            " sprank: " << sprank(temp.sparsity()) << " <-> " << temp.size2() << endl;
        if (verbose_) {
          ss << "Sparsity of the linear system: " << endl;
          sp_.disp(ss, true); // print detailed
        }
        throw CasadiException(ss.str());
      } else {
        stringstream ss;
        ss << "CsparseInterface::prepare: factorization failed, check if Jacobian is singular"
           << endl;
        if (verbose_) {
          ss << "Sparsity of the linear system: " << endl;
          sp_.disp(ss, true); // print detailed
        }
        throw CasadiException(ss.str());
      }
    }
    casadi_assert_dev(m->N!=nullptr);
    return 0;
  }

  int CsparseInterface::solve(void* mem, const double* A, double* x,
      casadi_int nrhs, bool tr) const {
    auto m = static_cast<CsparseMemory*>(mem);
    casadi_assert_dev(m->N!=nullptr);

    double *t = &m->temp_.front();

    for (casadi_int k=0; k<nrhs; ++k) {
      if (tr) {
        cs_pvec(m->S->q, x, t, m->A.n) ;       // t = P2*b
        casadi_assert_dev(m->N->U!=nullptr);
        cs_utsolve(m->N->U, t) ;              // t = U'\t
        cs_ltsolve(m->N->L, t) ;              // t = L'\t
        cs_pvec(m->N->pinv, t, x, m->A.n) ;    // x = P1*t
      } else {
        cs_ipvec(m->N->pinv, x, t, m->A.n) ;   // t = P1\b
        cs_lsolve(m->N->L, t) ;               // t = L\t
        cs_usolve(m->N->U, t) ;               // t = U\t
        cs_ipvec(m->S->q, t, x, m->A.n) ;      // x = P2\t
      }
      x += ncol();
    }
    return 0;
  }

} // namespace casadi
