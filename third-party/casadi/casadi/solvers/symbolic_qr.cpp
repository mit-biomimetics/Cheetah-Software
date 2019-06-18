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


#include "symbolic_qr.hpp"

#ifdef WITH_DL
#include <cstdlib>
#endif // WITH_DL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_SYMBOLICQR_EXPORT
  casadi_register_linsol_symbolicqr(LinsolInternal::Plugin* plugin) {
    plugin->creator = SymbolicQr::creator;
    plugin->name = "symbolicqr";
    plugin->doc = SymbolicQr::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SymbolicQr::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_SYMBOLICQR_EXPORT casadi_load_linsol_symbolicqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_symbolicqr);
  }

  SymbolicQr::SymbolicQr(const std::string& name, const Sparsity& sp) :
    LinsolInternal(name, sp) {
  }

  SymbolicQr::~SymbolicQr() {
    clear_mem();
  }

  Options SymbolicQr::options_
  = {{&FunctionInternal::options_},
    {{"fopts",
      {OT_DICT,
       "Options to be passed to generated function objects"}}
     }
  };

  void SymbolicQr::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="fopts") {
        fopts_ = op.second;
      }
    }

    // Symbolic expression for A
    SX A = SX::sym("A", sp_);

    // BTF factorization
    vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    sp_.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Get the inverted column permutation
    std::vector<casadi_int> inv_colperm(colperm.size());
    for (casadi_int k=0; k<colperm.size(); ++k)
      inv_colperm[colperm[k]] = k;

    // Get the inverted row permutation
    std::vector<casadi_int> inv_rowperm(rowperm.size());
    for (casadi_int k=0; k<rowperm.size(); ++k)
      inv_rowperm[rowperm[k]] = k;

    // Permute the linear system
    SX Aperm = A(rowperm, colperm);

    // Generate the QR factorization function
    SX Q1, R1;
    qr(Aperm, Q1, R1);
    factorize_ = Function("QR_fact", {A}, {Q1, R1}, fopts_);

    // Symbolic expressions for solve function
    SX Q = SX::sym("Q", Q1.sparsity());
    SX R = SX::sym("R", R1.sparsity());
    SX b = SX::sym("b", sp_.size2(), 1);

    // Solve non-transposed
    // We have Pb' * Q * R * Px * x = b <=> x = Px' * inv(R) * Q' * Pb * b

    // Permute the right hand sides
    SX bperm = b(rowperm, Slice());

    // Solve the factorized system
    SX xperm = SX::solve(R, mtimes(Q.T(), bperm));

    // Permute back the solution
    SX x = xperm(inv_colperm, Slice());

    // Generate the QR solve function
    vector<SX> solv_in = {Q, R, b};
    solve_ = Function("QR_solv", solv_in, {x}, fopts_);

    // Solve transposed
    // We have (Pb' * Q * R * Px)' * x = b
    // <=> Px' * R' * Q' * Pb * x = b
    // <=> x = Pb' * Q * inv(R') * Px * b

    // Permute the right hand side
    bperm = b(colperm, Slice());

    // Solve the factorized system
    xperm = mtimes(Q, SX::solve(R.T(), bperm));

    // Permute back the solution
    x = xperm(inv_rowperm, Slice());

    // Mofify the QR solve function
    solveT_ = Function("QR_solv_T", solv_in, {x}, fopts_);
  }

  int SymbolicQr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<SymbolicQrMemory*>(mem);

    m->alloc(solveT_);
    m->alloc(solve_);
    m->alloc(factorize_);

    // Temporary storage
    m->w.resize(m->w.size() + sp_.size1());

    // Allocate storage for QR factorization
    m->q.resize(factorize_.nnz_out(0));
    m->r.resize(factorize_.nnz_out(1));
    return 0;
  }

  int SymbolicQr::nfact(void* mem, const double* A) const {
    auto m = static_cast<SymbolicQrMemory*>(mem);

    // Factorize
    fill_n(get_ptr(m->arg), factorize_.n_in(), nullptr);
    m->arg[0] = A;
    fill_n(get_ptr(m->res), factorize_.n_out(), nullptr);
    m->res[0] = get_ptr(m->q);
    m->res[1] = get_ptr(m->r);
    if (factorize_(get_ptr(m->arg), get_ptr(m->res), get_ptr(m->iw), get_ptr(m->w))) return 1;
    return 0;
  }

  int SymbolicQr::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<SymbolicQrMemory*>(mem);

    // Select solve function
    const Function& solv = tr ? solveT_ : solve_;

    // Solve for all right hand sides
    fill_n(get_ptr(m->arg), solv.n_in(), nullptr);
    m->arg[0] = get_ptr(m->q);
    m->arg[1] = get_ptr(m->r);
    fill_n(get_ptr(m->res), solv.n_out(), nullptr);
    for (casadi_int i=0; i<nrhs; ++i) {
      copy_n(x, nrow(), get_ptr(m->w)); // Copy x to a temporary
      m->arg[2] = get_ptr(m->w);
      m->res[0] = x;
      if (solv(get_ptr(m->arg), get_ptr(m->res),
               get_ptr(m->iw), get_ptr(m->w)+nrow(), 0)) return 1;
      x += nrow();
    }
    return 0;
  }

  void SymbolicQr::linsol_eval_sx(const SXElem** arg, SXElem** res,
                                  casadi_int* iw, SXElem* w, void* mem,
                                  bool tr, casadi_int nrhs) const {
    //auto m = static_cast<SymbolicQrMemory*>(mem);
    casadi_assert_dev(arg[0]!=nullptr);
    casadi_assert_dev(arg[1]!=nullptr);
    casadi_assert_dev(res[0]!=nullptr);

    // Get A and factorize it
    SX A = SX::zeros(sp_);
    copy(arg[1], arg[1]+A.nnz(), A->begin());
    vector<SX> v = factorize_(A);

    // Select solve function
    const Function& solv = tr ? solveT_ : solve_;

    // Solve for every right hand side
    v.push_back(SX::zeros(A.size1()));
    const SXElem* a=arg[0];
    SXElem* r=res[0];
    for (casadi_int i=0; i<nrhs; ++i) {
      copy(a, a+v[2].nnz(), v[2]->begin());
      SX rr = solv(v).at(0);
      copy(rr->begin(), rr->end(), r);
      r += rr.nnz();
    }
  }

  void SymbolicQrMemory::alloc(const Function& f) {
    arg.resize(max(arg.size(), f.sz_arg()));
    res.resize(max(res.size(), f.sz_res()));
    iw.resize(max(iw.size(), f.sz_iw()));
    w.resize(max(w.size(), f.sz_w()));
  }

} // namespace casadi
