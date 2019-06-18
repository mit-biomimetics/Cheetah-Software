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


#ifndef CASADI_SOLVE_IMPL_HPP
#define CASADI_SOLVE_IMPL_HPP

#include "solve.hpp"
#include "linsol_internal.hpp"

using namespace std;

namespace casadi {

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A, const Linsol& linear_solver) :
      linsol_(linear_solver) {
    casadi_assert(r.size1() == A.size2(),
      "Solve::Solve: dimension mismatch. Got r " + r.dim() + " and A " + A.dim());
    set_dep(r, A);
    set_sparsity(r.sparsity());
  }

  template<bool Tr>
  std::string Solve<Tr>::disp(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    ss << "(" << arg.at(1);
    if (Tr) ss << "'";
    ss << "\\" << arg.at(0) << ")";
    return ss.str();
  }

  template<bool Tr>
  int Solve<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
    scoped_checkout<Linsol> mem(linsol_);
    if (linsol_.sfact(arg[1], mem)) return 1;
    if (linsol_.nfact(arg[1], mem)) return 1;
    if (linsol_.solve(arg[1], res[0], dep(0).size2(), Tr, mem)) return 1;
    return 0;
  }

  template<bool Tr>
  int Solve<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    linsol_->linsol_eval_sx(arg, res, iw, w, linsol_->memory(0), Tr, dep(0).size2());
    return 0;
  }

  template<bool Tr>
  void Solve<Tr>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    if (arg[0].is_zero()) {
      res[0] = MX(arg[0].size());
    } else {
      res[0] = linsol_.solve(arg[1], arg[0], Tr);
    }
  }

  template<bool Tr>
  void Solve<Tr>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(n_dep());
    for (casadi_int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (casadi_int i=0; i<res.size(); ++i) res[i] = get_output(i);

    // Number of derivatives
    casadi_int nfwd = fseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nfwd);
    std::vector<casadi_int> col_offset(nfwd+1, 0);
    for (casadi_int d=0; d<nfwd; ++d) {
      const MX& B_hat = fseed[d][0];
      const MX& A_hat = fseed[d][1];
      rhs[d] = Tr ? B_hat - mtimes(A_hat.T(), X) : B_hat - mtimes(A_hat, X);
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_.solve(A, horzcat(rhs), Tr), col_offset);

    // Fetch result
    fsens.resize(nfwd);
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d].resize(1);
      fsens[d][0] = rhs[d];
    }
  }

  template<bool Tr>
  void Solve<Tr>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    // Nondifferentiated inputs and outputs
    vector<MX> arg(n_dep());
    for (casadi_int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (casadi_int i=0; i<res.size(); ++i) res[i] = get_output(i);

    // Number of derivatives
    casadi_int nadj = aseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nadj);
    std::vector<casadi_int> col_offset(nadj+1, 0);
    for (casadi_int d=0; d<nadj; ++d) {
      rhs[d] = aseed[d][0];
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_.solve(A, horzcat(rhs), !Tr), col_offset);

    // Collect sensitivities
    asens.resize(nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d].resize(2);

      // Propagate to A
      MX a;
      if (!Tr) {
        a = -mac(rhs[d], X.T(), MX::zeros(A.sparsity()));
      } else {
        a = -mac(X, rhs[d].T(), MX::zeros(A.sparsity()));
      }
      if (asens[d][1].is_empty(true)) {
        asens[d][1] = a;
      } else {
        asens[d][1] += a;
      }

      // Propagate to B
      if (asens[d][0].is_empty(true)) {
        asens[d][0] = rhs[d];
      } else {
        asens[d][0] += rhs[d];
      }
    }
  }

  template<bool Tr>
  int Solve<Tr>::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Number of right-hand-sides
    casadi_int nrhs = dep(0).size2();

    // Sparsities
    const Sparsity& A_sp = dep(1).sparsity();
    const casadi_int* A_colind = A_sp.colind();
    const casadi_int* A_row = A_sp.row();
    casadi_int n = A_sp.size1();

    // Get pointers to data
    const bvec_t *B=arg[0], *A = arg[1];
    bvec_t* X = res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (casadi_int r=0; r<nrhs; ++r) {
      // Copy B to a temporary vector
      copy(B, B+n, tmp);

      // Add A_hat contribution to tmp
      for (casadi_int cc=0; cc<n; ++cc) {
        for (casadi_int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          casadi_int rr = A_row[k];
          tmp[Tr ? cc : rr] |= A[k];
        }
      }

      // Propagate to X
      std::fill(X, X+n, 0);
      A_sp.spsolve(X, tmp, Tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
    return 0;
  }

  template<bool Tr>
  int Solve<Tr>::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Number of right-hand-sides
    casadi_int nrhs = dep(0).size2();

    // Sparsities
    const Sparsity& A_sp = dep(1).sparsity();
    const casadi_int* A_colind = A_sp.colind();
    const casadi_int* A_row = A_sp.row();
    casadi_int n = A_sp.size1();

    // Get pointers to data
    bvec_t *B=arg[0], *A=arg[1], *X=res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (casadi_int r=0; r<nrhs; ++r) {
      // Solve transposed
      std::fill(tmp, tmp+n, 0);
      A_sp.spsolve(tmp, X, !Tr);

      // Clear seeds
      std::fill(X, X+n, 0);

      // Propagate to B
      for (casadi_int i=0; i<n; ++i) B[i] |= tmp[i];

      // Propagate to A
      for (casadi_int cc=0; cc<n; ++cc) {
        for (casadi_int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          casadi_int rr = A_row[k];
          A[k] |= tmp[Tr ? cc : rr];
        }
      }

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
    return 0;
  }

  template<bool Tr>
  size_t Solve<Tr>::sz_w() const {
    return sparsity().size1();
  }

  template<bool Tr>
  void Solve<Tr>::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = dep(0).size2();

    // Array for x
    g.local("rr", "casadi_real", "*");
    g << "rr = " << g.work(res[0], nnz()) << ";\n";

    // Array for A
    g.local("ss", "casadi_real", "*");
    g << "ss = " << g.work(arg[1], dep(1).nnz()) << ";\n";

    // Copy b to x if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], nnz()), nnz(), "rr") << '\n';
    }
    // Solver specific codegen
    linsol_->generate(g, "ss", "rr", nrhs, Tr);
  }

} // namespace casadi

#endif // CASADI_SOLVE_IMPL_HPP
