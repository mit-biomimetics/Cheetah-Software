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


#include "rank1.hpp"
using namespace std;
namespace casadi {

  Rank1::Rank1(const MX& A, const MX& alpha, const MX& x, const MX& y) {
    set_dep({A, alpha, x, y});
    set_sparsity(A.sparsity());
  }

  std::string Rank1::disp(const std::vector<std::string>& arg) const {
    return "rank1(" + arg.at(0) + ", " + arg.at(1)
      + ", " + arg.at(2) + ", " + arg.at(3) + ")";
  }

  void Rank1::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = rank1(arg[0], arg[1], arg[2], arg[3]);
  }

  void Rank1::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      MX v = project(fseed[d][0], sparsity());
      v = rank1(v, fseed[d][1], dep(2), dep(3));
      v = rank1(v, dep(1), fseed[d][2], dep(3));
      v = rank1(v, dep(1), dep(2), fseed[d][3]);
      fsens[d][0] = v;
    }
  }

  void Rank1::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][1] += bilin(aseed[d][0], dep(2), dep(3));
      asens[d][2] += dep(1) * mtimes(aseed[d][0], dep(3));
      asens[d][3] += dep(1) * mtimes(aseed[d][0].T(), dep(2));
      asens[d][0] += aseed[d][0];
    }
  }

  int Rank1::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Rank1::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Rank1::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    if (arg[0]!=res[0]) casadi_copy(arg[0], dep(0).nnz(), res[0]);
    casadi_rank1(res[0], sparsity(), *arg[1], arg[2], arg[3]);
    return 0;
  }

  int Rank1::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /* If not inline, copy to result */
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+dep(0).nnz(), res[0]);

    /* Get sparsities */
    casadi_int ncol_A = sparsity().size2();
    const casadi_int *colind_A = sparsity().colind(), *row_A = sparsity().row();

    /* Loop over the columns of A */
    casadi_int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get row */
        rr = row_A[el];

        /* Add the multiple */
        res[0][el] |= *arg[1] | arg[2][rr] | arg[3][cc];
      }
    }
    return 0;
  }

  int Rank1::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /* Get sparsities */
    casadi_int ncol_A = sparsity().size2();
    const casadi_int *colind_A = sparsity().colind(), *row_A = sparsity().row();

    /* Loop over the columns of A */
    casadi_int cc, rr, el;
    for (cc=0; cc<ncol_A; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
        /* Get row */
        rr = row_A[el];

        /* Add the multiple */
        *arg[1] |= res[0][el];
        arg[2][rr] |= res[0][el];
        arg[3][cc] |= res[0][el];
      }
    }

    // Clear seeds
    copy_rev(arg[0], res[0], nnz());
    return 0;
  }

  void Rank1::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << "\n";
    }

    // Perform operation inplace
    g << g.rank1(g.work(arg[0], dep(0).nnz()), sparsity(), g.workel(arg[1]),
                         g.work(arg[2], dep(2).nnz()), g.work(arg[3], dep(3).nnz())) << "\n";
  }

} // namespace casadi
