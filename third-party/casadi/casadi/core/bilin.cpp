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


#include "bilin.hpp"

using namespace std;
namespace casadi {

  Bilin::Bilin(const MX& A, const MX& x, const MX& y) {
    casadi_assert(x.is_column(), "Dimension mismatch");
    casadi_assert(y.is_column(), "Dimension mismatch");
    set_dep(A, densify(x), densify(y));
    set_sparsity(Sparsity::scalar());
  }

  std::string Bilin::disp(const std::vector<std::string>& arg) const {
    return "bilin(" + arg.at(0) + ", " + arg.at(1) + ", " + arg.at(2) + ")";
  }

  void Bilin::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = bilin(arg[0], arg[1], arg[2]);
  }

  void Bilin::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0]
        = bilin(fseed[d][0], dep(1), dep(2))
        + bilin(dep(0), fseed[d][1], dep(2))
        + bilin(dep(0), dep(1), fseed[d][2]);
    }
  }

  void Bilin::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] = rank1(project(asens[d][0], dep(0).sparsity()),
                          aseed[d][0], dep(1), dep(2));
      asens[d][1] += aseed[d][0] * mtimes(dep(0), dep(2));
      asens[d][2] += aseed[d][0] * mtimes(dep(0).T(), dep(1));
    }
  }

  int Bilin::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Bilin::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Bilin::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_bilin(arg[0], dep(0).sparsity(), arg[1], arg[2]);
    return 0;
  }

  int Bilin::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /* Return value */
    bvec_t r=0;

    /* Loop over the columns of A */
    SparsityStruct sp_A = dep(0).sparsity();
    casadi_int cc, rr, el;
    for (cc=0; cc<sp_A.ncol; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=sp_A.colind[cc]; el<sp_A.colind[cc+1]; ++el) {
        /* Get the row */
        rr=sp_A.row[el];

        /* Add contribution */
        r |= arg[1][rr] | arg[0][el] | arg[2][cc];
      }
    }
    *res[0] = r;
    return 0;
  }

  int Bilin::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    /* Seed */
    bvec_t s_r=res[0][0]; res[0][0] = 0;

    /* Loop over the columns of A */
    SparsityStruct sp_A = dep(0).sparsity();
    casadi_int cc, rr, el;
    for (cc=0; cc<sp_A.ncol; ++cc) {
      /* Loop over the nonzeros of A */
      for (el=sp_A.colind[cc]; el<sp_A.colind[cc+1]; ++el) {
        /* Get the row */
        rr=sp_A.row[el];

        /* Add contribution */
        arg[0][el] |= s_r;
        arg[1][rr] |= s_r;
        arg[2][cc] |= s_r;
      }
    }
    return 0;
  }

  void Bilin::generate(CodeGenerator& g,
                       const std::vector<casadi_int>& arg,
                       const std::vector<casadi_int>& res) const {
    g << g.workel(res[0]) << " = "
      << g.bilin(g.work(arg[0], dep(0).nnz()), dep(0).sparsity(),
                 g.work(arg[1], dep(1).nnz()),
                 g.work(arg[2], dep(2).nnz())) << ";\n";
  }

} // namespace casadi
