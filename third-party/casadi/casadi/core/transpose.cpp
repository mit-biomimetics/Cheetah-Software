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


#include "transpose.hpp"

using namespace std;

namespace casadi {

  Transpose::Transpose(const MX& x) {
    set_dep(x);
    set_sparsity(x.sparsity().T());
  }

  int Transpose::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int DenseTranspose::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Transpose::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  int DenseTranspose::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Transpose::eval_gen(const T* const* arg, T* const* res,
                          casadi_int* iw, T* w) const {
    // Get sparsity patterns
    //const vector<casadi_int>& x_colind = input[0]->colind();
    const casadi_int* x_row = dep(0).row();
    casadi_int x_sz = dep(0).nnz();
    const casadi_int* xT_colind = sparsity().colind();
    casadi_int xT_ncol = sparsity().size2();

    const T* x = arg[0];
    T* xT = res[0];

    // Transpose
    copy(xT_colind, xT_colind+xT_ncol+1, iw);
    for (casadi_int el=0; el<x_sz; ++el) {
      xT[iw[x_row[el]]++] = x[el];
    }
    return 0;
  }

  template<typename T>
  int DenseTranspose::eval_gen(const T* const* arg, T* const* res,
                               casadi_int* iw, T* w) const {
    // Get sparsity patterns
    casadi_int x_nrow = dep().size1();
    casadi_int x_ncol = dep().size2();

    const T* x = arg[0];
    T* xT = res[0];
    for (casadi_int i=0; i<x_ncol; ++i) {
      for (casadi_int j=0; j<x_nrow; ++j) {
        xT[i+j*x_ncol] = x[j+i*x_nrow];
      }
    }
    return 0;
  }

  int Transpose::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Shortands
    const bvec_t *x = arg[0];
    bvec_t *xT = res[0];

    // Get sparsity
    casadi_int nz = nnz();
    const casadi_int* x_row = dep().row();
    const casadi_int* xT_colind = sparsity().colind();
    casadi_int xT_ncol = sparsity().size2();

    // Loop over the nonzeros of the argument
    copy(xT_colind, xT_colind+xT_ncol+1, iw);
    for (casadi_int el=0; el<nz; ++el) {
      xT[iw[*x_row++]++] = *x++;
    }
    return 0;
  }

  int Transpose::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Shortands
    bvec_t *x = arg[0];
    bvec_t *xT = res[0];

    // Get sparsity
    casadi_int nz = nnz();
    const casadi_int* x_row = dep().row();
    const casadi_int* xT_colind = sparsity().colind();
    casadi_int xT_ncol = sparsity().size2();

    // Loop over the nonzeros of the argument
    copy(xT_colind, xT_colind+xT_ncol+1, iw);
    for (casadi_int el=0; el<nz; ++el) {
      casadi_int elT = iw[*x_row++]++;
      *x++ |= xT[elT];
      xT[elT] = 0;
    }
    return 0;
  }

  int DenseTranspose::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Shorthands
    const bvec_t *x = arg[0];
    bvec_t *xT = res[0];
    casadi_int x_nrow = dep().size1();
    casadi_int x_ncol = dep().size2();

    // Loop over the elements
    for (casadi_int rr=0; rr<x_nrow; ++rr) {
      for (casadi_int cc=0; cc<x_ncol; ++cc) {
        *xT++ = x[rr+cc*x_nrow];
      }
    }
    return 0;
  }

  int DenseTranspose::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Shorthands
    bvec_t *x = arg[0];
    bvec_t *xT = res[0];
    casadi_int x_nrow = dep().size1();
    casadi_int x_ncol = dep().size2();

    // Loop over the elements
    for (casadi_int rr=0; rr<x_nrow; ++rr) {
      for (casadi_int cc=0; cc<x_ncol; ++cc) {
        x[rr+cc*x_nrow] |= *xT;
        *xT++ = 0;
      }
    }
    return 0;
  }

  std::string Transpose::disp(const std::vector<std::string>& arg) const {
    return arg.at(0) + "'";
  }

  void Transpose::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0].T();
  }

  void Transpose::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = fseed[d][0].T();
    }
  }

  void Transpose::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0].T();
    }
  }

  void Transpose::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    g << g.trans(g.work(arg[0], nnz()), dep().sparsity(),
                 g.work(res[0], nnz()), sparsity(), "iw") <<  ";\n";
  }

  void DenseTranspose::generate(CodeGenerator& g,
                                const std::vector<casadi_int>& arg,
                                const std::vector<casadi_int>& res) const {
    g.local("cs", "const casadi_real", "*");
    g.local("rr", "casadi_real", "*");
    g.local("i", "casadi_int");
    g.local("j", "casadi_int");
    g << "for (i=0, rr=" << g.work(res[0], nnz()) << ", "
      << "cs=" << g.work(arg[0], nnz()) << "; i<" << dep().size2() << "; ++i) "
      << "for (j=0; j<" << dep().size1() << "; ++j) "
      << "rr[i+j*" << dep().size2() << "] = *cs++;\n";
  }

} // namespace casadi
