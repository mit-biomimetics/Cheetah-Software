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


#include "dot.hpp"
using namespace std;
namespace casadi {

  Dot::Dot(const MX& x, const MX& y) {
    casadi_assert_dev(x.sparsity()==y.sparsity());
    set_dep(x, y);
    set_sparsity(Sparsity::scalar());
  }

  std::string Dot::disp(const std::vector<std::string>& arg) const {
    return "dot(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void Dot::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_dot(arg[1]);
  }

  void Dot::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->get_dot(fseed[d][1])
        + fseed[d][0]->get_dot(dep(1));
    }
  }

  void Dot::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0] * dep(1);
      asens[d][1] += aseed[d][0] * dep(0);
    }
  }

  int Dot::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Dot::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Dot::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_dot(dep(0).nnz(), arg[0], arg[1]);
    return 0;
  }

  int Dot::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a0=arg[0], *a1=arg[1];
    bvec_t* r = res[0];
    const casadi_int n = dep(0).nnz();
    *r = 0;
    for (casadi_int i=0; i<n; ++i) {
      *r |= *a0++ | *a1++;
    }
    return 0;
  }

  int Dot::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a0=arg[0], *a1=arg[1], *r=res[0];
    const casadi_int n = dep(0).nnz();
    for (casadi_int i=0; i<n; ++i) {
      *a0++ |= *r;
      *a1++ |= *r;
    }
    *r = 0;
    return 0;
  }

  void Dot::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    g << g.workel(res[0]) << " = "
      << g.dot(dep().nnz(), g.work(arg[0], dep(0).nnz()), g.work(arg[1], dep(1).nnz()))
      << ";\n";
  }

} // namespace casadi
