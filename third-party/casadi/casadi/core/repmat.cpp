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


#include "repmat.hpp"
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  HorzRepmat::HorzRepmat(const MX& x, casadi_int n) : n_(n) {
    set_dep(x);
    set_sparsity(repmat(x.sparsity(), 1, n));
  }

  std::string HorzRepmat::disp(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    ss << "repmat("  << arg.at(0) << ", " << n_ << ")";
    return ss.str();
  }

  template<typename T>
  int HorzRepmat::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    casadi_int nnz = dep(0).nnz();
    for (casadi_int i=0; i<n_; ++i) {
      std::copy(arg[0], arg[0]+nnz, res[0]+i*nnz);
    }
    return 0;
  }

  int HorzRepmat::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int HorzRepmat::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  void HorzRepmat::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_repmat(1, n_);
  }

  static bvec_t Orring(bvec_t x, bvec_t y) { return x | y; }

  int HorzRepmat::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nnz = dep(0).nnz();
    std::fill(res[0], res[0]+nnz, 0);
    return eval_gen<bvec_t>(arg, res, iw, w);
  }

  int HorzRepmat::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nnz = dep(0).nnz();
    for (casadi_int i=0;i<n_;++i) {
      std::transform(res[0]+i*nnz, res[0]+(i+1)*nnz, arg[0], arg[0], &Orring);
    }
    std::fill(res[0], res[0]+nnz, 0);
    return 0;
  }

  void HorzRepmat::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = fseed[d][0]->get_repmat(1, n_);
    }
  }

  void HorzRepmat::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<asens.size(); ++d) {
      asens[d][0] += aseed[d][0]->get_repsum(1, n_);
    }
  }

  void HorzRepmat::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    casadi_int nnz = dep(0).nnz();
    g.local("i", "casadi_int");
    g << "for (i=0;i<" << n_ << ";++i) {\n"
      << "    " << g.copy(g.work(arg[0], dep(0).nnz()), nnz,
                          g.work(res[0], sparsity().nnz()) + "+ i*" + str(nnz)) << "\n"
      << "  }\n";
  }

  HorzRepsum::HorzRepsum(const MX& x, casadi_int n) : n_(n) {
    casadi_assert_dev(x.size2() % n == 0);
    std::vector<Sparsity> sp = horzsplit(x.sparsity(), x.size2()/n);
    Sparsity block = sp[0];
    for (casadi_int i=1;i<sp.size();++i) {
      block = block+sp[i];
    }
    Sparsity goal = repmat(block, 1, n);
    set_dep(project(x, goal));
    set_sparsity(block);
  }

  std::string HorzRepsum::disp(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    ss << "repsum("  << arg.at(0) << ", " << n_ << ")";
    return ss.str();
  }

  template<typename T, typename R>
  int HorzRepsum::eval_gen(const T** arg, T** res, casadi_int* iw, T* w,
                           R reduction) const {
    casadi_int nnz = sparsity().nnz();
    fill_n(res[0], nnz, 0);
    for (casadi_int i=0;i<n_;++i) {
      std::transform(arg[0]+i*nnz, arg[0]+(i+1)*nnz, res[0], res[0], reduction);
    }
    return 0;
  }

  int HorzRepsum::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w, std::plus<double>());
  }

  int HorzRepsum::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w, std::plus<SXElem>());
  }

  void HorzRepsum::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_repsum(1, n_);
  }

  int HorzRepsum::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nnz = sparsity().nnz();
    std::fill(res[0], res[0]+nnz, 0);
    return eval_gen<bvec_t>(arg, res, iw, w, &Orring);
  }

  int HorzRepsum::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nnz = sparsity().nnz();
    for (casadi_int i=0;i<n_;++i) {
      std::transform(res[0], res[0]+nnz, arg[0]+i*nnz, arg[0]+i*nnz, &Orring);
    }
    std::fill(res[0], res[0]+nnz, 0);
    return 0;
  }

  void HorzRepsum::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = fseed[d][0]->get_repsum(1, n_);
    }
  }

  void HorzRepsum::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<asens.size(); ++d) {
      asens[d][0] += aseed[d][0]->get_repmat(1, n_);
    }
  }

  void HorzRepsum::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    casadi_int nnz = sparsity().nnz();
    g.local("i", "casadi_int");
    g.local("j", "casadi_int");
    g << g.fill(g.work(res[0], nnz), nnz, "0") << "\n"
      << "  for (i=0;i<" << n_ << ";++i) {\n"
      << "    for (j=0;j<" << nnz << ";++j) {\n"
      << "      " << g.work(res[0], nnz)<< "[j] += "
      << g.work(arg[0], dep(0).nnz()) << "[j+i*" << nnz << "];\n"
      << "    }\n"
      << "  }\n";
  }

} // namespace casadi
