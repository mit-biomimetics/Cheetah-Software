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


#include "norm.hpp"

using namespace std;
namespace casadi {

  Norm::Norm(const MX& x) {
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  std::string NormF::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_F";
  }

  int NormF::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int NormF::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    eval_gen<SXElem>(arg, res, iw, w);
    return 0;
  }

  template<typename T>
  int NormF::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    *res[0] = casadi_norm_2(dep().nnz(), arg[0]);
    return 0;
  }

  void NormF::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_norm_fro();
  }

  void NormF::ad_forward(const std::vector<std::vector<MX> >& fseed,
                      std::vector<std::vector<MX> >& fsens) const {
    MX self = shared_from_this<MX>();
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = dep(0)->get_dot(fseed[d][0]) / self;
    }
  }

  void NormF::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                      std::vector<std::vector<MX> >& asens) const {
    MX self = shared_from_this<MX>();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += (aseed[d][0]/self) * dep(0);
    }
  }

  void NormF::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    g << g.workel(res[0]) << " = sqrt("
      << g.dot(dep().nnz(), g.work(arg[0], dep(0).nnz()), g.work(arg[0], dep(0).nnz()))
      << ");\n";
  }

  std::string Norm2::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_2";
  }

  std::string Norm1::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_1";
  }

  std::string NormInf::disp(const std::vector<std::string>& arg) const {
    return "||" + arg.at(0) + "||_inf";
  }

} // namespace casadi
