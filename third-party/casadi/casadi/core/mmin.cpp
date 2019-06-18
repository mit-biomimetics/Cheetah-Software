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


#include "mmin.hpp"

using namespace std;
namespace casadi {

  MMin::MMin(const MX& x) {
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  std::string MMin::disp(const std::vector<std::string>& arg) const {
    return "min(" + arg.at(0) + ")";
  }

  MMax::MMax(const MX& x) {
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  std::string MMax::disp(const std::vector<std::string>& arg) const {
    return "max(" + arg.at(0) + ")";
  }

  int MMin::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (!res[0]) return 0;
    res[0][0] = casadi_mmin(arg[0], dep(0).nnz(), dep(0).is_dense());
    return 0;
  }

  int MMax::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (!res[0]) return 0;
    res[0][0] = casadi_mmax(arg[0], dep(0).nnz(), dep(0).is_dense());

    return 0;
  }

  int MMin::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (!res[0]) return 0;
    res[0][0] = casadi_mmin(arg[0], dep(0).nnz(), dep(0).is_dense());

    return 0;
  }

  int MMax::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (!res[0]) return 0;
    res[0][0] = casadi_mmax(arg[0], dep(0).nnz(), dep(0).is_dense());

    return 0;
  }

  void MMin::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = mmin(arg[0]);
  }

  void MMax::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = mmax(arg[0]);
  }

} // namespace casadi
