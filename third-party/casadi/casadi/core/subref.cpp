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


#include "subref.hpp"

using namespace std;

namespace casadi {

  SubRef::SubRef(const MX& x, const Slice& i, const Slice& j) : i_(i), j_(j) {
    set_dep(x);
  }

  int SubRef::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int SubRef::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int SubRef::eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const {
    casadi_error("not ready");
    return 1;
  }

  int SubRef::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_error("not ready");
    return 1;
  }

  int SubRef::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_error("not ready");
    return 1;
  }

  std::string SubRef::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << i_ << ", " << j_ << "]";
    return ss.str();
  }

  void SubRef::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    casadi_error("not ready");
  }

  void SubRef::ad_forward(const std::vector<std::vector<MX> >& fseed,
                       std::vector<std::vector<MX> >& fsens) const {
    casadi_error("not ready");
  }

  void SubRef::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                       std::vector<std::vector<MX> >& asens) const {
    casadi_error("not ready");
  }

  void SubRef::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    casadi_error("not ready");
  }

  Dict SubRef::info() const {
    return {{"i", i_.info()}, {"j", j_.info()}};
  }

} // namespace casadi
