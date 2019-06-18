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


#include "symbolic_mx.hpp"
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  SymbolicMX::SymbolicMX(const std::string& name, casadi_int nrow, casadi_int ncol) : name_(name) {
    set_sparsity(Sparsity::dense(nrow, ncol));
  }

  SymbolicMX::SymbolicMX(const std::string& name, const Sparsity & sp) : name_(name) {
    set_sparsity(sp);
  }

  std::string SymbolicMX::disp(const std::vector<std::string>& arg) const {
    return name_;
  }

  int SymbolicMX::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return 0;
  }

  int SymbolicMX::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return 0;
  }

  void SymbolicMX::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
  }

  void SymbolicMX::ad_forward(const std::vector<std::vector<MX> >& fseed,
                           std::vector<std::vector<MX> >& fsens) const {
  }

  void SymbolicMX::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                           std::vector<std::vector<MX> >& asens) const {
  }

  const std::string& SymbolicMX::name() const {
    return name_;
  }

  int SymbolicMX::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  int SymbolicMX::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  bool SymbolicMX::has_duplicates() const {
    if (this->temp!=0) {
      casadi_warning("Duplicate expression: " + name());
      return true;
    } else {
      this->temp = 1;
      return false;
    }
  }

  void SymbolicMX::reset_input() const {
    this->temp = 0;
  }

} // namespace casadi
