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


#include "constant_mx.hpp"
#include <vector>
#include <algorithm>
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  ConstantMX::ConstantMX(const Sparsity& sp) {
    set_sparsity(sp);
  }

  ConstantMX::~ConstantMX() {
  }

  bool ConstantMX::is_valid_input() const {
    return nnz()==0;
  }

  casadi_int ConstantMX::n_primitives() const {
    if (nnz()==0) {
      return 0;
    } else {
      return MXNode::n_primitives();
    }
  }

  void ConstantMX::primitives(std::vector<MX>::iterator& it) const {
    if (nnz()!=0) {
      MXNode::primitives(it);
    }
  }

  void ConstantMX::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    if (nnz()!=0) {
      MXNode::split_primitives(x, it);
    }
  }

  MX ConstantMX::join_primitives(std::vector<MX>::const_iterator& it) const {
    if (nnz()==0) {
      return MX(sparsity());
    } else {
      return MXNode::join_primitives(it);
    }
  }

  void ConstantMX::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = shared_from_this<MX>();
  }

 void ConstantMX::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
   MX zero_sens(size1(), size2());
   for (casadi_int d=0; d<fsens.size(); ++d) {
     fsens[d][0] = zero_sens;
   }
 }

  void ConstantMX::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                           std::vector<std::vector<MX> >& asens) const {
  }

  int ConstantMX::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  int ConstantMX::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    fill_n(res[0], nnz(), 0);
    return 0;
  }

  void ConstantDM::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    // Print the constant
    string ind = g.constant(x_.nonzeros());

    // Copy the constant to the work vector
    g << g.copy(ind, nnz(), g.work(res[0], nnz())) << '\n';
  }

  bool ConstantMX::__nonzero__() const {
    if (numel()!=1) casadi_error("Can only determine truth value of scalar MX.");
    if (nnz()!=1) casadi_error("Can only determine truth value of dense scalar MX.");
    return !is_zero();
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, casadi_int val) {
    if (sp.is_empty(true)) {
      return ZeroByZero::getInstance();
    } else {
      switch (val) {
      case 0: return new Constant<CompiletimeConst<0> >(sp);
      case 1: return new Constant<CompiletimeConst<1> >(sp);
      case -1: return new Constant<CompiletimeConst<(-1)> >(sp);
      default: return new Constant<RuntimeConst<casadi_int> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Sparsity& sp, double val) {
    if (sp.is_empty(true)) {
      return ZeroByZero::getInstance();
    } else {
      casadi_int intval = static_cast<casadi_int>(val);
      if (static_cast<double>(intval)-val==0) {
        return create(sp, intval);
      } else {
        return new Constant<RuntimeConst<double> >(sp, val);
      }
    }
  }

  ConstantMX* ConstantMX::create(const Matrix<double>& val) {
    if (val.nnz()==0) {
      return create(val.sparsity(), 0);
    } else if (val.is_scalar()) {
      return create(val.sparsity(), val.scalar());
    } else {
      // Check if all values are the same
      const vector<double> vdata = val.nonzeros();
      double v = vdata[0];
      for (auto&& i : vdata) {
        if (i!=v) {
          // Values not all the same
          return new ConstantDM(val);
        }
      }

      // All values identical if reached this point
      return create(val.sparsity(), v);
    }
  }

  bool ConstantDM::is_zero() const {
    return x_.is_zero();
  }

  bool ConstantDM::is_one() const {
    return x_.is_one();
  }

  bool ConstantDM::is_minus_one() const {
    return x_.is_minus_one();
  }

  bool ConstantDM::is_eye() const {
    return x_.is_eye();
  }

  // MX ConstantMX::get_mac(const MX& y) const {
  //   if (y.is_constant()) {
  //     // Constant folding
  //     DM xv = get_DM();
  //     DM yv = y->get_DM();
  //     return mul(xv, yv);
  //   } else {
  //     return MXNode::get_mac(y);
  //   }
  // }

  MX ConstantMX::get_dot(const MX& y) const {
    if (y.is_constant()) {
      // Constant folding
      DM xv = get_DM();
      DM yv = y->get_DM();
      return dot(xv, yv);
    } else {
      return MXNode::get_dot(y);
    }
  }

  bool ConstantDM::is_equal(const MXNode* node, casadi_int depth) const {
    // Check if same node
    const ConstantDM* n = dynamic_cast<const ConstantDM*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (!std::equal(x_->begin(), x_->end(), n->x_->begin())) return false;

    return true;
  }

  std::string ZeroByZero::disp(const std::vector<std::string>& arg) const {
    return "0x0";
  }

  MX ZeroByZero::get_project(const Sparsity& sp) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    casadi_assert_dev(nz.empty());
    return MX::zeros(sp);
  }

  MX ZeroByZero::get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_transpose() const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_unary(casadi_int op) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::_get_binary(casadi_int op, const MX& y, bool ScX, bool ScY) const {
    return shared_from_this<MX>();
  }

  MX ZeroByZero::get_reshape(const Sparsity& sp) const {
    casadi_assert_dev(sp.is_empty());
    return MX::zeros(sp);
  }

} // namespace casadi
