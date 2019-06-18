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


#include "sx_elem.hpp"
#include "matrix.hpp"
#include <stack>
#include <cassert>
#include "calculus.hpp"
#include "constant_sx.hpp"
#include "symbolic_sx.hpp"
#include "unary_sx.hpp"
#include "binary_sx.hpp"
#include "global_options.hpp"
#include "sx_function.hpp"

using namespace std;
namespace casadi {



  // Allocate storage for the caching
  CACHING_MAP<casadi_int, IntegerSX*> IntegerSX::cached_constants_;
  CACHING_MAP<double, RealtypeSX*> RealtypeSX::cached_constants_;

  SXElem::SXElem() {
    node = casadi_limits<SXElem>::nan.node;
    node->count++;
  }

  SXElem::SXElem(SXNode* node_, bool dummy) : node(node_) {
    node->count++;
  }

  SXElem SXElem::create(SXNode* node) {
    return SXElem(node, false);
  }

  SXElem::SXElem(const SXElem& scalar) {
    node = scalar.node;
    node->count++;
  }

  SXElem::SXElem(double val) {
    // Only ints fit here, not casadi_int
    int intval = static_cast<int>(val);
    if (val-static_cast<double>(intval) == 0) { // check if integer
      if (intval == 0)             node = casadi_limits<SXElem>::zero.node;
      else if (intval == 1)        node = casadi_limits<SXElem>::one.node;
      else if (intval == 2)        node = casadi_limits<SXElem>::two.node;
      else if (intval == -1)       node = casadi_limits<SXElem>::minus_one.node;
      else                        node = IntegerSX::create(intval);
      node->count++;
    } else {
      if (isnan(val))              node = casadi_limits<SXElem>::nan.node;
      else if (isinf(val))         node = val > 0 ? casadi_limits<SXElem>::inf.node :
                                      casadi_limits<SXElem>::minus_inf.node;
      else                        node = RealtypeSX::create(val);
      node->count++;
    }
  }

  SXElem SXElem::sym(const std::string& name) {
    return create(new SymbolicSX(name));
  }

  SXElem::~SXElem() {
    if (--node->count == 0) delete node;
  }

  SXElem& SXElem::operator=(const SXElem &scalar) {
    // quick return if the old and new pointers point to the same object
    if (node == scalar.node) return *this;

    // decrease the counter and delete if this was the last pointer
    if (--node->count == 0) delete node;

    // save the new pointer
    node = scalar.node;
    node->count++;
    return *this;
  }

  void SXElem::assignIfDuplicate(const SXElem& scalar, casadi_int depth) {
    casadi_assert_dev(depth>=1);
    if (!is_equal(*this, scalar, 0) && is_equal(*this, scalar, depth)) {
      *this = scalar;
    }
  }

  SXNode* SXElem::assignNoDelete(const SXElem& scalar) {
    // Return value
    SXNode* ret = node;

    // quick return if the old and new pointers point to the same object
    if (node == scalar.node) return ret;

    // decrease the counter but do not delete if this was the last pointer
    --node->count;

    // save the new pointer
    node = scalar.node;
    node->count++;

    // Return a pointer to the old node
    return ret;
  }

  SXElem& SXElem::operator=(double scalar) {
    return *this = SXElem(scalar);
  }

  void SXElem::disp(std::ostream& stream, bool more) const {
    node->disp(stream, more);
  }

  SXElem SXElem::operator-() const {
    if (is_op(OP_NEG))
      return dep();
    else if (is_zero())
      return 0;
    else if (is_minus_one())
      return 1;
    else if (is_one())
      return -1;
    else
      return UnarySX::create(OP_NEG, *this);
  }

  bool SXElem::__nonzero__() const {
    if (is_constant()) return !is_zero();
    casadi_error("Cannot compute the truth value of a CasADi SXElem symbolic expression.");
  }

  bool SXElem::is_doubled() const {
    return is_op(OP_ADD) && is_equal(dep(0), dep(1), SXNode::eq_depth_);
  }

  SXElem SXElem::inv() const {
    if (is_op(OP_INV)) {
      return dep(0);
    } else {
      return UnarySX::create(OP_INV, *this);
    }
  }

  SXNode* SXElem::get() const {
    return node;
  }

  const SXNode* SXElem::operator->() const {
    return node;
  }

  SXNode* SXElem::operator->() {
    return node;
  }

  SXElem SXElem::binary(casadi_int op, const SXElem& x, const SXElem& y) {
    // Simplifications
    if (GlobalOptions::simplification_on_the_fly) {
      switch (op) {
      case OP_ADD:
        if (x.is_zero())
          return y;
        else if (y->is_zero()) // term2 is zero
          return x;
        else if (y.is_op(OP_NEG))  // x + (-y) -> x - y
          return x - (-y);
        else if (x.is_op(OP_NEG)) // (-x) + y -> y - x
          return y - x.dep();
        else if (x.is_op(OP_MUL) && y.is_op(OP_MUL) &&
                 x.dep(0).is_constant() && static_cast<double>(x.dep(0))==0.5 &&
                 y.dep(0).is_constant() && static_cast<double>(y.dep(0))==0.5 &&
                 is_equal(y.dep(1), x.dep(1), SXNode::eq_depth_)) // 0.5x+0.5x = x
          return x.dep(1);
        else if (x.is_op(OP_DIV) && y.is_op(OP_DIV) &&
                 x.dep(1).is_constant() && static_cast<double>(x.dep(1))==2 &&
                 y.dep(1).is_constant() && static_cast<double>(y.dep(1))==2 &&
                 is_equal(y.dep(0), x.dep(0), SXNode::eq_depth_)) // x/2+x/2 = x
          return x.dep(0);
        else if (x.is_op(OP_SUB) && is_equal(x.dep(1), y, SXNode::eq_depth_))
          return x.dep(0);
        else if (y.is_op(OP_SUB) && is_equal(x, y.dep(1), SXNode::eq_depth_))
          return y.dep(0);
        else if (x.is_op(OP_SQ) && y.is_op(OP_SQ) &&
                 ((x.dep().is_op(OP_SIN) && y.dep().is_op(OP_COS))
                  || (x.dep().is_op(OP_COS) && y.dep().is_op(OP_SIN)))
                 && is_equal(x.dep().dep(), y.dep().dep(), SXNode::eq_depth_))
          return 1; // sin^2 + cos^2 -> 1
        break;
      case OP_SUB:
        if (y->is_zero()) // term2 is zero
          return x;
        if (x.is_zero()) // term1 is zero
          return -y;
        if (is_equal(x, y, SXNode::eq_depth_)) // the terms are equal
          return 0;
        else if (y.is_op(OP_NEG)) // x - (-y) -> x + y
          return x + y.dep();
        else if (x.is_op(OP_ADD) && is_equal(x.dep(1), y, SXNode::eq_depth_))
          return x.dep(0);
        else if (x.is_op(OP_ADD) && is_equal(x.dep(0), y, SXNode::eq_depth_))
          return x.dep(1);
        else if (y.is_op(OP_ADD) && is_equal(x, y.dep(1), SXNode::eq_depth_))
          return -y.dep(0);
        else if (y.is_op(OP_ADD) && is_equal(x, y.dep(0), SXNode::eq_depth_))
          return -y.dep(1);
        else if (x.is_op(OP_NEG))
          return -(x.dep() + y);
        break;
      case OP_MUL:
        if (is_equal(y, x, SXNode::eq_depth_))
          return sq(x);
        else if (!x.is_constant() && y.is_constant())
          return y * x;
        // Make sure NaN does not propagate through an inactive branch
        // On demand by Deltares, July 2016
        else if (x.is_op(OP_IF_ELSE_ZERO))
          return if_else_zero(x.dep(0), x.dep(1)*y);
        else if (y.is_op(OP_IF_ELSE_ZERO))
          return if_else_zero(y.dep(0), y.dep(1)*x);
        else if (x.is_zero() || y->is_zero()) // one of the terms is zero
          return 0;
        else if (x.is_one()) // term1 is one
          return y;
        else if (y->is_one()) // term2 is one
          return x;
        else if (y->is_minus_one())
          return -x;
        else if (x.is_minus_one())
          return -y;
        else if (y.is_op(OP_INV))
          return x/y.inv();
        else if (x.is_op(OP_INV))
          return y / x.inv();
        else if (x.is_constant() && y.is_op(OP_MUL) && y.dep(0).is_constant() &&
                 static_cast<double>(x)*static_cast<double>(y.dep(0))==1) // 5*(0.2*x) = x
          return y.dep(1);
        else if (x.is_constant() && y.is_op(OP_DIV) && y.dep(1).is_constant() &&
                 static_cast<double>(x)==static_cast<double>(y.dep(1))) // 5*(x/5) = x
          return y.dep(0);
        else if (x.is_op(OP_DIV) && is_equal(x.dep(1), y, SXNode::eq_depth_)) // ((2/x)*x)
          return x.dep(0);
        else if (y.is_op(OP_DIV) &&
                 is_equal(y.dep(1), x, SXNode::eq_depth_)) // ((2/x)*x)
          return y.dep(0);
        else if (x.is_op(OP_NEG))
          return -(x.dep() * y);
        else if (y.is_op(OP_NEG))
          return -(x * y.dep());
        break;
      case OP_DIV:
        if (y->is_zero()) // term2 is zero
          return casadi_limits<SXElem>::nan;
        else if (x.is_zero()) // term1 is zero
          return 0;
        else if (y->is_one()) // term2 is one
          return x;
        else if (y->is_minus_one())
          return -x;
        else if (is_equal(x, y, SXNode::eq_depth_)) // terms are equal
          return 1;
        else if (x.is_doubled() && is_equal(y, 2))
          return x.dep(0);
        else if (x.is_op(OP_MUL) && is_equal(y, x.dep(0), SXNode::eq_depth_))
          return x.dep(1);
        else if (x.is_op(OP_MUL) && is_equal(y, x.dep(1), SXNode::eq_depth_))
          return x.dep(0);
        else if (x.is_one())
          return y.inv();
        else if (y.is_op(OP_INV))
          return x*y.inv();
        else if (x.is_doubled() && y.is_doubled())
          return x.dep(0) / y->dep(0);
        else if (y.is_constant() && x.is_op(OP_DIV) && x.dep(1).is_constant() &&
                 static_cast<double>(y)*static_cast<double>(x.dep(1))==1) // (x/5)/0.2
          return x.dep(0);
        else if (y.is_op(OP_MUL) &&
                 is_equal(y.dep(1), x, SXNode::eq_depth_)) // x/(2*x) = 1/2
          return BinarySX::create(OP_DIV, 1, y.dep(0));
        else if (x.is_op(OP_NEG) &&
                 is_equal(x.dep(0), y, SXNode::eq_depth_))      // (-x)/x = -1
          return -1;
        else if (y.is_op(OP_NEG) &&
                 is_equal(y.dep(0), x, SXNode::eq_depth_))      // x/(-x) = 1
          return -1;
        else if (y.is_op(OP_NEG) && x.is_op(OP_NEG) &&
                 is_equal(x.dep(0), y.dep(0), SXNode::eq_depth_))  // (-x)/(-x) = 1
          return 1;
        else if (x.is_op(OP_DIV) && is_equal(y, x.dep(0), SXNode::eq_depth_))
          return x.dep(1).inv();
        else if (x.is_op(OP_NEG))
          return -(x.dep() / y);
        else if (y.is_op(OP_NEG))
          return -(x / y.dep());
        break;
      case OP_POW:
        if (y->is_constant()) {
          if (y->is_integer()) {
            casadi_int nn = y->to_int();
            if (nn == 0) {
              return 1;
            } else if (nn>100 || nn<-100) { // maximum depth
              return binary(OP_CONSTPOW, x, nn);
            } else if (nn<0) { // negative power
              return 1/pow(x, -nn);
            } else if (nn%2 == 1) { // odd power
              return x*pow(x, nn-1);
            } else { // even power
              SXElem rt = pow(x, nn/2);
              return rt*rt;
            }
          } else if (y->to_double()==0.5) {
            return sqrt(x);
          } else {
            return binary(OP_CONSTPOW, x, y);
          }
        }
        break;
        case OP_LE:
        if ((y-x).is_nonnegative())
          return 1;
        break;
      case OP_LT:
        if (((x)-y).is_nonnegative())
          return 0;
        break;
      case OP_EQ:
        if (is_equal(x, y))
          return 1;
        break;
      case OP_NE:
        if (is_equal(x, y))
          return 0;
        break;
      case OP_IF_ELSE_ZERO:
        if (y->is_zero()) {
          return y;
        } else if (x.is_constant()) {
          if (static_cast<double>(x)!=0) {
            return y;
          } else {
            return 0;
          }
        }
      }
    }
    return BinarySX::create(Operation(op), x, y);
  }

  SXElem SXElem::unary(casadi_int op, const SXElem& x) {
    // Simplifications
    if (GlobalOptions::simplification_on_the_fly) {
      switch (op) {
        case OP_SQ:
          if (x.is_op(OP_SQRT))
            return x.dep();
          else if (x.is_op(OP_NEG))
            return sq(x.dep());
          break;
        case OP_FABS:
          if (x.is_op(OP_FABS) || x.is_op(OP_SQ))
            return x;
          break;
        case OP_NOT:
          if (x.is_op(OP_NOT))
            return x.dep();
          break;
        case OP_SINH:
        case OP_TANH:
        case OP_ATANH:
        case OP_ACOSH:
        case OP_ASINH:
          if (x.is_zero())
            return 0;
          break;
        case OP_COSH:
          if (x.is_zero())
            return 1;
          break;
        case OP_SQRT:
          if (x.is_op(OP_SQ))
            return fabs(x.dep());
          break;
      }
    }
    return UnarySX::create(Operation(op), x);
  }

  bool SXElem::is_leaf() const {
    if (!node) return true;
    return is_constant() || is_symbolic();
  }

  bool SXElem::is_commutative() const {
    casadi_assert(n_dep(), "SX::is_commutative: must be binary");
    return operation_checker<CommChecker>(op());
  }

  bool SXElem::is_constant() const {
    return node->is_constant();
  }

  bool SXElem::is_integer() const {
    return node->is_integer();
  }

  bool SXElem::is_symbolic() const {
    return node->is_symbolic();
  }

  bool SXElem::is_zero() const {
    return node->is_zero();
  }

  bool SXElem::is_almost_zero(double tol) const {
    return node->is_almost_zero(tol);
  }

  bool SXElem::is_one() const {
    return node->is_one();
  }

  bool SXElem::is_minus_one() const {
    return node->is_minus_one();
  }

  bool SXElem::is_nan() const {
    return node->is_nan();
  }

  bool SXElem::is_inf() const {
    return node->is_inf();
  }

  bool SXElem::is_minus_inf() const {
    return node->is_minus_inf();
  }

  const std::string& SXElem::name() const {
    return node->name();
  }

  casadi_int SXElem::op() const {
    return node->op();
  }

  bool SXElem::is_op(casadi_int op) const {
    return node->is_op(op);
  }

  bool SXElem::is_equal(const SXElem& x, const SXElem& y, casadi_int depth) {
    SXNode *x_node = x.get(), *y_node = y.get();
    if (x_node==y_node) {
      return true;
    } else if (depth>0) {
      return x_node->is_equal(y_node, depth);
    } else {
      return false;
    }
  }

  bool SXElem::is_nonnegative() const {
    if (is_constant()) {
      return static_cast<double>(*this)>=0;
    } else if (is_op(OP_SQ) || is_op(OP_FABS)) {
      return true;
    } else {
      return false;
    }
  }

  SXElem::operator double() const {
    return node->to_double();
  }

  SXElem::operator casadi_int() const {
    return node->to_int();
  }

  SXElem SXElem::dep(casadi_int ch) const {
    casadi_assert_dev(ch==0 || ch==1);
    return node->dep(ch);
  }

  casadi_int SXElem::n_dep() const {
    return node->n_dep();
  }

  casadi_int SXElem::__hash__() const {
    return reinterpret_cast<casadi_int>(node);
  }

  // node corresponding to a constant 0
  const SXElem casadi_limits<SXElem>::zero(new ZeroSX(), false);
  // node corresponding to a constant 1
  const SXElem casadi_limits<SXElem>::one(new OneSX(), false);
  // node corresponding to a constant 2
  const SXElem casadi_limits<SXElem>::two(IntegerSX::create(2), false);
  // node corresponding to a constant -1
  const SXElem casadi_limits<SXElem>::minus_one(new MinusOneSX(), false);
  const SXElem casadi_limits<SXElem>::nan(new NanSX(), false);
  const SXElem casadi_limits<SXElem>::inf(new InfSX(), false);
  const SXElem casadi_limits<SXElem>::minus_inf(new MinusInfSX(), false);

  bool casadi_limits<SXElem>::is_zero(const SXElem& val) {
    return val.is_zero();
  }

  bool casadi_limits<SXElem>::is_equal(const SXElem& x, const SXElem& y, casadi_int depth) {
    return SXElem::is_equal(x, y, depth);
  }

  bool casadi_limits<SXElem>::is_almost_zero(const SXElem& val, double tol) {
    return val.is_almost_zero(tol);
  }

  bool casadi_limits<SXElem>::is_one(const SXElem& val) {
    return val.is_one();
  }

  bool casadi_limits<SXElem>::is_minus_one(const SXElem& val) {
    return val.is_minus_one();
  }

  bool casadi_limits<SXElem>::is_constant(const SXElem& val) {
    return val.is_constant();
  }

  bool casadi_limits<SXElem>::is_integer(const SXElem& val) {
    return val.is_integer();
  }

  bool casadi_limits<SXElem>::is_inf(const SXElem& val) {
    return val.is_inf();
  }

  bool casadi_limits<SXElem>::is_minus_inf(const SXElem& val) {
    return val.is_minus_inf();
  }

  bool casadi_limits<SXElem>::is_nan(const SXElem& val) {
    return val.is_nan();
  }

  SXElem::operator SX() const {
    return SX(Sparsity::scalar(), *this, false);
  }

  int SXElem::get_temp() const {
    return (*this)->temp;
  }

  void SXElem::set_temp(int t) const {
    (*this)->temp = t;
  }

  bool SXElem::marked() const {
    return (*this)->marked();
  }

  void SXElem::mark() const {
    (*this)->mark();
  }

  bool SXElem::is_regular() const {
    if (is_constant()) {
      return !(is_nan() || is_inf() || is_minus_inf());
    } else {
      casadi_error("Cannot check regularity for symbolic SXElem");
    }
  }

} // namespace casadi

using namespace casadi;
namespace std {
  SXElem numeric_limits<SXElem>::infinity() throw() {
    return casadi::casadi_limits<SXElem>::inf;
  }

  SXElem numeric_limits<SXElem>::quiet_NaN() throw() {
    return casadi::casadi_limits<SXElem>::nan;
  }

  SXElem numeric_limits<SXElem>::min() throw() {
    return SXElem(numeric_limits<double>::min());
  }

  SXElem numeric_limits<SXElem>::max() throw() {
    return SXElem(numeric_limits<double>::max());
  }

  SXElem numeric_limits<SXElem>::epsilon() throw() {
    return SXElem(numeric_limits<double>::epsilon());
  }

  SXElem numeric_limits<SXElem>::round_error() throw() {
    return SXElem(numeric_limits<double>::round_error());
  }

} // namespace std
