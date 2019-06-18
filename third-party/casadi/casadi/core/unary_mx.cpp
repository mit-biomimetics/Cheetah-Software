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


#include "unary_mx.hpp"
#include <vector>
#include <sstream>
#include "casadi_misc.hpp"
#include "global_options.hpp"

using namespace std;

namespace casadi {

  UnaryMX::UnaryMX(Operation op, MX x) : op_(op) {
    // Put a densifying node in between if necessary
    if (!operation_checker<F00Checker>(op_)) {
      x = densify(x);
    }

    set_dep(x);
    set_sparsity(x->sparsity());
  }

  std::string UnaryMX::disp(const std::vector<std::string>& arg) const {
    return casadi_math<double>::print(op_, arg.at(0));
  }

  int UnaryMX::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    double dummy = numeric_limits<double>::quiet_NaN();
    casadi_math<double>::fun(op_, arg[0], dummy, res[0], nnz());
    return 0;
  }

  int UnaryMX::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    SXElem dummy = 0;
    casadi_math<SXElem>::fun(op_, arg[0], dummy, res[0], nnz());
    return 0;
  }

  void UnaryMX::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    MX dummy;
    casadi_math<MX>::fun(op_, arg[0], dummy, res[0]);
  }

  void UnaryMX::ad_forward(const std::vector<std::vector<MX> >& fseed,
                     std::vector<std::vector<MX> >& fsens) const {
    // Get partial derivatives
    MX pd[2];
    MX dummy; // Function value, dummy second argument
    casadi_math<MX>::der(op_, dep(), dummy, shared_from_this<MX>(), pd);

    // Propagate forward seeds
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = pd[0]*fseed[d][0];
    }
  }

  void UnaryMX::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                     std::vector<std::vector<MX> >& asens) const {
    // Get partial derivatives
    MX pd[2];
    MX dummy; // Function value, dummy second argument
    casadi_math<MX>::der(op_, dep(), dummy, shared_from_this<MX>(), pd);

    // Propagate adjoint seeds
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += pd[0]*aseed[d][0];
    }
  }

  int UnaryMX::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_fwd(arg[0], res[0], nnz());
    return 0;
  }

  int UnaryMX::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_rev(arg[0], res[0], nnz());
    return 0;
  }

  void UnaryMX::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const {
    string r, x;
    if (nnz()==1) {
      // Scalar assignment
      r = g.workel(res[0]);
      x = g.workel(arg[0]);
    } else {
      // Vector assignment
      g.local("cs", "const casadi_real", "*");
      g.local("rr", "casadi_real", "*");
      g.local("i", "casadi_int");
      g << "for (i=0, rr=" << g.work(res[0], nnz()) << ", cs=" << g.work(arg[0], nnz())
        << "; i<" << sparsity().nnz() << "; ++i) ";
      r = "*rr++";
      x = "*cs++";
    }

    // Output the operation
    g << r << " = " << g.print_op(op_, " " + x + " ") << ";\n";
  }

  MX UnaryMX::get_unary(casadi_int op) const {
    if (!GlobalOptions::simplification_on_the_fly) return MXNode::get_unary(op);

    switch (op_) {
    case OP_NEG:
      if (op==OP_NEG) return dep();
      else if (op==OP_SQ) return dep()->get_unary(OP_SQ);
      else if (op==OP_FABS) return dep()->get_unary(OP_FABS);
      else if (op==OP_COS) return dep()->get_unary(OP_COS);
      break;
    case OP_SQRT:
      if (op==OP_SQ) return dep();
      else if (op==OP_FABS) return shared_from_this<MX>();
      break;
    case OP_SQ:
      if (op==OP_SQRT) return dep()->get_unary(OP_FABS);
      else if (op==OP_FABS) return shared_from_this<MX>();
      break;
    case OP_EXP:
      if (op==OP_LOG) return dep();
      else if (op==OP_FABS) return shared_from_this<MX>();
      break;
    case OP_LOG:
      if (op==OP_EXP) return dep();
      break;
    case OP_FABS:
      if (op==OP_FABS) return shared_from_this<MX>();
      else if (op==OP_SQ) return dep()->get_unary(OP_SQ);
      else if (op==OP_COS) return dep()->get_unary(OP_COS);
      break;
    case OP_INV:
      if (op==OP_INV) return dep();
      break;
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::get_unary(op);
  }

  MX UnaryMX::_get_binary(casadi_int op, const MX& y, bool scX, bool scY) const {
    switch (op_) {
    case OP_NEG:
      if (op==OP_ADD) return y->_get_binary(OP_SUB, dep(), scY, scX);
      else if (op==OP_MUL) return -dep()->_get_binary(OP_MUL, y, scX, scY);
      else if (op==OP_DIV) return -dep()->_get_binary(OP_DIV, y, scX, scY);
      break;
    case OP_INV:
      if (op==OP_MUL) return y->_get_binary(OP_DIV, dep(), scY, scX);
      break;
    case OP_TWICE:
      if (op==OP_SUB && MX::is_equal(y, dep(), maxDepth())) return dep();
      break;
    case OP_SQ:
      if (op==OP_ADD && y.op()==OP_SQ) /*sum of squares:*/
        if ((dep().op()==OP_SIN && y->dep().op()==OP_COS) ||
           (dep().op()==OP_COS && y->dep()->op()==OP_SIN)) /* sin^2(x)+sin^2(y) */
          if (MX::is_equal(dep()->dep(), y->dep()->dep(), maxDepth())) /*sin^2(x) + cos^2(x) */
            return MX::ones(y.sparsity());
      break;
    default: break; // no rule
    }

    // Fallback to default implementation
    return MXNode::_get_binary(op, y, scX, scY);
  }

} // namespace casadi
