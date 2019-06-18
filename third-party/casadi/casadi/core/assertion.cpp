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


#include "assertion.hpp"

using namespace std;

namespace casadi {

  Assertion::Assertion(const MX& x, const MX& y, const std::string & fail_message)
      : fail_message_(fail_message) {
    casadi_assert(y.is_scalar(),
      "Assertion:: assertion expression y must be scalar, but got " + y.dim());
    set_dep(x, y);
    set_sparsity(x.sparsity());
  }

  std::string Assertion::disp(const std::vector<std::string>& arg) const {
    return "assertion(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void Assertion::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0].attachAssert(arg[1], fail_message_);
  }

  void Assertion::ad_forward(const std::vector<std::vector<MX> >& fseed,
                     std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = fseed[d][0];
    }
  }

  void Assertion::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                     std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0];
    }
  }

  int Assertion::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Assertion::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[1][0]!=1) {
      casadi_error("Assertion error: " + fail_message_);
      return 1;
    }

    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Assertion::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    if (arg[0]!=res[0]) {
      copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Assertion::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    casadi_int n = nnz();
    if (a != r) {
      for (casadi_int i=0; i<n; ++i) {
        *a++ |= *r;
        *r++ = 0;
      }
    }
    return 0;
  }

  void Assertion::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    // Generate assertion
    g << "if (" << g.workel(arg[1]) << "!=1.) {\n"
      << "    /* " << fail_message_ << " */\n"
      << "    return 1;\n"
      << "  }\n";

    // Copy if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << '\n';
    }
  }

} // namespace casadi
