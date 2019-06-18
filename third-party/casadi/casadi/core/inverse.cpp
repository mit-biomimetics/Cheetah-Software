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


#include "inverse.hpp"

using namespace std;

namespace casadi {

  Inverse::Inverse(const MX& x) {
    casadi_assert(x.size1()==x.size2(),
    "Inverse: matrix must be square, but you supllied " + x.dim());
    set_dep(x);
    set_sparsity(Sparsity::dense(x.size1(), x.size2()));
  }

  std::string Inverse::disp(const std::vector<std::string>& arg) const {
    return "inv(" + arg.at(0) + ")";
  }

  void Inverse::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = inv(arg[0]);
  }

  void Inverse::ad_forward(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) const {
    MX inv_X = shared_from_this<MX>();
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = -mtimes(inv_X, mtimes(fseed[d][0], inv_X));
    }
  }

  void Inverse::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    MX inv_X = shared_from_this<MX>();
    MX trans_inv_X = inv_X.T();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] -= mtimes(trans_inv_X, mtimes(aseed[d][0], trans_inv_X));
    }
  }

} // namespace casadi
