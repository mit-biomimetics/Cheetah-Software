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


#include "determinant.hpp"

using namespace std;

namespace casadi {

  Determinant::Determinant(const MX& x) {
    casadi_assert(x.is_square(), "Dimension mismatch. Matrix must be square, "
      "but got " + x.dim() + " instead.");
    set_dep(x);
    set_sparsity(Sparsity::dense(1, 1));
  }

  std::string Determinant::disp(const std::vector<std::string>& arg) const {
    return "det(" + arg.at(0) + ")";
  }

  void Determinant::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = det(arg[0]);
  }

  void Determinant::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {
    const MX& X = dep();
    MX det_X = shared_from_this<MX>();
    MX trans_inv_X = inv(X).T();
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = det_X * dot(trans_inv_X, fseed[d][0]);
    }
  }

  void Determinant::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    const MX& X = dep();
    MX det_X = shared_from_this<MX>();
    MX trans_inv_X = inv(X).T();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += aseed[d][0]*det_X * trans_inv_X;
    }
  }

} // namespace casadi
