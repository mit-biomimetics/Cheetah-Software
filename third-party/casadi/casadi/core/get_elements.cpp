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


#include "get_elements.hpp"
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  MX GetElements::create(const MX& ind, const MX& x) {
    return MX::create(new GetElements(ind, x));
  }

  GetElements::GetElements(const MX& ind, const MX& x) {
    set_sparsity(ind.sparsity());
    set_dep(ind, x);
  }

  int GetElements::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    // Get input and output arguments
    const double* ind = arg[0];
    const double* x = arg[1];
    double* ret = res[0];
    // Dimensions
    casadi_int nnz = dep(0).nnz();
    casadi_int max_ind = dep(1).nnz();
    // If not in-place, copy
    if (ind != ret) casadi_copy(ind, nnz, ret);
    // Get elements
    for (casadi_int i=0; i<nnz; ++i) {
      // Get index
      casadi_int index = static_cast<casadi_int>(*ret);
      // Make assignment if in bounds, else NaN
      *ret++ = index>=0 && index<max_ind ? x[index] : nan;
    }
    return 0;
  }

  std::string GetElements::disp(const std::vector<std::string>& arg) const {
    return arg.at(1) + "(" + arg.at(0) + ")";
  }

} // namespace casadi
