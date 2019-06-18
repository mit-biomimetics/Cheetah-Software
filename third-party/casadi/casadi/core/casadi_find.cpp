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


#include "casadi_find.hpp"

using namespace std;

namespace casadi {

  Find::Find(const MX& x) {
    casadi_assert_dev(x.is_column());
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  std::string Find::disp(const std::vector<std::string>& arg) const {
    return "find(" + arg.at(0) + ")";
  }

  int Find::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    const double* x = arg[0];
    casadi_int nnz = dep(0).nnz();
    casadi_int k=0;
    while (k<nnz && *x++ == 0) k++;
    res[0][0] = k<nnz ? static_cast<double>(dep(0).row(k)) : static_cast<double>(dep(0).size1());
    return 0;
  }

  void Find::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = find(arg[0]);
  }

  void Find::ad_forward(const std::vector<std::vector<MX> >& fseed,
                     std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = 0;
    }
  }

  void Find::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                     std::vector<std::vector<MX> >& asens) const {
  }

  int Find::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    res[0][0] = 0; // pw constant
    return 0;
  }

  int Find::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    res[0][0] = 0; // pw constant
    return 0;
  }

  void Find::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    casadi_int nnz = dep(0).nnz();
    g.local("i", "casadi_int");
    g.local("cr", "const casadi_real", "*");
    g << "for (i=0, cr=" << g.work(arg[0], nnz) << "; i<" << nnz
      << " && *cr++==0; ++i) {}\n"
      << g.workel(res[0]) << " = ";
    if (dep(0).is_dense()) {
      g << "i;\n";
    } else {
      // The row is in position 1+1+2+i (colind has length 2)
      g << "i<" << nnz << " ? " << g.sparsity(dep(0).sparsity()) << "[4+i] : "
        << dep(0).size1() << "\n";
    }
  }

} // namespace casadi
