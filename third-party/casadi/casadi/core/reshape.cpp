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


#include "reshape.hpp"
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  Reshape::Reshape(const MX& x, Sparsity sp) {
    casadi_assert_dev(x.nnz()==sp.nnz());
    set_dep(x);
    set_sparsity(sp);
  }

  int Reshape::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Reshape::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Reshape::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    if (arg[0]!=res[0]) copy(arg[0], arg[0]+nnz(), res[0]);
    return 0;
  }

  int Reshape::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_fwd(arg[0], res[0], nnz());
    return 0;
  }

  int Reshape::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_rev(arg[0], res[0], nnz());
    return 0;
  }

  std::string Reshape::disp(const std::vector<std::string>& arg) const {
    // For vectors, reshape is also a transpose
    if (dep().is_vector() && sparsity().is_vector()) {
      // Print as transpose: X'
      return arg.at(0) + "'";
    } else {
      // Print as reshape(X) or vec(X)
      if (sparsity().is_column()) {
        return "vec(" + arg.at(0) + ")";
      } else {
        return "reshape(" + arg.at(0) + ")";
      }
    }
  }

  void Reshape::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = reshape(arg[0], size());
  }

  void Reshape::ad_forward(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d = 0; d<fsens.size(); ++d) {
      fsens[d][0] = reshape(fseed[d][0], size());
    }
  }

  void Reshape::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += reshape(aseed[d][0], dep().size());
    }
  }

  void Reshape::generate(CodeGenerator& g,
                         const std::vector<casadi_int>& arg,
                         const std::vector<casadi_int>& res) const {
    if (arg[0]==res[0]) return;
    g << g.copy(g.work(arg[0], nnz()), nnz(), g.work(res[0], nnz())) << "\n";
  }

  MX Reshape::get_reshape(const Sparsity& sp) const {
    return reshape(dep(0), sp);
  }

  MX Reshape::get_transpose() const {
    // For vectors, reshape is also a transpose
    if (dep().is_vector() && sparsity().is_vector()) {
      return dep();
    } else {
      return MXNode::get_transpose();
    }
  }

  bool Reshape::is_valid_input() const {
    if (!dep()->is_valid_input()) return false;
    return true;
  }

  casadi_int Reshape::n_primitives() const {
    return dep()->n_primitives();
  }

  void Reshape::primitives(std::vector<MX>::iterator& it) const {
    dep()->primitives(it);
  }

  void Reshape::split_primitives(const MX& x, std::vector<MX>::iterator& it) const {
    dep()->split_primitives(reshape(x, dep().size()), it);
  }

  MX Reshape::join_primitives(std::vector<MX>::const_iterator& it) const {
    return reshape(dep()->join_primitives(it), size());
  }

  bool Reshape::has_duplicates() const {
    return dep()->has_duplicates();
  }

  void Reshape::reset_input() const {
    dep()->reset_input();
  }

} // namespace casadi
