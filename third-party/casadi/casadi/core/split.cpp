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


#include "split.hpp"
#include "casadi_misc.hpp"
#include "global_options.hpp"

using namespace std;

namespace casadi {

  Split::Split(const MX& x, const std::vector<casadi_int>& offset) : offset_(offset) {
    set_dep(x);
    set_sparsity(Sparsity::scalar());
  }

  Split::~Split() {
  }

  int Split::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int Split::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int Split::eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    // Number of derivatives
    casadi_int nx = offset_.size()-1;

    for (casadi_int i=0; i<nx; ++i) {
      casadi_int nz_first = offset_[i];
      casadi_int nz_last = offset_[i+1];
      if (res[i]!=nullptr) {
        copy(arg[0]+nz_first, arg[0]+nz_last, res[i]);
      }
    }
    return 0;
  }

  int Split::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nx = offset_.size()-1;
    for (casadi_int i=0; i<nx; ++i) {
      if (res[i]!=nullptr) {
        const bvec_t *arg_ptr = arg[0] + offset_[i];
        casadi_int n_i = sparsity(i).nnz();
        bvec_t *res_i_ptr = res[i];
        for (casadi_int k=0; k<n_i; ++k) {
          *res_i_ptr++ = *arg_ptr++;
        }
      }
    }
    return 0;
  }

  int Split::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    casadi_int nx = offset_.size()-1;
    for (casadi_int i=0; i<nx; ++i) {
      if (res[i]!=nullptr) {
        bvec_t *arg_ptr = arg[0] + offset_[i];
        casadi_int n_i = sparsity(i).nnz();
        bvec_t *res_i_ptr = res[i];
        for (casadi_int k=0; k<n_i; ++k) {
          *arg_ptr++ |= *res_i_ptr;
          *res_i_ptr++ = 0;
        }
      }
    }
    return 0;
  }

  void Split::generate(CodeGenerator& g,
                        const std::vector<casadi_int>& arg,
                        const std::vector<casadi_int>& res) const {
    casadi_int nx = nout();
    for (casadi_int i=0; i<nx; ++i) {
      casadi_int nz_first = offset_[i];
      casadi_int nz_last = offset_[i+1];
      casadi_int nz = nz_last-nz_first;
      if (res[i]>=0 && nz>0) { // if anything to assign
        if (nz==1) { // assign scalar
          g << g.workel(res[i]) << " = ";
          if (dep(0).nnz()==1) {
            // rhs is also scalar
            casadi_assert_dev(nz_first==0);
            g << g.workel(arg[0]) << ";\n";
          } else {
            // rhs is an element in a vector
            g << g.work(arg[0], dep(0).nnz()) << "[" << nz_first << "];\n";
          }
        } else {
          // assign vector
          std::string r = g.work(arg[0], dep(0).nnz());
          if (nz_first!=0) r = r + "+" + str(nz_first);
          g << g.copy(r, nz, g.work(res[i], nnz(i))) << "\n";
        }
      }
    }
  }

  Dict Split::info() const {
    std::vector<MX> arg;
    for (auto& sp : output_sparsity_)
      arg.push_back(MX::sym("x", sp));
    Function output("output", std::vector<MX>{}, arg);
    return {{"offset", offset_}, {"output", output}};
  }

  Horzsplit::Horzsplit(const MX& x, const std::vector<casadi_int>& offset) : Split(x, offset) {

    // Split up the sparsity pattern
    output_sparsity_ = horzsplit(x.sparsity(), offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }
  }

  std::string Horzsplit::disp(const std::vector<std::string>& arg) const {
    return "horzsplit(" + arg.at(0) + ")";
  }

  void Horzsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Get column offsets
    vector<casadi_int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    res = horzsplit(arg[0], col_offset);
  }

  void Horzsplit::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    casadi_int nfwd = fsens.size();

    // Get column offsets
    vector<casadi_int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    // Non-differentiated output and forward sensitivities
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d] = horzsplit(fseed[d][0], col_offset);
    }
  }

  void Horzsplit::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    casadi_int nadj = aseed.size();

    // Get column offsets
    vector<casadi_int> col_offset;
    col_offset.reserve(offset_.size());
    col_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      col_offset.push_back(col_offset.back() + s.size2());
    }

    for (casadi_int d=0; d<nadj; ++d) {
      asens[d][0] += horzcat(aseed[d]);
    }
  }

  Diagsplit::Diagsplit(const MX& x,
    const std::vector<casadi_int>& offset1,
    const std::vector<casadi_int>& offset2) : Split(x, offset1) {

    // Split up the sparsity pattern
    output_sparsity_ = diagsplit(x.sparsity(), offset1, offset2);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }

    casadi_assert(offset_.back()==x.nnz(),
      "DiagSplit:: the presence of nonzeros outside the diagonal blocks in unsupported.");
  }

  std::string Diagsplit::disp(const std::vector<std::string>& arg) const {
    return "diagsplit(" + arg.at(0) + ")";
  }

  void Diagsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Get offsets
    vector<casadi_int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<casadi_int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    res = diagsplit(arg[0], offset1, offset2);
  }

  void Diagsplit::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    casadi_int nfwd = fsens.size();
    // Get offsets
    vector<casadi_int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<casadi_int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    // Non-differentiated output and forward sensitivities
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d] = diagsplit(fseed[d][0], offset1, offset2);
    }
  }

  void Diagsplit::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    casadi_int nadj = asens.size();

    // Get offsets
    vector<casadi_int> offset1;
    offset1.reserve(offset_.size());
    offset1.push_back(0);
    vector<casadi_int> offset2;
    offset2.reserve(offset_.size());
    offset2.push_back(0);
    for (auto&& s : output_sparsity_) {
      offset1.push_back(offset1.back() + s.size1());
      offset2.push_back(offset2.back() + s.size2());
    }

    for (casadi_int d=0; d<nadj; ++d) {
      asens[d][0] += diagcat(aseed[d]);
    }
  }

  Vertsplit::Vertsplit(const MX& x, const std::vector<casadi_int>& offset) : Split(x, offset) {

    // Split up the sparsity pattern
    output_sparsity_ = vertsplit(x.sparsity(), offset_);

    // Have offset_ refer to the nonzero offsets instead of column offsets
    offset_.resize(1);
    for (auto&& s : output_sparsity_) {
      offset_.push_back(offset_.back() + s.nnz());
    }
  }

  std::string Vertsplit::disp(const std::vector<std::string>& arg) const {
    return "vertsplit(" + arg.at(0) + ")";
  }

  void Vertsplit::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Get row offsets
    vector<casadi_int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    res = vertsplit(arg[0], row_offset);
  }

  void Vertsplit::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    casadi_int nfwd = fsens.size();

    // Get row offsets
    vector<casadi_int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d] = vertsplit(fseed[d][0], row_offset);
    }
  }

  void Vertsplit::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    casadi_int nadj = aseed.size();

    // Get row offsets
    vector<casadi_int> row_offset;
    row_offset.reserve(offset_.size());
    row_offset.push_back(0);
    for (auto&& s : output_sparsity_) {
      row_offset.push_back(row_offset.back() + s.size1());
    }

    for (casadi_int d=0; d<nadj; ++d) {
      asens[d][0] += vertcat(aseed[d]);
    }
  }

  MX Horzsplit::get_horzcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::get_horzcat(x);
    }

    // Check x content
    for (casadi_int i=0; i<x.size(); ++i) {
      if (!(x[i]->is_output() && x[i]->which_output()==i && x[i]->dep().get()==this)) {
        return MXNode::get_horzcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

  MX Vertsplit::get_vertcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::get_vertcat(x);
    }

    // Check x content
    for (casadi_int i=0; i<x.size(); ++i) {
      if (!(x[i]->is_output() && x[i]->which_output()==i && x[i]->dep().get()==this)) {
        return MXNode::get_vertcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

  MX Diagsplit::get_diagcat(const std::vector<MX>& x) const {
    // Check x length
    if (x.size()!=nout()) {
      return MXNode::get_diagcat(x);
    }

    // Check x content
    for (casadi_int i=0; i<x.size(); ++i) {
      if (!(x[i]->is_output() && x[i]->which_output()==i && x[i]->dep().get()==this)) {
        return MXNode::get_diagcat(x);
      }
    }

    // OK if reached this point
    return dep();
  }

} // namespace casadi
