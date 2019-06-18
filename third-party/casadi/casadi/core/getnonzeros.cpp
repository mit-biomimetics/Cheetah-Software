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


#include "getnonzeros.hpp"
#include "casadi_misc.hpp"

using namespace std;

namespace casadi {

  MX GetNonzeros::create(const Sparsity& sp, const MX& x, const std::vector<casadi_int>& nz) {
    // No elements at all
    if (nz.size()==0) return MX::zeros(sp);
    // Simplify to slice
    if (is_slice(nz)) return create(sp, x, to_slice(nz));
    // Simplify to slice2
    if (is_slice2(nz)) {
      pair<Slice, Slice> sl = to_slice2(nz);
      return create(sp, x, sl.first, sl.second);
    }
    return MX::create(new GetNonzerosVector(sp, x, nz));
  }

  MX GetNonzeros::create(const Sparsity& sp, const MX& x, const Slice& s) {
    // Simplify identity assignments
    if (sp==x.sparsity() && s.start==0 && s.step==1 && s.stop==x.nnz()) return x;
    return MX::create(new GetNonzerosSlice(sp, x, s));
  }

  MX GetNonzeros::create(const Sparsity& sp, const MX& x,
                         const Slice& inner, const Slice& outer) {
    return MX::create(new GetNonzerosSlice2(sp, x, inner, outer));
  }

  GetNonzeros::GetNonzeros(const Sparsity& sp, const MX& y) {
    set_sparsity(sp);
    set_dep(y);
  }

  int GetNonzerosVector::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int GetNonzerosVector::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int GetNonzerosVector::
  eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const {
    const T* idata = arg[0];
    T* odata = res[0];
    for (auto&& k : nz_) {
      *odata++ = k>=0 ? idata[k] : 0;
    }
    return 0;
  }

  int GetNonzerosSlice::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int GetNonzerosSlice::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int GetNonzerosSlice::eval_gen(const T* const* arg, T* const* res,
                                 casadi_int* iw, T* w) const {
    const T* idata = arg[0] + s_.start;
    const T* idata_stop = arg[0] + s_.stop;
    T* odata = res[0];
    for (; idata != idata_stop; idata += s_.step) {
      *odata++ = *idata;
    }
    return 0;
  }

  int GetNonzerosSlice2::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  int GetNonzerosSlice2::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<typename T>
  int GetNonzerosSlice2::
  eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const {
    const T* outer = arg[0] + outer_.start;
    const T* outer_stop = arg[0] + outer_.stop;
    T* odata = res[0];
    for (; outer != outer_stop; outer += outer_.step) {
      for (const T* inner = outer+inner_.start;
          inner != outer+inner_.stop;
          inner += inner_.step) {
        *odata++ = *inner;
      }
    }
    return 0;
  }

  int GetNonzerosVector::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (auto&& k : nz_) *r++ = k>=0 ? a[k] : 0;
    return 0;
  }

  int GetNonzerosVector::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (auto&& k : nz_) {
      if (k>=0) a[k] |= *r;
      *r++ = 0;
    }
    return 0;
  }

  int GetNonzerosSlice::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      *r++ = a[k];
    }
    return 0;
  }

  int GetNonzerosSlice::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      a[k] |= *r;
      *r++ = 0;
    }
    return 0;
  }

  int GetNonzerosSlice2::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        *r++ = a[k2];
      }
    }
    return 0;
  }

  int GetNonzerosSlice2::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        a[k2] |= *r;
        *r++ = 0;
      }
    }
    return 0;
  }

  std::string GetNonzerosVector::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << nz_;
    return ss.str();
  }

  std::string GetNonzerosSlice::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << s_ << "]";
    return ss.str();
  }

  std::string GetNonzerosSlice2::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << arg.at(0) << "[" << outer_ << ";" << inner_ << "]";
    return ss.str();
  }

  void GetNonzeros::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Output sparsity
    const Sparsity& osp = sparsity();
    const casadi_int* orow = osp.row();
    vector<casadi_int> ocol = osp.get_col();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<casadi_int>& irow = isp.row();
    vector<casadi_int> icol = isp.get_col();

    // Get all input elements
    vector<casadi_int> el_input;
    isp.find(el_input);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    // Get the matching nonzeros
    r_ind.resize(el_input.size());
    copy(el_input.begin(), el_input.end(), r_ind.begin());
    arg[0].sparsity().get_nz(r_ind);

    // Sparsity pattern for the result
    r_colind.resize(osp.size2()+1); // Col count
    fill(r_colind.begin(), r_colind.end(), 0);
    r_row.clear();

    // Perform the assignments
    r_nz.clear();
    for (casadi_int k=0; k<nz.size(); ++k) {

      // Get the corresponding nonzero for the input
      casadi_int el = nz[k];

      // Skip if zero assignment
      if (el==-1) continue;

      // Get the corresponding nonzero in the argument
      casadi_int el_arg = r_ind[el];

      // Skip if no argument
      if (el_arg==-1) continue;

      // Save the assignment
      r_nz.push_back(el_arg);

      // Get the corresponding element
      casadi_int i=ocol[k], j=orow[k];

      // Add to sparsity pattern
      r_row.push_back(j);
      r_colind[1+i]++;
    }

    // col count -> col offset
    for (casadi_int i=1; i<r_colind.size(); ++i) r_colind[i] += r_colind[i-1];

    // Create a sparsity pattern from vectors
    if (r_nz.size()==0) {
      res[0] = MX(osp.size());
    } else {
      Sparsity f_sp(osp.size1(), osp.size2(), r_colind, r_row);
      res[0] = arg[0]->get_nzref(f_sp, r_nz);
    }
  }

  void GetNonzeros::ad_forward(const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens) const {

    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Number of derivative directions
    casadi_int nfwd = fsens.size();

    // Output sparsity
    const Sparsity& osp = sparsity();
    const casadi_int* orow = osp.row();
    vector<casadi_int> ocol = osp.get_col();

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<casadi_int>& irow = isp.row();
    vector<casadi_int> icol;

    // Get all input elements
    vector<casadi_int> el_input;

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = fseed[d][0];
      MX& res = fsens[d][0];

      if (arg.sparsity()==isp) { // Matching sparsity
        if (nz.size()==0) {
          res = MX(osp.size());
        } else {
          res = arg->get_nzref(osp, nz);
        }
      } else {
        // Expensive operations
        if (el_input.empty()) isp.find(el_input);
        if (icol.empty()) icol = isp.get_col();

        // Get the matching nonzeros
        r_ind.resize(el_input.size());
        copy(el_input.begin(), el_input.end(), r_ind.begin());
        arg.sparsity().get_nz(r_ind);

        // Sparsity pattern for the result
        r_colind.resize(osp.size2()+1); // Col count
        fill(r_colind.begin(), r_colind.end(), 0);
        r_row.clear();

        // Perform the assignments
        r_nz.clear();
        for (casadi_int k=0; k<nz.size(); ++k) {

          // Get the corresponding nonzero for the input
          casadi_int el = nz[k];

          // Skip if zero assignment
          if (el==-1) continue;

          // Get the corresponding nonzero in the argument
          casadi_int el_arg = r_ind[el];

          // Skip if no argument
          if (el_arg==-1) continue;

          // Save the assignment
          r_nz.push_back(el_arg);

          // Get the corresponding element
          casadi_int i=ocol[k], j=orow[k];

          // Add to sparsity pattern
          r_row.push_back(j);
          r_colind[1+i]++;
        }

        // col count -> col offset
        for (casadi_int i=1; i<r_colind.size(); ++i) r_colind[i] += r_colind[i-1];

        // Create a sparsity pattern from vectors
        if (r_nz.size()==0) {
          res = MX(osp.size());
        } else {
          Sparsity f_sp(osp.size1(), osp.size2(), r_colind, r_row);
          res = arg->get_nzref(f_sp, r_nz);
        }
      }
    }
  }

  void GetNonzeros::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens) const {
    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Number of derivative directions
    casadi_int nadj = aseed.size();

    // Output sparsity
    const Sparsity& osp = sparsity();
    vector<casadi_int> ocol;

    // Input sparsity
    const Sparsity& isp = dep().sparsity();
    //const vector<casadi_int>& irow = isp.row();
    vector<casadi_int> icol;

    // Get all input elements
    vector<casadi_int> el_input;

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    // Adjoint sensitivities
    for (casadi_int d=0; d<nadj; ++d) {

      // Get an owning references to the seeds and sensitivities
      // and clear the seeds for the next run
      MX aseed0 = aseed[d][0];
      MX asens0 = asens[d][0]; // Sensitivity before addition

      if (aseed0.sparsity()==osp && asens0.sparsity().nnz()==0) { // Matching sparsity
        asens[d][0] = aseed0->get_nzadd(DM::zeros(isp), nz);
      } else {
        // Expensive operations
        if (el_input.empty()) isp.find(el_input);
        if (icol.empty()) icol = isp.get_col();
        if (ocol.empty()) ocol = osp.get_col();

        // Get the corresponding nz locations in the output sparsity pattern
        aseed0.sparsity().find(r_nz);
        osp.get_nz(r_nz);

        // Filter out ignored entries and check if there is anything to add at all
        bool elements_to_add = false;
        for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
          if (*k>=0) {
            if (nz[*k]>=0) {
              elements_to_add = true;
            } else {
              *k = -1;
            }
          }
        }

        // Quick continue of no elements to add
        if (!elements_to_add) continue;

        // Get the nz locations in the adjoint sensitivity corresponding to the inputs
        r_ind.resize(el_input.size());
        copy(el_input.begin(), el_input.end(), r_ind.begin());
        asens0.sparsity().get_nz(r_ind);

        // Enlarge the sparsity pattern of the sensitivity if not all additions fit
        for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
          if (*k>=0 && r_ind[nz[*k]]<0) {

            // Create a new pattern which includes both the the previous seed and the addition
            Sparsity sp = asens0.sparsity().unite(dep().sparsity());
            asens0 = asens0->get_project(sp);

            // Recalculate the nz locations in the adjoint sensitivity corresponding to the inputs
            copy(el_input.begin(), el_input.end(), r_ind.begin());
            asens0.sparsity().get_nz(r_ind);

            break;
          }
        }

        // Have r_nz point to locations in the sensitivity instead of the output
        for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
          if (*k>=0) {
            *k = r_ind[nz[*k]];
          }
        }

        asens[d][0] = aseed0->get_nzadd(asens0, r_nz);
      }
    }
  }

  Matrix<casadi_int> GetNonzeros::mapping() const {
    vector<casadi_int> nz = all();
    return Matrix<casadi_int>(sparsity(), nz, false);
  }

  void GetNonzerosVector::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
    // Codegen the indices
    string ind = g.constant(nz_);

    // Codegen the assignments
    g.local("cii", "const casadi_int", "*");
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g << "for (cii=" << ind << ", rr=" << g.work(res[0], nnz())
      << ", ss=" << g.work(arg[0], dep(0).nnz())
      << "; cii!=" << ind << "+" << nz_.size()
      << "; ++cii) *rr++ = *cii>=0 ? ss[*cii] : 0;\n";
  }

  MX GetNonzeros::get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const {
    // Get all the nonzeros
    vector<casadi_int> nz_all = all();

    // Eliminate recursive calls
    vector<casadi_int> nz_new(nz);
    for (vector<casadi_int>::iterator i=nz_new.begin(); i!=nz_new.end(); ++i) {
      if (*i>=0) *i = nz_all[*i];
    }
    return dep()->get_nzref(sp, nz_new);
  }

  void GetNonzerosSlice::generate(CodeGenerator& g,
                                  const std::vector<casadi_int>& arg,
                                  const std::vector<casadi_int>& res) const {
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g << "for (rr=" << g.work(res[0], nnz()) << ", ss=" << g.work(arg[0], dep(0).nnz())
      << "+" << s_.start << "; ss!=" << g.work(arg[0], dep(0).nnz()) << "+" << s_.stop
      << "; ss+=" << s_.step << ") *rr++ = *ss;\n";
  }

  void GetNonzerosSlice2::generate(CodeGenerator& g,
                                    const std::vector<casadi_int>& arg,
                                    const std::vector<casadi_int>& res) const {
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g.local("tt", "casadi_real", "*");
    g << "for (rr=" << g.work(res[0], nnz()) << ", ss=" << g.work(arg[0], dep(0).nnz())
      << "+" << outer_.start << "; ss!=" << g.work(arg[0], dep(0).nnz()) << "+"
      << outer_.stop << "; ss+=" << outer_.step << ") "
      << "for (tt=ss+" << inner_.start << "; tt!=ss+" << inner_.stop
      << "; tt+=" << inner_.step << ") *rr++ = *tt;\n";
  }

  bool GetNonzerosVector::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosVector* n = dynamic_cast<const GetNonzerosVector*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->nz_.size()!=n->nz_.size()) return false;
    if (!std::equal(this->nz_.begin(), this->nz_.end(), n->nz_.begin())) return false;

    return true;
  }

  bool GetNonzerosSlice::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosSlice* n = dynamic_cast<const GetNonzerosSlice*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->s_ != n->s_) return false;

    return true;
  }

  bool GetNonzerosSlice2::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const GetNonzerosSlice2* n = dynamic_cast<const GetNonzerosSlice2*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->inner_ != n->inner_ || this->outer_!=n->outer_) return false;

    return true;
  }


} // namespace casadi
