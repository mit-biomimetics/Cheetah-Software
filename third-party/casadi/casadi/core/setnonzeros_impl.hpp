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


#ifndef CASADI_SETNONZEROS_IMPL_HPP
#define CASADI_SETNONZEROS_IMPL_HPP

#include "setnonzeros.hpp"
#include "casadi_misc.hpp"

/// \cond INTERNAL

using namespace std;

namespace casadi {

  template<bool Add>
  MX SetNonzeros<Add>::create(const MX& y, const MX& x, const std::vector<casadi_int>& nz) {
    if (is_slice(nz)) return create(y, x, to_slice(nz));
    if (is_slice2(nz)) {
      pair<Slice, Slice> sl = to_slice2(nz);
      return create(y, x, sl.first, sl.second);
    }
    return MX::create(new SetNonzerosVector<Add>(y, x, nz));
  }

  template<bool Add>
  MX SetNonzeros<Add>::create(const MX& y, const MX& x, const Slice& s) {
    // Simplify if assignment
    if (y.sparsity()==x.sparsity() && s.start==0 && s.step==1 && s.stop==x.nnz()) {
      if (Add) {
        return y + x;
      } else {
        return x;
      }
    }
    return MX::create(new SetNonzerosSlice<Add>(y, x, s));
  }

  template<bool Add>
  MX SetNonzeros<Add>::create(const MX& y, const MX& x, const Slice& inner, const Slice& outer) {
    return MX::create(new SetNonzerosSlice2<Add>(y, x, inner, outer));
  }

  template<bool Add>
  SetNonzeros<Add>::SetNonzeros(const MX& y, const MX& x) {
    this->set_sparsity(y.sparsity());
    this->set_dep(y, x);
  }

  template<bool Add>
  SetNonzerosVector<Add>::SetNonzerosVector(const MX& y, const MX& x,
      const std::vector<casadi_int>& nz) : SetNonzeros<Add>(y, x), nz_(nz) {
    // Ignore duplicate assignments
    if (!Add) {
      vector<bool> already_set(this->nnz(), false);
      for (vector<casadi_int>::reverse_iterator i=nz_.rbegin(); i!=nz_.rend(); ++i) {
        if (*i>=0) {
          if (already_set[*i]) {
            *i = -1;
          } else {
            already_set[*i] = true;
          }
        }
      }
    }
  }

  template<bool Add>
  SetNonzeros<Add>:: ~SetNonzeros() {
  }

  template<bool Add>
  void SetNonzeros<Add>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const casadi_int* orow = osp.row();
    vector<casadi_int> ocol = osp.get_col();

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    vector<casadi_int> icol = isp.get_col();

    // We next need to resort the assignment vector by outputs instead of inputs
    // Start by counting the number of output nonzeros corresponding to each input nonzero
    vector<casadi_int> onz_count(osp.nnz()+2, 0);
    for (vector<casadi_int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
      onz_count[*it+2]++;
    }

    // Cumsum to get index offset for output nonzero
    for (casadi_int i=0; i<onz_count.size()-1; ++i) {
      onz_count[i+1] += onz_count[i];
    }

    // Get the order of assignments
    vector<casadi_int> nz_order(nz.size());
    for (casadi_int k=0; k<nz.size(); ++k) {
      // Save the new index
      nz_order[onz_count[1+nz[k]]++] = k;
    }

    // Find out which elements are being set
    vector<casadi_int>& with_duplicates = onz_count; // Reuse memory
    onz_count.resize(nz.size());
    for (casadi_int k=0; k<nz.size(); ++k) {
      // Get output nonzero
      casadi_int onz_k = nz[nz_order[k]];

      // Get element (note: may contain duplicates)
      if (onz_k>=0) {
        with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
      } else {
        with_duplicates[k] = -1;
      }
    }

    // Get all output elements (this time without duplicates)
    vector<casadi_int> el_output;
    osp.find(el_output);

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    // Get references to arguments and results
    res[0] = arg[0];

    // Entries in res with elements zero'ed out
    if (!Add) {

      // Get the nz locations in res corresponding to the output sparsity pattern
      r_nz.resize(with_duplicates.size());
      copy(with_duplicates.begin(), with_duplicates.end(), r_nz.begin());
      res[0].sparsity().get_nz(r_nz);

      // Zero out the corresponding entries
      res[0] = MX::zeros(isp)->get_nzassign(res[0], r_nz);
    }

    // Get the nz locations of the elements in arg corresponding to the argument sparsity pattern
    arg[1].sparsity().find(r_nz);
    isp.get_nz(r_nz);

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

    // Quick continue of no elements to set/add
    if (!elements_to_add) return;

    // Get the nz locations in the argument corresponding to the inputs
    r_ind.resize(el_output.size());
    copy(el_output.begin(), el_output.end(), r_ind.begin());
    res[0].sparsity().get_nz(r_ind);

    // Enlarge the sparsity pattern of the arguments if not all assignments fit
    for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
      if (*k>=0 && nz[*k]>=0 && r_ind[nz[*k]]<0) {

        // Create a new pattern which includes both the the previous seed
        // and the addition/assignment
        Sparsity sp = res[0].sparsity().unite(osp);
        res[0] = res[0]->get_project(sp);

        // Recalculate the nz locations in the arguments corresponding to the inputs
        copy(el_output.begin(), el_output.end(), r_ind.begin());
        res[0].sparsity().get_nz(r_ind);

        break;
      }
    }

    // Have r_nz point to locations in the result instead of the output
    for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
      if (*k>=0) {
        *k = r_ind[nz[*k]];
      }
    }

    // Add to the element to the sensitivity, if any
    res[0] = arg[1]->get_nzadd(res[0], r_nz);
  }

  template<bool Add>
  void SetNonzeros<Add>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                 std::vector<std::vector<MX> >& fsens) const {
    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Number of derivative directions
    casadi_int nfwd = fsens.size();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const casadi_int* orow = osp.row();
    vector<casadi_int> ocol;

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    vector<casadi_int> icol;

    bool first_run = true;

    vector<casadi_int> onz_count;

    vector<casadi_int> nz_order;

    // Find out which elements are being set
    vector<casadi_int>& with_duplicates = onz_count;

    // Get all output elements (this time without duplicates)
    vector<casadi_int> el_output;

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    // Nondifferentiated function and forward sensitivities
    for (casadi_int d=0; d<nfwd; ++d) {

      // Get references to arguments and results
      const MX& arg = fseed[d][1];
      const MX& arg0 = fseed[d][0];

      MX& res = fsens[d][0];
      res = arg0;

      if (isp==arg.sparsity() && arg0.sparsity()==osp) {
        /*
          dep(0) <-> y
          dep(1) <-> x

          y[nz]+=x

          dot(y)[nz]+=dot(x)

          dot(x)->get_nzadd(dot(y), nz)

        */
        if (!Add) {
          // Zero out the corresponding entries
          res = MX::zeros(isp)->get_nzassign(res, nz);
        }

        res = arg->get_nzadd(res, nz);
      } else {
        if (first_run) {
          osp.find(el_output);
          ocol = osp.get_col();
          icol = isp.get_col();

          // We next need to resort the assignment vector by outputs instead of inputs
          // Start by counting the number of output nonzeros corresponding to each input nonzero
          onz_count.resize(osp.nnz()+2, 0);
          for (vector<casadi_int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
            onz_count[*it+2]++;
          }

          // Cumsum to get index offset for output nonzero
          for (casadi_int i=0; i<onz_count.size()-1; ++i) {
            onz_count[i+1] += onz_count[i];
          }

          // Get the order of assignments
          nz_order.resize(nz.size());
          for (casadi_int k=0; k<nz.size(); ++k) {
            // Save the new index
            nz_order[onz_count[1+nz[k]]++] = k;
          }

          onz_count.resize(nz.size());
          for (casadi_int k=0; k<nz.size(); ++k) {
            // Get output nonzero
            casadi_int onz_k = nz[nz_order[k]];

            // Get element (note: may contain duplicates)
            if (onz_k>=0) {
              with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
            } else {
              with_duplicates[k] = -1;
            }
          }


          first_run = false;
        }

        // Entries in res with elements zero'ed out
        if (!Add) {

          // Get the nz locations in res corresponding to the output sparsity pattern
          r_nz.resize(with_duplicates.size());
          copy(with_duplicates.begin(), with_duplicates.end(), r_nz.begin());
          res.sparsity().get_nz(r_nz);

          // Zero out the corresponding entries
          res = MX::zeros(isp)->get_nzassign(res, r_nz);
        }

        // Get the nz locations of the elements in arg corresponding to the argument
        // sparsity pattern
        arg.sparsity().find(r_nz);
        isp.get_nz(r_nz);

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

        // Quick continue of no elements to set/add
        if (!elements_to_add) continue;

        // Get the nz locations in the argument corresponding to the inputs
        r_ind.resize(el_output.size());
        copy(el_output.begin(), el_output.end(), r_ind.begin());
        res.sparsity().get_nz(r_ind);

        // Enlarge the sparsity pattern of the arguments if not all assignments fit
        for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
          if (*k>=0 && nz[*k]>=0 && r_ind[nz[*k]]<0) {

            // Create a new pattern which includes both the the previous seed
            // and the addition/assignment
            Sparsity sp = res.sparsity().unite(osp);
            res = res->get_project(sp);

            // Recalculate the nz locations in the arguments corresponding to the inputs
            copy(el_output.begin(), el_output.end(), r_ind.begin());
            res.sparsity().get_nz(r_ind);

            break;
          }
        }

        // Have r_nz point to locations in the result instead of the output
        for (vector<casadi_int>::iterator k=r_nz.begin(); k!=r_nz.end(); ++k) {
          if (*k>=0) {
            *k = r_ind[nz[*k]];
          }
        }

        // Add to the element to the sensitivity, if any
        res = arg->get_nzadd(res, r_nz);
      }
    }
  }

  template<bool Add>
  void SetNonzeros<Add>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                 std::vector<std::vector<MX> >& asens) const {
    // Get all the nonzeros
    vector<casadi_int> nz = all();

    // Number of derivative directions
    casadi_int nadj = aseed.size();

    // Output sparsity
    const Sparsity &osp = sparsity();
    const casadi_int* orow = osp.row();
    vector<casadi_int> ocol;

    // Input sparsity (first input same as output)
    const Sparsity &isp = dep(1).sparsity();
    const casadi_int* irow = isp.row();
    vector<casadi_int> icol;

    vector<casadi_int> onz_count;

    // Get the order of assignments
    vector<casadi_int> nz_order;

    vector<casadi_int>& with_duplicates = onz_count; // Reuse memory

    // Get all output elements (this time without duplicates)
    vector<casadi_int> el_output;

    bool first_run = true;

    // Sparsity pattern being formed and corresponding nonzero mapping
    vector<casadi_int> r_colind, r_row, r_nz, r_ind;

    for (casadi_int d=0; d<nadj; ++d) {
      if (osp==aseed[d][0].sparsity()) {
        /*
          dep(0) <-> y
          dep(1) <-> x

          z: y[nz]+=x

          bar(x) += bar(z)[nz]
          bar(y) += bar(z)
        */
        asens[d][1] += aseed[d][0]->get_nzref(isp, nz);
        if (!Add) {
          asens[d][0] += MX::zeros(isp)->get_nzassign(aseed[d][0], nz);
        } else {
          asens[d][0] += aseed[d][0];
        }
      } else {
        if (first_run) {
          ocol = osp.get_col();
          icol = isp.get_col();
          // We next need to resort the assignment vector by outputs instead of inputs
          // Start by counting the number of output nonzeros corresponding to each input nonzero
          onz_count.resize(osp.nnz()+2, 0);
          for (vector<casadi_int>::const_iterator it=nz.begin(); it!=nz.end(); ++it) {
            onz_count[*it+2]++;
          }

          // Cumsum to get index offset for output nonzero
          for (casadi_int i=0; i<onz_count.size()-1; ++i) {
            onz_count[i+1] += onz_count[i];
          }

          // Get the order of assignments
          nz_order.resize(nz.size());
          for (casadi_int k=0; k<nz.size(); ++k) {
            // Save the new index
            nz_order[onz_count[1+nz[k]]++] = k;
          }

          // Find out which elements are being set
          onz_count.resize(nz.size());
          for (casadi_int k=0; k<nz.size(); ++k) {
            // Get output nonzero
            casadi_int onz_k = nz[nz_order[k]];

            // Get element (note: may contain duplicates)
            if (onz_k>=0) {
              with_duplicates[k] = ocol[onz_k]*osp.size1() + orow[onz_k];
            } else {
              with_duplicates[k] = -1;
            }
          }

          osp.find(el_output);
          first_run = false;
        }

        // Get the matching nonzeros
        r_ind.resize(el_output.size());
        copy(el_output.begin(), el_output.end(), r_ind.begin());
        aseed[d][0].sparsity().get_nz(r_ind);

        // Sparsity pattern for the result
        r_colind.resize(isp.size2()+1); // Col count
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
          casadi_int i=icol[k], j=irow[k];

          // Add to sparsity pattern
          r_row.push_back(j);
          r_colind[1+i]++;
        }

        // col count -> col offset
        for (casadi_int i=1; i<r_colind.size(); ++i) r_colind[i] += r_colind[i-1];

        // If anything to set/add
        if (!r_nz.empty()) {
          // Create a sparsity pattern from vectors
          Sparsity f_sp(isp.size1(), isp.size2(), r_colind, r_row);
          asens[d][1] += aseed[d][0]->get_nzref(f_sp, r_nz);
          if (!Add) {
            asens[d][0] += MX::zeros(f_sp)->get_nzassign(aseed[d][0], r_nz);
          } else {
            asens[d][0] += aseed[d][0];
          }
        } else {
          asens[d][0] += aseed[d][0];
        }
      }
    }
  }

  template<bool Add>
  int SetNonzerosVector<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  template<bool Add>
  int SetNonzerosVector<Add>::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<bool Add>
  template<typename T>
  int SetNonzerosVector<Add>::
  eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    const T* idata0 = arg[0];
    const T* idata = arg[1];
    T* odata = res[0];
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    for (vector<casadi_int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++idata) {
      if (Add) {
        if (*k>=0) odata[*k] += *idata;
      } else {
        if (*k>=0) odata[*k] = *idata;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  template<bool Add>
  int SetNonzerosSlice<Add>::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<bool Add>
  template<typename T>
  int SetNonzerosSlice<Add>::
  eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    const T* idata0 = arg[0];
    const T* idata = arg[1];
    T* odata = res[0];
    if (idata0 != odata) {
      copy(idata0, idata0+this->dep(0).nnz(), odata);
    }
    T* odata_stop = odata + s_.stop;
    for (odata += s_.start; odata != odata_stop; odata += s_.step) {
      if (Add) {
        *odata += *idata++;
      } else {
        *odata = *idata++;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice2<Add>::
  eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    return eval_gen<double>(arg, res, iw, w);
  }

  template<bool Add>
  int SetNonzerosSlice2<Add>::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    return eval_gen<SXElem>(arg, res, iw, w);
  }

  template<bool Add>
  template<typename T>
  int SetNonzerosSlice2<Add>::
  eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
    const T* idata0 = arg[0];
    const T* idata = arg[1];
    T* odata = res[0];
    if (idata0 != odata) {
      copy(idata0, idata0 + this->dep(0).nnz(), odata);
    }
    T* outer_stop = odata + outer_.stop;
    T* outer = odata + outer_.start;
    for (; outer != outer_stop; outer += outer_.step) {
      for (T* inner = outer+inner_.start;
          inner != outer+inner_.stop;
          inner += inner_.step) {
        if (Add) {
          *inner += *idata++;
        } else {
          *inner = *idata++;
        }
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosVector<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    casadi_int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (vector<casadi_int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++a) {
      if (Add) {
        if (*k>=0) r[*k] |= *a;
      } else {
        if (*k>=0) r[*k] = *a;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosVector<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (vector<casadi_int>::const_iterator k=this->nz_.begin(); k!=this->nz_.end(); ++k, ++a) {
      if (*k>=0) {
        *a |= r[*k];
        if (!Add) {
          r[*k] = 0;
        }
      }
    }
    MXNode::copy_rev(arg[0], r, this->nnz());
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    casadi_int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      if (Add) {
        r[k] |= *a++;
      } else {
        r[k] = *a++;
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (casadi_int k=s_.start; k!=s_.stop; k+=s_.step) {
      *a++ |= r[k];
      if (!Add) {
        r[k] = 0;
      }
    }
    MXNode::copy_rev(arg[0], r, this->nnz());
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice2<Add>::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    const bvec_t *a0 = arg[0];
    const bvec_t *a = arg[1];
    bvec_t *r = res[0];
    casadi_int n = this->nnz();

    // Propagate sparsity
    if (r != a0) copy(a0, a0+n, r);
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        if (Add) {
          r[k2] |= *a++;
        } else {
          r[k2] = *a++;
        }
      }
    }
    return 0;
  }

  template<bool Add>
  int SetNonzerosSlice2<Add>::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[1];
    bvec_t *r = res[0];
    for (casadi_int k1=outer_.start; k1!=outer_.stop; k1+=outer_.step) {
      for (casadi_int k2=k1+inner_.start; k2!=k1+inner_.stop; k2+=inner_.step) {
        *a++ |= r[k2];
        if (!Add) {
          r[k2] = 0;
        }
      }
    }
    MXNode::copy_rev(arg[0], r, this->nnz());
    return 0;
  }

  template<bool Add>
  std::string SetNonzerosVector<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << nz_ << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosSlice<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << s_ << "]" << (Add ? " += " : " = ") << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  std::string SetNonzerosSlice2<Add>::disp(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << "(" << arg.at(0) << "[" << outer_ << ";" << inner_ << "]" << (Add ? " += " : " = ")
       << arg.at(1) << ")";
    return ss.str();
  }

  template<bool Add>
  Matrix<casadi_int> SetNonzeros<Add>::mapping() const {
    vector<casadi_int> nz = all();
    return Matrix<casadi_int>(this->dep(1).sparsity(), nz, false);
  }

  template<bool Add>
  bool SetNonzerosVector<Add>::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosVector<Add>* n = dynamic_cast<const SetNonzerosVector<Add>*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->nz_.size()!=n->nz_.size()) return false;
    if (!std::equal(this->nz_.begin(), this->nz_.end(), n->nz_.begin())) return false;

    return true;
  }

  template<bool Add>
  bool SetNonzerosSlice<Add>::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosSlice<Add>* n = dynamic_cast<const SetNonzerosSlice<Add>*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->s_ != n->s_) return false;

    return true;
  }

  template<bool Add>
  bool SetNonzerosSlice2<Add>::is_equal(const MXNode* node, casadi_int depth) const {
    // Check dependencies
    if (!this->sameOpAndDeps(node, depth)) return false;

    // Check if same node
    const SetNonzerosSlice2<Add>* n = dynamic_cast<const SetNonzerosSlice2<Add>*>(node);
    if (n==nullptr) return false;

    // Check sparsity
    if (this->sparsity()!=node->sparsity()) return false;

    // Check indices
    if (this->inner_ != n->inner_ || this->outer_!=n->outer_) return false;

    return true;
  }

  template<bool Add>
  void SetNonzerosVector<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->dep(0).nnz()), this->nnz(),
                          g.work(res[0], this->nnz())) << '\n';
    }

    // Condegen the indices
    std::string ind = g.constant(this->nz_);

    // Perform the operation inplace
    g.local("cii", "const casadi_int", "*");
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g << "for (cii=" << ind << ", rr=" << g.work(res[0], this->nnz()) << ", "
      << "ss=" << g.work(arg[1], this->dep(1).nnz()) << "; cii!=" << ind
      << "+" << this->nz_.size() << "; ++cii, ++ss)"
      << " if (*cii>=0) rr[*cii] " << (Add?"+=":"=") << " *ss;\n";
  }

  template<bool Add>
  void SetNonzerosSlice<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->dep(0).nnz()), this->nnz(),
                          g.work(res[0], this->nnz())) << '\n';
    }

    // Perform the operation inplace
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g << "for (rr=" << g.work(res[0], this->nnz()) << "+" << s_.start << ", ss="
      << g.work(arg[1], this->dep(1).nnz()) << "; rr!="
      << g.work(res[0], this->nnz()) << "+" << s_.stop
      << "; rr+=" << s_.step << ")"
      << " *rr " << (Add?"+=":"=") << " *ss++;\n";
  }

  template<bool Add>
  void SetNonzerosSlice2<Add>::
  generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res) const {
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->dep(0).nnz()), this->nnz(),
                          g.work(res[0], this->nnz())) << '\n';
    }

    // Perform the operation inplace
    g.local("rr", "casadi_real", "*");
    g.local("ss", "casadi_real", "*");
    g.local("tt", "casadi_real", "*");
    g << "for (rr=" << g.work(res[0], this->nnz()) << "+" << outer_.start
      << ", ss=" << g.work(arg[1], this->dep(1).nnz()) << "; rr!="
      << g.work(res[0], this->nnz()) << "+" << outer_.stop
      << "; rr+=" << outer_.step << ")"
      << " for (tt=rr+" << inner_.start << "; tt!=rr+" << inner_.stop
      << "; tt+=" << inner_.step << ")"
      << " *tt " << (Add?"+=":"=") << " *ss++;\n";
  }

} // namespace casadi

/// \endcond

#endif // CASADI_SETNONZEROS_IMPL_HPP
