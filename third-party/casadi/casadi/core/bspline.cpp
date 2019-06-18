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


#include "bspline.hpp"
#include "function_internal.hpp"
#include "interpolant_impl.hpp"
#include "casadi_misc.hpp"
#include "mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  Options BSplineCommon::options_
  = {{&FunctionInternal::options_},
      {{"lookup_mode",
        {OT_STRINGVECTOR,
         "Specifies, for each grid dimenion, the lookup algorithm used to find the correct index. "
         "'linear' uses a for-loop + break; (default when #knots<=100), "
         "'exact' uses floored division (only for uniform grids), "
         "'binary' uses a binary search. (default when #knots>100)."}}
      }
  };


  void BSplineCommon::init(const Dict& opts) {
    casadi::FunctionInternal::init(opts);

    std::vector<std::string> lookup_mode;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_mode = op.second;
      }
    }

    lookup_mode_ = Interpolant::interpret_lookup_mode(lookup_mode, knots_,
      offset_, degree_, degree_);

    casadi_int n_dims = degree_.size();

    casadi_int k;
    for (k=0;k<n_dims-1;++k) {
      alloc_w(degree_[k]+1, true); // boor
    }
    alloc_w(2*degree_[k]+1, true);
    alloc_w(n_dims+1, true);
    alloc_iw(n_dims+1, true); // boor_offset
    alloc_iw(n_dims, true); // starts
    alloc_iw(n_dims, true); // index
    alloc_iw(n_dims+1, true); // coeff_offset

    coeffs_dims_.resize(n_dims+1);
    coeffs_dims_[0] = m_;
    for (casadi_int i=0;i<n_dims;++i) coeffs_dims_[i+1] = offset_[i+1]-offset_[i]-degree_[i]-1;

    // Prepare strides
    strides_.resize(n_dims);
    strides_[0] = m_;
    for (casadi_int i=0;i<n_dims-1;++i) {
      strides_[i+1] = strides_[i]*coeffs_dims_[i+1];
    }

  }


  void BSplineCommon::from_knots(const std::vector< std::vector<double> >& knots,
    std::vector<casadi_int>& offset, std::vector<double>& stacked) {

    // Get offset for each input dimension
    offset.clear();
    offset.reserve(knots.size()+1);
    offset.push_back(0);
    for (auto&& g : knots) offset.push_back(offset.back()+g.size());

    // Stack input grids
    stacked.clear();
    stacked.reserve(offset.back());
    for (auto&& g : knots) stacked.insert(stacked.end(), g.begin(), g.end());

  }

  Function BSpline::create(const std::string &name,
    const std::vector< std::vector<double> >& knots,
    const vector<double>& coeffs, const vector<casadi_int>& degree,
    casadi_int m, const Dict& opts) {

      vector<casadi_int> offset;
      vector<double> stacked;
      BSplineCommon::from_knots(knots, offset, stacked);
      return Function::create(new BSpline(name, stacked, offset, coeffs, degree, m), opts);
    }

    BSplineCommon::BSplineCommon(const std::string &name, const std::vector<double>& knots,
      const std::vector<casadi_int>& offset, const vector<casadi_int>& degree, casadi_int m) :
        casadi::FunctionInternal(name), knots_(knots), offset_(offset), degree_(degree), m_(m) {

      coeffs_size_= m;
      for (casadi_int i=0;i<degree_.size();++i) coeffs_size_*= offset_[i+1]-offset_[i]-degree_[i]-1;
    }

    BSpline::BSpline(const std::string &name, const std::vector<double>& knots,
        const std::vector<casadi_int>& offset, const vector<double>& coeffs,
        const vector<casadi_int>& degree,
        casadi_int m) :
        casadi::BSplineCommon(name, knots, offset, degree, m), coeffs_(coeffs) {}

    void BSpline::init(const Dict& opts) {
      casadi::BSplineCommon::init(opts);

      casadi_assert(coeffs_size_==coeffs_.size(),
        "Expected coefficient size " + str(coeffs_size_) + ", "
        "got " + str(coeffs_.size()) + " instead.");
    }

    int BSpline::eval(const double** arg, double** res,
                      casadi_int* iw, double* w, void* mem) const {
      if (!res[0]) return 0;

      casadi_fill(res[0], m_, 0.0);
      casadi_nd_boor_eval(res[0], degree_.size(), get_ptr(knots_), get_ptr(offset_),
        get_ptr(degree_), get_ptr(strides_), get_ptr(coeffs_), m_, arg[0], get_ptr(lookup_mode_),
        false, iw, w);
      return 0;
    }

    void BSpline::codegen_body(CodeGenerator& g) const {
      casadi_int n_dims = offset_.size()-1;

      g.add_auxiliary(CodeGenerator::AUX_ND_BOOR_EVAL);
      g.add_auxiliary(CodeGenerator::AUX_FILL);
      g << "  if (res[0]) " << g.fill("res[0]", m_, "0.0") << "\n";

      // Input and output buffers
      g << "  if (res[0]) CASADI_PREFIX(nd_boor_eval)(res[0]," << n_dims << ","
        << g.constant(knots_) << "," << g.constant(offset_) << "," <<  g.constant(degree_)
        << "," << g.constant(strides_) << "," <<  g.constant(coeffs_) << "," << m_  << ",arg[0],"
        <<  g.constant(lookup_mode_) << ", 0, iw, w);\n";
    }

    Function BSpline::get_forward(casadi_int nfwd, const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
      MX x = MX::sym("x", degree_.size());
      MX J = jac(x);

      std::vector<MX> seed = MX::sym("seed", degree_.size(), 1, nfwd);
      std::vector<MX> sens;
      for (casadi_int i=0;i<nfwd;++i) {
        sens.push_back(mtimes(J, seed[i]));
      }

      MX dummy = MX(m_, 1);
      return Function(name, {x, dummy, horzcat(seed)}, {horzcat(sens)}, opts);
    }

    Function BSpline::get_reverse(casadi_int nadj, const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {

      MX x = MX::sym("x", degree_.size());
      MX JT = jac(x).T();

      std::vector<MX> seed = MX::sym("seed", m_, 1, nadj);
      std::vector<MX> sens;
      for (casadi_int i=0;i<nadj;++i) {
        sens.push_back(mtimes(JT, seed[i]));
      }

      MX dummy = MX(m_, 1);
      return Function(name, {x, dummy, horzcat(seed)}, {horzcat(sens)}, opts);

    }

    MX BSpline::jac(const MX& x) const {

      casadi_int n_dims = degree_.size();
      std::vector<MX> parts;

      // Loop over dimensions
      for (casadi_int k=0;k<n_dims;++k) {
        std::vector< std::vector<double> > knots;
        std::vector< casadi_int> degree;
        for (casadi_int i=0;i<degree_.size();++i) {
          if (i==k) {
            knots.push_back(
              std::vector<double>(get_ptr(knots_)+offset_[i]+1, get_ptr(knots_)+offset_[i+1]-1));
            degree.push_back(degree_[i]-1);
          } else {
            knots.push_back(
              std::vector<double>(get_ptr(knots_)+offset_[i], get_ptr(knots_)+offset_[i+1]));
            degree.push_back(degree_[i]);
          }
        }
        Dict opts;
        opts["lookup_mode"] = Interpolant::lookup_mode_from_enum(lookup_mode_);
        Function d = Function::bspline("jac_helper", knots, derivative_coeff(k), degree, m_, opts);
        parts.push_back(d(std::vector<MX>{x})[0]);
      }

      return horzcat(parts);
    }

    std::vector<double> BSpline::derivative_coeff(casadi_int i) const {
      casadi_int n_dims = degree_.size();

      std::vector<double> coeffs;
      casadi_int n_knots = offset_[i+1]-offset_[i];
      casadi_int n = n_knots-degree_[i]-1;
      DM knots = std::vector<double>(get_ptr(knots_)+offset_[i], get_ptr(knots_)+offset_[i+1]);
      DM delta_knots = knots(range(1+degree_[i], n_knots-1))
           - knots(range(1, n_knots-degree_[i]-1));
      Sparsity sp_diag = vertsplit(Sparsity::diag(n), {0, n-1, n})[0];
      Sparsity sp_band = vertsplit(Sparsity::band(n, -1), {0, n-1, n})[0];
      DM delta_knots_inv = 1/delta_knots;
      DM T = DM(sp_diag, -delta_knots_inv) + DM(sp_band, delta_knots_inv);
      T*= degree_[i];

      std::vector<casadi_int> coeffs_dims_new = coeffs_dims_;
      coeffs_dims_new[i+1] = T.size1();

      // Apply transformation T on axis i

      // Bring axis i to the back
      std::vector<casadi_int> order = range(n_dims+1);
      std::swap(order.back(), order[i+1]);
      std::vector<casadi_int> mapping = tensor_permute_mapping(coeffs_dims_, order);
      DM coeff_matrix = DM(coeffs_).nz(mapping);

      // Cast as matrix
      coeff_matrix = reshape(coeff_matrix, -1, T.size2());

      // Apply the transformation matrix from the right
      coeff_matrix = mtimes(coeff_matrix, T.T());

      // Bring axis i back to the original place
      mapping = tensor_permute_mapping(permute(coeffs_dims_new, order), order);
      coeff_matrix = coeff_matrix.nz(mapping);

      // Return the flat vector
      return coeff_matrix.nonzeros();
    }


    Function BSpline::get_jacobian(const std::string& name,
          const std::vector<std::string>& inames,
          const std::vector<std::string>& onames, const Dict& opts) const {
      MX x = MX::sym(inames.at(0), degree_.size());
      MX dummy = MX::sym(inames.at(1), Sparsity(m_, 1));
      return Function(name, {x, dummy}, {jac(x)}, opts);
    }

    Sparsity BSplineDual::get_sparsity_in(casadi_int i) {
      if (reverse_) return Sparsity::dense(m_, N_);
      return Sparsity::dense(coeffs_size_);
    }
    Sparsity BSplineDual::get_sparsity_out(casadi_int i) {
      if (reverse_) return Sparsity::dense(coeffs_size_);
      return Sparsity::dense(m_, N_);
    }

    Function BSplineDual::create(const std::string &name,
        const std::vector< std::vector<double> >& knots,
        const vector<double>& x, const vector<casadi_int>& degree,
        casadi_int m, bool reverse, const Dict& opts) {

      vector<casadi_int> offset;
      vector<double> stacked;
      BSplineCommon::from_knots(knots, offset, stacked);
      return Function::create(new BSplineDual(name, stacked, offset, x, degree, m, reverse), opts);
    }

    BSplineDual::BSplineDual(const std::string &name, const std::vector<double>& knots,
        const std::vector<casadi_int>& offset, const vector<double>& x,
        const vector<casadi_int>& degree, casadi_int m,
        bool reverse) :
        casadi::BSplineCommon(name, knots, offset, degree, m), x_(x), reverse_(reverse) {

      casadi_int n_dims = degree_.size();
      N_ = x_.size()/n_dims;
      casadi_assert_dev(N_*n_dims==x_.size());
    }

    void BSplineDual::init(const Dict& opts) {
      BSplineCommon::init(opts);
    }

    int BSplineDual::eval(const double** arg, double** res,
                          casadi_int* iw, double* w, void* mem) const {
      if (!res[0]) return 0;
      casadi_fill(res[0], reverse_? coeffs_size_: m_*N_, 0.0);

      casadi_int n_dims = degree_.size();

      for (casadi_int i=0;i<N_;++i) {
        casadi_nd_boor_eval(res[0]+(reverse_? 0 : i*m_), n_dims, get_ptr(knots_), get_ptr(offset_),
        get_ptr(degree_), get_ptr(strides_), arg[0]+(reverse_? i*m_ : 0), m_, get_ptr(x_)+i*n_dims,
        get_ptr(lookup_mode_), reverse_, iw, w);
      }
      return 0;
    }

    void BSplineDual::codegen_body(CodeGenerator& g) const {
      casadi_int n_dims = offset_.size()-1;

      g.add_auxiliary(CodeGenerator::AUX_ND_BOOR_EVAL);
      g.add_auxiliary(CodeGenerator::AUX_FILL);
      g << "  if (res[0]) " <<
                g.fill("res[0]", reverse_? coeffs_size_: m_*N_, "0.0") << "\n";

      // Input and output buffers
      g << "  if (res[0]) ";
      g << "for (casadi_int i=0;i<" << N_ << ";++i) CASADI_PREFIX(nd_boor_eval)(res[0]"
        << (reverse_? "" : "+i*" + str(m_)) << "," << n_dims << "," << g.constant(knots_)
        << "," << g.constant(offset_) << "," <<  g.constant(degree_)
        << "," << g.constant(strides_) << ",arg[0]" << (reverse_? "i*" + str(m_) : "")
        << "," << m_  << "," << g.constant(x_) <<"+i*" << n_dims << "," <<  g.constant(lookup_mode_)
        << ", 0, iw, w);\n";
    }

    Function BSplineDual::get_forward(casadi_int nfwd, const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {

      MX C = MX::sym("C", sparsity_in_.at(0));
      MX dummy = MX(size_out(0));
      std::vector<MX> seed = MX::sym("seed", sparsity_in_.at(0), nfwd);
      std::vector<MX> sens;
      Function self = shared_from_this<Function>();
      for (casadi_int i=0;i<nfwd;++i) {
        sens.push_back(self(seed[i])[0]);
      }
      Function ret = Function(name, {C, dummy, horzcat(seed)}, {horzcat(sens)}, inames, onames);
      return ret;

    }

    Function BSplineDual::get_reverse(casadi_int nadj, const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {

      MX C = MX::sym("C", sparsity_in_.at(0));
      MX dummy = MX(size_out(0));
      std::vector<MX> seed = MX::sym("seed", sparsity_out_.at(0), nadj);
      std::vector<MX> sens;
      Function rev = Function::create(new BSplineDual(name, knots_, offset_, x_,
                                                      degree_, m_, !reverse_), opts);
      for (casadi_int i=0;i<nadj;++i) {
        sens.push_back(rev(seed[i])[0]);
      }
      Function ret = Function(name, {C, dummy, horzcat(seed)}, {horzcat(sens)}, inames, onames);
      return ret;

    }

    void nd_boor_eval_sp(bvec_t* ret, casadi_int n_dims, const double* all_knots,
        const casadi_int* offset, const casadi_int* all_degree, const casadi_int* strides,
        const bvec_t* c, casadi_int m, const double* all_x,
        const casadi_int* lookup_mode, casadi_int reverse, casadi_int* iw, bvec_t* w) {
      casadi_int* boor_offset = iw; iw+=n_dims+1;
      casadi_int* starts = iw; iw+=n_dims;
      casadi_int* index = iw; iw+=n_dims;
      casadi_int* coeff_offset = iw;

      boor_offset[0] = 0;
      coeff_offset[n_dims] = 0;

      casadi_int n_iter = 1;
      for (casadi_int k=0;k<n_dims;++k) {
        casadi_int degree = all_degree[k];
        const double* knots = all_knots + offset[k];
        casadi_int n_knots = offset[k+1]-offset[k];
        casadi_int n_b = n_knots-degree-1;

        double x = all_x[k];
        casadi_int L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[k]);

        casadi_int start = L;
        if (start>n_b-degree-1) start = n_b-degree-1;

        starts[k] = start;
        boor_offset[k+1] = boor_offset[k] + degree+1;
        n_iter*= degree+1;
      }

      casadi_fill(index, n_dims, casadi_int(0));

      // Prepare offset
      for (casadi_int pivot=n_dims-1;pivot>=0;--pivot) {
        coeff_offset[pivot] = starts[pivot]*strides[pivot]+coeff_offset[pivot+1];
      }

      for (casadi_int k=0;k<n_iter;++k) {

        // accumulate result
        for (casadi_int i=0;i<m;++i) {
          if (reverse) {
            ret[coeff_offset[0]+i] |= c[i];
          } else {
            ret[i] |= c[coeff_offset[0]+i];
          }
        }

        // Increment index
        index[0]++;
        casadi_int pivot = 0;

        // Handle index overflow
        {
          // increment next index (forward)
          while (index[pivot]==boor_offset[pivot+1]-boor_offset[pivot]) {
            index[pivot] = 0;
            if (pivot==n_dims-1) break;
            index[++pivot]++;
          }

          // update cumulative structures (reverse)
          while (pivot>0) {
            // Compute offset
            coeff_offset[pivot] = (starts[pivot]+index[pivot])*strides[pivot]+coeff_offset[pivot+1];
            pivot--;
          }
        }

        // Compute offset
        coeff_offset[0] = (starts[0]+index[0])*m+coeff_offset[1];

      }
    }

    int BSplineDual::
    sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
      if (!res[0]) return 0;
      casadi_fill(res[0], reverse_? coeffs_size_: m_*N_, bvec_t(0));

      casadi_int n_dims = degree_.size();
      for (casadi_int i=0;i<N_;++i) {
        nd_boor_eval_sp(res[0]+(reverse_? 0 : i*m_), n_dims, get_ptr(knots_), get_ptr(offset_),
          get_ptr(degree_), get_ptr(strides_), arg[0]+(reverse_? i*m_ : 0), m_,
          get_ptr(x_)+i*n_dims, get_ptr(lookup_mode_), reverse_, iw, w);
      }
      return 0;
    }

    int BSplineDual::
    sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
      if (!res[0]) return 0;
      casadi_int n_dims = degree_.size();
      for (casadi_int i=0;i<N_;++i) {
        nd_boor_eval_sp(arg[0]+(!reverse_? 0 : i*m_), n_dims, get_ptr(knots_), get_ptr(offset_),
          get_ptr(degree_), get_ptr(strides_), res[0]+(!reverse_? i*m_ : 0), m_,
          get_ptr(x_)+i*n_dims, get_ptr(lookup_mode_), !reverse_, iw, w);
      }
      return 0;
    }

} // namespace casadi
