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


#include "bspline_interpolant.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTERPOLANT_BSPLINE_EXPORT
  casadi_register_interpolant_bspline(Interpolant::Plugin* plugin) {
    plugin->creator = BSplineInterpolant::creator;
    plugin->name = "bspline";
    plugin->doc = BSplineInterpolant::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &BSplineInterpolant::options_;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_BSPLINE_EXPORT casadi_load_interpolant_bspline() {
    Interpolant::registerPlugin(casadi_register_interpolant_bspline);
  }

  Options BSplineInterpolant::options_
  = {{&Interpolant::options_},
      {{"degree",
       {OT_INTVECTOR,
        "Sets, for each grid dimension, the degree of the spline."}},
       {"linear_solver",
        {OT_STRING,
         "Solver used for constructing the coefficient tensor."}},
       {"algorithm",
        {OT_STRING,
         "Algorithm used for fitting the data: 'not_a_knot' (default, same as Matlab),"
        " 'smooth_linear'."}},
       {"smooth_linear_frac",
        {OT_DOUBLE,
         "When 'smooth_linear' algorithm is active, determines sharpness between"
         " 0 (sharp, as linear interpolation) and 0.5 (smooth)."
         "Default value is 0.1."}}
     }
  };

  BSplineInterpolant::
  BSplineInterpolant(const string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const vector<double>& values,
                    casadi_int m)
                    : Interpolant(name, grid, offset, values, m) {

  }

  BSplineInterpolant::~BSplineInterpolant() {
  }

  std::vector<double> meshgrid(const std::vector< std::vector<double> >& grid) {
    std::vector<casadi_int> cnts(grid.size()+1, 0);
    std::vector<casadi_int> sizes(grid.size(), 0);
    for (casadi_int k=0;k<grid.size();++k) sizes[k]= grid[k].size();

    casadi_int total_iter = 1;
    for (casadi_int k=0;k<grid.size();++k) total_iter*= sizes[k];

    casadi_int n_dims = grid.size();

    std::vector<double> ret(total_iter*n_dims);
    for (casadi_int i=0;i<total_iter;++i) {

      for (casadi_int j=0;j<grid.size();++j) {
        ret[i*n_dims+j] = grid[j][cnts[j]];
      }

      cnts[0]++;
      casadi_int j = 0;
      while (j<n_dims && cnts[j]==sizes[j]) {
        cnts[j] = 0;
        j++;
        cnts[j]++;
      }

    }

    return ret;
  }

  std::vector<double> not_a_knot(const std::vector<double>& x, casadi_int k) {
    std::vector<double> ret;
    if (k%2) {
      casadi_int m = (k-1)/2;
      casadi_assert(x.size()>=2*m+2, "Need more data points");
      for (casadi_int i=0;i<k+1;++i) ret.push_back(x[0]);
      for (casadi_int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      for (casadi_int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    } else {
      casadi_error("Not implemented");
      //for (casadi_int i=0;i<k+1;++i) ret.push_back(x[0]);
      //for (casadi_int i=0;i<x.size()-2*m-2;++i) ret.push_back(x[m+1+i]);
      //for (casadi_int i=0;i<k+1;++i) ret.push_back(x[x.size()-1]);
    }
    return ret;
  }

  void BSplineInterpolant::init(const Dict& opts) {
    // Call the base class initializer
    Interpolant::init(opts);

    degree_  = std::vector<casadi_int>(offset_.size()-1, 3);

    linear_solver_ = "lsqr";
    algorithm_ = ALG_NOT_A_KNOT;
    smooth_linear_frac_ = 0.1;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="degree") {
        degree_ = op.second;
      } else if (op.first=="linear_solver") {
        linear_solver_ = op.second.to_string();
      } else if (op.first=="algorithm") {
        std::string alg = op.second.to_string();
        if (alg=="not_a_knot") {
          algorithm_ = ALG_NOT_A_KNOT;
        } else if (alg=="smooth_linear") {
          algorithm_ = ALG_SMOOTH_LINEAR;
        } else {
          casadi_error("Algorithm option invalid: " + get_options().info("algorithm"));
        }
      } else if (op.first=="smooth_linear_frac") {
        smooth_linear_frac_ = op.second;
        casadi_assert(smooth_linear_frac_>0 && smooth_linear_frac_<0.5,
          "smooth_linear_frac must be in ]0,0.5[");
      }
    }

    casadi_assert_dev(degree_.size()==offset_.size()-1);

    std::vector< std::vector<double> > grid;
    for (casadi_int k=0;k<degree_.size();++k) {
      std::vector<double> local_grid(grid_.begin()+offset_[k], grid_.begin()+offset_[k+1]);
      grid.push_back(local_grid);
    }

    Dict opts_bspline;
    opts_bspline["lookup_mode"] = lookup_modes_;

    switch (algorithm_) {
      case ALG_NOT_A_KNOT:
        {
          std::vector< std::vector<double> > knots;
          for (casadi_int k=0;k<degree_.size();++k)
            knots.push_back(not_a_knot(grid[k], degree_[k]));
          Dict opts_dual;
          opts_dual["ad_weight_sp"] = 0;
          opts_dual["lookup_mode"] = lookup_modes_;

          Function B = Function::bspline_dual("spline", knots, meshgrid(grid), degree_, 1, false,
            opts_dual);

          Function Jf = B.jacobian_old(0, 0);

          MX C = MX::sym("C", B.size_in(0));

          MX Js = Jf(std::vector<MX>{C})[0];
          Function temp = Function("J", {C}, {Js});
          DM J = temp(std::vector<DM>{0})[0];

          casadi_assert_dev(J.size1()==J.size2());

          DM V = DM::reshape(DM(values_), m_, -1).T();
          DM C_opt = solve(J, V, linear_solver_);

          double fit = static_cast<double>(norm_1(mtimes(J, C_opt) - V));

          if (verbose_) casadi_message("Lookup table fitting error: " + str(fit));

          S_ = Function::bspline("spline", knots, C_opt.T().nonzeros(), degree_, m_, opts_bspline);
        }
        break;
      case ALG_SMOOTH_LINEAR:
        {
          casadi_int n_dim = degree_.size();
          // Linear fit
          Function linear = interpolant("linear", "linear", grid, values_);

          std::vector< std::vector<double> > egrid;
          std::vector< std::vector<double> > new_grid;

          for (casadi_int k=0;k<n_dim;++k) {
            casadi_assert(degree_[k]==3, "Only degree 3 supported for 'smooth_linear'.");

            // Add extra knots
            const std::vector<double>& g = grid[k];

            // Determine smallest gap.
            double m = inf;
            for (casadi_int i=0;i<g.size()-1;++i) {
              double delta = g[i+1]-g[i];
              if (delta<m) m = delta;
            }
            double step = smooth_linear_frac_*m;

            // Add extra knots
            std::vector<double> new_g;
            new_g.push_back(g.front());
            new_g.push_back(g.front()+step);
            for (casadi_int i=1;i<g.size()-1;++i) {
              new_g.push_back(g[i]-step);
              new_g.push_back(g[i]);
              new_g.push_back(g[i]+step);
            }
            new_g.push_back(g.back()-step);
            new_g.push_back(g.back());
            new_grid.push_back(new_g);

            // Correct multiplicity
            double v1 = new_g.front();
            double vend = new_g.back();
            new_g.insert(new_g.begin(), degree_[k], v1);
            new_g.insert(new_g.end(), degree_[k], vend);

            grid[k] = new_g;

            // Compute greville points
            egrid.push_back(greville_points(new_g, degree_[k]));
          }

          // Evaluate linear interpolation on greville grid
          std::vector<double> mg = meshgrid(egrid);
          std::vector<double> Z(m_*mg.size()/n_dim);

          // Work vectors
          vector<casadi_int> iw(linear.sz_iw());
          vector<double> w(linear.sz_w());

          std::vector<const double*> arg(1);
          std::vector<double*> res(1);

          for (int i=0;i<Z.size();++i) {
            arg[0] = get_ptr(mg)+n_dim*i;
            res[0] = get_ptr(Z)+m_*i;
            linear(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
          }
          S_ = Function::bspline("spline", grid, Z, degree_, m_, opts_bspline);
        }
        break;
    }

    alloc_w(S_->sz_w(), true);
    alloc_iw(S_->sz_iw(), true);

  }



  std::vector<double> BSplineInterpolant::greville_points(const std::vector<double>& x,
                                                          casadi_int deg) {
    casadi_int dim = x.size()-deg-1;
    std::vector<double> ret(dim);
    for (casadi_int i = 0; i < dim; ++i) {
      ret[i] = 0;
      for (casadi_int j = 0; j < deg; j++) {
        ret[i] += x[i+1+j];
      }
      ret[i] = ret[i] / deg;
    }
    return ret;
  }

  int BSplineInterpolant::eval(const double** arg, double** res,
                                casadi_int* iw, double* w, void* mem) const {
    return S_->eval(arg, res, iw, w, mem);
  }

  void BSplineInterpolant::codegen_body(CodeGenerator& g) const {
    S_->codegen_body(g);
  }

  Function BSplineInterpolant::
  get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    return S_->get_jacobian(name, inames, onames, opts);
  }

} // namespace casadi
