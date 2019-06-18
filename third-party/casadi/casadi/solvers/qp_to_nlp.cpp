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


#include "qp_to_nlp.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_NLPSOL_EXPORT
  casadi_register_conic_nlpsol(Conic::Plugin* plugin) {
    plugin->creator = QpToNlp::creator;
    plugin->name = "nlpsol";
    plugin->doc = QpToNlp::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &QpToNlp::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_NLPSOL_EXPORT casadi_load_conic_nlpsol() {
    Conic::registerPlugin(casadi_register_conic_nlpsol);
  }

  QpToNlp::QpToNlp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  QpToNlp::~QpToNlp() {
  }

  Options QpToNlp::options_
  = {{&Conic::options_},
     {{"nlpsol",
       {OT_STRING,
        "Name of solver."}},
      {"nlpsol_options",
       {OT_DICT,
        "Options to be passed to solver."}}
     }
  };

  void QpToNlp::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    string nlpsol_plugin;
    Dict nlpsol_options;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="nlpsol") {
        nlpsol_plugin = op.second.to_string();
      } else if (op.first=="nlpsol_options") {
        nlpsol_options = op.second;
      }
    }

    // Create a symbolic matrix for the decision variables
    SX X = SX::sym("X", nx_, 1);

    // Parameters to the problem
    SX H = SX::sym("H", H_);
    SX G = SX::sym("G", nx_);
    SX A = SX::sym("A", A_);

    // Put parameters in a vector
    std::vector<SX> par;
    par.push_back(H.nonzeros());
    par.push_back(G.nonzeros());
    par.push_back(A.nonzeros());

    // The nlp looks exactly like a mathematical description of the NLP
    SXDict nlp = {{"x", X}, {"p", vertcat(par)},
                  {"f", mtimes(G.T(), X) + 0.5*mtimes(mtimes(X.T(), H), X)},
                  {"g", mtimes(A, X)}};

    // Create an Nlpsol instance
    casadi_assert(!nlpsol_plugin.empty(), "'nlpsol' option has not been set");
    solver_ = nlpsol("nlpsol", nlpsol_plugin, nlp, nlpsol_options);
    alloc(solver_);

    // Allocate storage for NLP solver  parameters
    alloc_w(solver_.nnz_in(NLPSOL_P), true);
  }

  int QpToNlp::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // Inputs
    const double *h_, *g_, *a_, *lba_, *uba_, *lbx_, *ubx_, *x0_;
    // Outputs
    double *x_, *f_, *lam_a_, *lam_x_;

    // Get input pointers
    h_ = arg[CONIC_H];
    g_ = arg[CONIC_G];
    a_ = arg[CONIC_A];
    lba_ = arg[CONIC_LBA];
    uba_ = arg[CONIC_UBA];
    lbx_ = arg[CONIC_LBX];
    ubx_ = arg[CONIC_UBX];
    x0_ = arg[CONIC_X0];

    // Get output pointers
    x_ = res[CONIC_X];
    f_ = res[CONIC_COST];
    lam_a_ = res[CONIC_LAM_A];
    lam_x_ = res[CONIC_LAM_X];

    // Buffers for calling the NLP solver
    const double** arg1 = arg + n_in_;
    double** res1 = res + n_out_;
    fill_n(arg1, static_cast<casadi_int>(NLPSOL_NUM_IN), nullptr);
    fill_n(res1, static_cast<casadi_int>(NLPSOL_NUM_OUT), nullptr);

    // NLP inputs
    arg1[NLPSOL_X0] = x0_;
    arg1[NLPSOL_LBG] = lba_;
    arg1[NLPSOL_UBG] = uba_;
    arg1[NLPSOL_LBX] = lbx_;
    arg1[NLPSOL_UBX] = ubx_;

    // NLP parameters
    arg1[NLPSOL_P] = w;

    // Quadratic term
    casadi_int nh = nnz_in(CONIC_H);
    if (h_) {
      copy_n(h_, nh, w);
    } else {
      fill_n(w, nh, 0);
    }
    w += nh;

    // Linear objective term
    casadi_int ng = nnz_in(CONIC_G);
    if (g_) {
      copy_n(g_, ng, w);
    } else {
      fill_n(w, ng, 0);
    }
    w += ng;

    // Linear constraints
    casadi_int na = nnz_in(CONIC_A);
    if (a_) {
      copy_n(a_, na, w);
    } else {
      fill_n(w, na, 0);
    }
    w += na;

    // Solution
    res1[NLPSOL_X] = x_;
    res1[NLPSOL_F] = f_;
    res1[NLPSOL_LAM_X] = lam_x_;
    res1[NLPSOL_LAM_G] = lam_a_;

    // Solve the NLP
    return solver_(arg1, res1, iw, w, 0);
  }

  Dict QpToNlp::get_stats(void* mem) const {
    Dict stats;
    Dict solver_stats = solver_.stats();
    stats["solver_stats"] = solver_stats;
    stats["success"] = solver_stats["success"];

    return stats;
  }

} // namespace casadi
