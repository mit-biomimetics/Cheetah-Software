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


#include "rootfinder_impl.hpp"
#include "mx_node.hpp"
#include <iterator>
#include "linsol.hpp"

#include "global_options.hpp"

using namespace std;
namespace casadi {

  vector<string> rootfinder_in() {
    vector<string> ret(rootfinder_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=rootfinder_in(i);
    return ret;
  }

  vector<string> rootfinder_out() {
    vector<string> ret(rootfinder_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=rootfinder_out(i);
    return ret;
  }

  string rootfinder_in(casadi_int ind) {
    switch (static_cast<RootfinderInput>(ind)) {
    case ROOTFINDER_X0:  return "x0";
    case ROOTFINDER_P:   return "p";
    case ROOTFINDER_NUM_IN: break;
    }
    return string();
  }

  string rootfinder_out(casadi_int ind) {
    switch (static_cast<RootfinderOutput>(ind)) {
    case ROOTFINDER_X:  return "x";
    case ROOTFINDER_NUM_OUT: break;
    }
    return string();
  }

  casadi_int rootfinder_n_in() {
    return ROOTFINDER_NUM_IN;
  }

  casadi_int rootfinder_n_out() {
    return ROOTFINDER_NUM_OUT;
  }

  std::vector<std::string> rootfinder_options(const std::string& name) {
    return Rootfinder::plugin_options(name).all();
  }

  std::string rootfinder_option_type(const std::string& name, const std::string& op) {
    return Rootfinder::plugin_options(name).type(op);
  }

  std::string rootfinder_option_info(const std::string& name, const std::string& op) {
    return Rootfinder::plugin_options(name).info(op);
  }

  bool has_rootfinder(const string& name) {
    return Rootfinder::has_plugin(name);
  }

  void load_rootfinder(const string& name) {
    Rootfinder::load_plugin(name);
  }

  string doc_rootfinder(const string& name) {
    return Rootfinder::getPlugin(name).doc;
  }

  Function rootfinder(const string& name, const string& solver,
                      const SXDict& rfp, const Dict& opts) {
    return rootfinder(name, solver, Rootfinder::create_oracle(rfp, opts), opts);
  }

  Function rootfinder(const string& name, const string& solver,
                      const MXDict& rfp, const Dict& opts) {
    return rootfinder(name, solver, Rootfinder::create_oracle(rfp, opts), opts);
  }

  template<typename XType>
  Function Rootfinder::create_oracle(const std::map<std::string, XType>& d,
                                 const Dict& opts) {
    std::vector<XType> rfp_in(RFP_NUM_IN), rfp_out(RFP_NUM_OUT);
    for (auto&& i : d) {
      if (i.first=="x") {
        rfp_in[RFP_X]=i.second;
      } else if (i.first=="p") {
        rfp_in[RFP_P]=i.second;
      } else if (i.first=="g") {
        rfp_out[RFP_G]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }

    // Options for the oracle
    Dict oracle_options;
    Dict::const_iterator it = opts.find("oracle_options");
    if (it!=opts.end()) {
      // "oracle_options" has been set
      oracle_options = it->second;
    } else if ((it=opts.find("verbose")) != opts.end()) {
      // "oracle_options" has not been set, but "verbose" has
      oracle_options["verbose"] = it->second;
    }

    // Create oracle
    return Function("rfp", rfp_in, rfp_out, {"x0", "p"}, {"x"}, oracle_options);
  }

  Function rootfinder(const std::string& name, const std::string& solver,
                   const Function& f, const Dict& opts) {
    // Make sure that residual function is sound
    if (f.has_free()) {
      casadi_error("Cannot create '" + name + "' since " + str(f.get_free()) + " are free.");
    }
    return Function::create(Rootfinder::instantiate(name, solver, f), opts);
  }

  Rootfinder::Rootfinder(const std::string& name, const Function& oracle)
    : OracleFunction(name, oracle) {

    // Default options
    iin_ = 0;
    iout_ = 0;
    // TODO(jgillis): remove hack in new release -- need uniform default.
    error_on_fail_ = name=="kinsol" ? true : false;
  }

  Rootfinder::~Rootfinder() {
  }

  Options Rootfinder::options_
  = {{&OracleFunction::options_},
     {{"linear_solver",
       {OT_STRING,
        "User-defined linear solver class. Needed for sensitivities."}},
      {"linear_solver_options",
       {OT_DICT,
        "Options to be passed to the linear solver."}},
      {"constraints",
       {OT_INTVECTOR,
        "Constrain the unknowns. 0 (default): no constraint on ui, "
        "1: ui >= 0.0, -1: ui <= 0.0, 2: ui > 0.0, -2: ui < 0.0."}},
      {"implicit_input",
       {OT_INT,
        "Index of the input that corresponds to the actual root-finding"}},
      {"implicit_output",
       {OT_INT,
        "Index of the output that corresponds to the actual root-finding"}},
      {"jacobian_function",
       {OT_FUNCTION,
        "Function object for calculating the Jacobian (autogenerated by default)"}},
      {"error_on_fail",
       {OT_BOOL,
        "When the numerical process returns unsuccessfully, raise an error (default false)."}}
     }
  };

  void Rootfinder::init(const Dict& opts) {

    // Default (temporary) options
    Dict linear_solver_options;
    string linear_solver = "qr";
    Function jac; // Jacobian of f with respect to z

    // Read options
    for (auto&& op : opts) {
      if (op.first=="implicit_input") {
        iin_ = op.second;
      } else if (op.first=="implicit_output") {
        iout_ = op.second;
      } else if (op.first=="jacobian_function") {
        jac = op.second;
      } else if (op.first=="linear_solver_options") {
        linear_solver_options = op.second;
      } else if (op.first=="linear_solver") {
        linear_solver = op.second.to_string();
      } else if (op.first=="constraints") {
        u_c_ = op.second;
      } else if (op.first=="error_on_fail") {
        error_on_fail_ = op.second;
      }
    }

    // Get the number of equations and check consistency
    casadi_assert(iin_>=0 && iin_<oracle_.n_in() && oracle_.n_in()>0,
                          "Implicit input not in range");
    casadi_assert(iout_>=0 && iout_<oracle_.n_out() && oracle_.n_out()>0,
                          "Implicit output not in range");
    casadi_assert(oracle_.sparsity_out(iout_).is_dense()
                          && oracle_.sparsity_out(iout_).is_column(),
                          "Residual must be a dense vector");
    casadi_assert(oracle_.sparsity_in(iin_).is_dense()
                          && oracle_.sparsity_in(iin_).is_column(),
                          "Unknown must be a dense vector");
    n_ = oracle_.nnz_out(iout_);
    casadi_assert(n_ == oracle_.nnz_in(iin_),
      "Dimension mismatch. Input size is " + str(oracle_.nnz_in(iin_)) + ", "
      "while output size is " + str(oracle_.nnz_out(iout_)));

    // Call the base class initializer
    OracleFunction::init(opts);

    // Generate Jacobian if not provided
    if (jac.is_null()) jac = oracle_.jacobian_old(iin_, iout_);
    set_function(jac, "jac_f_z");
    sp_jac_ = jac.sparsity_out(0);

    // Check for structural singularity in the Jacobian
    casadi_assert(!sp_jac_.is_singular(),
      "Rootfinder::init: singularity - the jacobian is structurally rank-deficient. "
      "sprank(J)=" + str(sprank(sp_jac_)) + " (instead of " + str(sp_jac_.size1()) + ")");

    // Get the linear solver creator function
    linsol_ = Linsol("linsol", linear_solver, sp_jac_, linear_solver_options);

    // Constraints
    casadi_assert(u_c_.size()==n_ || u_c_.empty(),
      "Constraint vector if supplied, must be of length n, but got "
      + str(u_c_.size()) + " and n = " + str(n_));

    // Allocate sufficiently large work vectors
    alloc(oracle_);
    size_t sz_w = oracle_.sz_w();
    if (!jac.is_null()) {
      sz_w = max(sz_w, jac.sz_w());
    }
    alloc_w(sz_w + 2*static_cast<size_t>(n_));
  }

  int Rootfinder::init_mem(void* mem) const {
    if (OracleFunction::init_mem(mem)) return 1;
    //auto m = static_cast<RootfinderMemory*>(mem);
    return 0;
  }

  int Rootfinder::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    // Reset the solver, prepare for solution
    setup(mem, arg, res, iw, w);

    // Solve the NLP
    int ret = solve(mem);
    auto m = static_cast<RootfinderMemory*>(mem);
    if (error_on_fail_ && !m->success)
      casadi_error("rootfinder process failed. "
                   "Set 'error_on_fail' option to false to ignore this error.");

    return ret;
  }

  void Rootfinder::set_work(void* mem, const double**& arg, double**& res,
                        casadi_int*& iw, double*& w) const {
    auto m = static_cast<RootfinderMemory*>(mem);

    // Problem has not been solved at this point
    m->success = false;

    // Get input pointers
    m->iarg = arg;
    arg += n_in_;

    // Get output pointers
    m->ires = res;
    res += n_out_;
  }

  Function Rootfinder
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Symbolic expression for the input
    vector<MX> arg = mx_in(), res = mx_out();
    vector<vector<MX>> fseed = fwd_seed<MX>(nfwd), fsens;
    arg[iin_] = MX::sym(arg[iin_].name(), Sparsity(arg[iin_].size()));
    for (auto&& e : fseed) e[iin_] = MX::sym(e[iin_].name(), e[iin_].size());
    ad_forward(arg, res, fseed, fsens, false, false);

    // Construct return function
    arg.insert(arg.end(), res.begin(), res.end());
    vector<MX> v(nfwd);
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = fseed[d][i];
      arg.push_back(horzcat(v));
    }
    res.clear();
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = fsens[d][i];
      res.push_back(horzcat(v));
    }
    return Function(name, arg, res, inames, onames, opts);
  }

  Function Rootfinder
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Symbolic expression for the input
    vector<MX> arg = mx_in();
    arg[iin_] = MX::sym(arg[iin_].name() + "_guess",
                        Sparsity(arg[iin_].size()));
    vector<MX> res = mx_out();
    vector<vector<MX> > aseed = symbolicAdjSeed(nadj, res), asens;
    ad_reverse(arg, res, aseed, asens, false, false);

    // Construct return function
    arg.insert(arg.end(), res.begin(), res.end());
    vector<MX> v(nadj);
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = aseed[d][i];
      arg.push_back(horzcat(v));
    }
    res.clear();
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = asens[d][i];
      res.push_back(horzcat(v));
    }
    return Function(name, arg, res, inames, onames, opts);
  }

  int Rootfinder::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    bvec_t* tmp1 = w; w += n_;
    bvec_t* tmp2 = w; w += n_;

    // Propagate dependencies through the function
    const bvec_t** arg1 = arg+n_in_;
    copy(arg, arg+n_in_, arg1);
    arg1[iin_] = nullptr;
    bvec_t** res1 = res+n_out_;
    fill_n(res1, n_out_, static_cast<bvec_t*>(nullptr));
    res1[iout_] = tmp1;
    oracle_(arg1, res1, iw, w, 0);

    // "Solve" in order to propagate to z
    fill_n(tmp2, n_, 0);
    sp_jac_.spsolve(tmp2, tmp1, false);
    if (res[iout_]) copy(tmp2, tmp2+n_, res[iout_]);

    // Propagate to auxiliary outputs
    if (n_out_>1) {
      arg1[iin_] = tmp2;
      copy(res, res+n_out_, res1);
      res1[iout_] = nullptr;
      oracle_(arg1, res1, iw, w, 0);
    }
    return 0;
  }

  int Rootfinder::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    bvec_t* tmp1 = w; w += n_;
    bvec_t* tmp2 = w; w += n_;

    // Get & clear seed corresponding to implicitly defined variable
    if (res[iout_]) {
      copy(res[iout_], res[iout_]+n_, tmp1);
      fill_n(res[iout_], n_, 0);
    } else {
      fill_n(tmp1, n_, 0);
    }

    // Propagate dependencies from auxiliary outputs to z
    bvec_t** res1 = res+n_out_;
    copy(res, res+n_out_, res1);
    res1[iout_] = nullptr;
    bvec_t** arg1 = arg+n_in_;
    copy(arg, arg+n_in_, arg1);
    arg1[iin_] = tmp1;
    if (n_out_>1) {
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    }

    // "Solve" in order to get seed
    fill_n(tmp2, n_, 0);
    sp_jac_.spsolve(tmp2, tmp1, true);

    // Propagate dependencies through the function
    for (casadi_int i=0; i<n_out_; ++i) res1[i] = nullptr;
    res1[iout_] = tmp2;
    arg1[iin_] = nullptr; // just a guess
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    return 0;
  }

  std::map<std::string, Rootfinder::Plugin> Rootfinder::solvers_;

  const std::string Rootfinder::infix_ = "rootfinder";

  void Rootfinder::
  ad_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
          const std::vector<std::vector<MX> >& fseed,
          std::vector<std::vector<MX> >& fsens,
          bool always_inline, bool never_inline) const {
    // Number of directional derivatives
    casadi_int nfwd = fseed.size();
    fsens.resize(nfwd);

    // Quick return if no seeds
    if (nfwd==0) return;

    // Propagate through f_
    vector<MX> f_arg(arg);
    f_arg.at(iin_) = res.at(iout_);
    vector<MX> f_res(res);
    f_res.at(iout_) = MX(size_in(iin_)); // zero residual
    std::vector<std::vector<MX> > f_fseed(fseed);
    for (casadi_int d=0; d<nfwd; ++d) {
      f_fseed[d].at(iin_) = MX(size_in(iin_)); // ignore seeds for guess
    }
    oracle_->call_forward(f_arg, f_res, f_fseed, fsens,
                          always_inline, never_inline);

    // Get expression of Jacobian
    Function jac = get_function("jac_f_z");
    MX J = jac(f_arg).front();

    // Solve for all the forward derivatives at once
    vector<MX> rhs(nfwd);
    for (casadi_int d=0; d<nfwd; ++d) rhs[d] = vec(fsens[d][iout_]);
    rhs = horzsplit(J->get_solve(-horzcat(rhs), false, linsol_));
    for (casadi_int d=0; d<nfwd; ++d) fsens[d][iout_] = reshape(rhs[d], size_in(iin_));

    // Propagate to auxiliary outputs
    if (n_out_>1) {
      for (casadi_int d=0; d<nfwd; ++d) f_fseed[d][iin_] = fsens[d][iout_];
      oracle_->call_forward(f_arg, f_res, f_fseed, fsens,
                            always_inline, never_inline);
      for (casadi_int d=0; d<nfwd; ++d) fsens[d][iout_] = f_fseed[d][iin_]; // Otherwise overwritten
    }
  }

  void Rootfinder::
  ad_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
          const std::vector<std::vector<MX> >& aseed,
          std::vector<std::vector<MX> >& asens,
          bool always_inline, bool never_inline) const {

    // Number of directional derivatives
    casadi_int nadj = aseed.size();
    asens.resize(nadj);

    // Quick return if no seeds
    if (nadj==0) return;

    // Get expression of Jacobian
    vector<MX> f_arg(arg);
    f_arg[iin_] = res.at(iout_);
    Function jac = get_function("jac_f_z");
    MX J = jac(f_arg).front();

    // Get adjoint seeds for calling f
    vector<MX> f_res(res);
    f_res[iout_] = MX(size_in(iin_)); // zero residual
    vector<vector<MX> > f_aseed(nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      f_aseed[d].resize(n_out_);
      for (casadi_int i=0; i<n_out_; ++i) f_aseed[d][i] = i==iout_ ? f_res[iout_] : aseed[d][i];
    }

    // Propagate dependencies from auxiliary outputs
    vector<MX> rhs(nadj);
    vector<vector<MX> > asens_aux;
    if (n_out_>1) {
      oracle_->call_reverse(f_arg, f_res, f_aseed, asens_aux, always_inline, never_inline);
      for (casadi_int d=0; d<nadj; ++d) rhs[d] = vec(asens_aux[d][iin_] + aseed[d][iout_]);
    } else {
      for (casadi_int d=0; d<nadj; ++d) rhs[d] = vec(aseed[d][iout_]);
    }

    // Solve for all the adjoint seeds at once
    rhs = horzsplit(J->get_solve(-horzcat(rhs), true, linsol_));
    for (casadi_int d=0; d<nadj; ++d) {
      for (casadi_int i=0; i<n_out_; ++i) {
        if (i==iout_) {
          f_aseed[d][i] = reshape(rhs[d], size_out(i));
        } else {
          // Avoid counting the auxiliary seeds twice
          f_aseed[d][i] = MX(size_out(i));
        }
      }
    }

    // No dependency on guess (1)
    vector<MX> tmp(nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d].resize(n_in_);
      tmp[d] = asens[d][iin_].is_empty(true) ? MX(size_in(iin_)) : asens[d][iin_];
    }

    // Propagate through f_
    oracle_->call_reverse(f_arg, f_res, f_aseed, asens, always_inline, never_inline);

    // No dependency on guess (2)
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d][iin_] = tmp[d];
    }

    // Add contribution from auxiliary outputs
    if (n_out_>1) {
      for (casadi_int d=0; d<nadj; ++d) {
        for (casadi_int i=0; i<n_in_; ++i) if (i!=iin_) asens[d][i] += asens_aux[d][i];
      }
    }
  }

  Dict Rootfinder::get_stats(void* mem) const {
    Dict stats = OracleFunction::get_stats(mem);
    auto m = static_cast<RootfinderMemory*>(mem);
    stats["success"] = m->success;
    return stats;
  }

} // namespace casadi
