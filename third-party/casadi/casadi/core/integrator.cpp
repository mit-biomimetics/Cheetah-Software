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


#include "integrator_impl.hpp"
#include "casadi_misc.hpp"

using namespace std;
namespace casadi {

  bool has_integrator(const string& name) {
    return Integrator::has_plugin(name);
  }

  void load_integrator(const string& name) {
    Integrator::load_plugin(name);
  }

  string doc_integrator(const string& name) {
    return Integrator::getPlugin(name).doc;
  }

  Function integrator(const string& name, const string& solver,
                      const SXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2oracle("dae", dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                      const MXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2oracle("dae", dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                      const Function& dae, const Dict& opts) {
    // Make sure that dae is sound
    if (dae.has_free()) {
      casadi_error("Cannot create '" + name + "' since " + str(dae.get_free()) + " are free.");
    }
    return Function::create(Integrator::getPlugin(solver).creator(name, dae), opts);
  }

  vector<string> integrator_in() {
    vector<string> ret(integrator_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_in(i);
    return ret;
  }

  vector<string> integrator_out() {
    vector<string> ret(integrator_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_out(i);
    return ret;
  }

  string integrator_in(casadi_int ind) {
    switch (static_cast<IntegratorInput>(ind)) {
    case INTEGRATOR_X0:  return "x0";
    case INTEGRATOR_P:   return "p";
    case INTEGRATOR_Z0:  return "z0";
    case INTEGRATOR_RX0: return "rx0";
    case INTEGRATOR_RP:  return "rp";
    case INTEGRATOR_RZ0: return "rz0";
    case INTEGRATOR_NUM_IN: break;
    }
    return string();
  }

  string integrator_out(casadi_int ind) {
    switch (static_cast<IntegratorOutput>(ind)) {
    case INTEGRATOR_XF:  return "xf";
    case INTEGRATOR_QF:  return "qf";
    case INTEGRATOR_ZF:  return "zf";
    case INTEGRATOR_RXF: return "rxf";
    case INTEGRATOR_RQF: return "rqf";
    case INTEGRATOR_RZF: return "rzf";
    case INTEGRATOR_NUM_OUT: break;
    }
    return string();
  }

  casadi_int integrator_n_in() {
    return INTEGRATOR_NUM_IN;
  }

  casadi_int integrator_n_out() {
    return INTEGRATOR_NUM_OUT;
  }

  Integrator::Integrator(const std::string& name, const Function& oracle)
    : OracleFunction(name, oracle) {

    // Negative number of parameters for consistancy checking
    np_ = -1;

    // Default options
    print_stats_ = false;
    output_t0_ = false;
    print_time_ = false;
  }

  Integrator::~Integrator() {
  }

  Sparsity Integrator::get_sparsity_in(casadi_int i) {
    switch (static_cast<IntegratorInput>(i)) {
    case INTEGRATOR_X0: return x();
    case INTEGRATOR_P: return p();
    case INTEGRATOR_Z0: return z();
    case INTEGRATOR_RX0: return repmat(rx(), 1, ntout_);
    case INTEGRATOR_RP: return repmat(rp(), 1, ntout_);
    case INTEGRATOR_RZ0: return repmat(rz(), 1, ntout_);
    case INTEGRATOR_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Integrator::get_sparsity_out(casadi_int i) {
    switch (static_cast<IntegratorOutput>(i)) {
    case INTEGRATOR_XF: return repmat(x(), 1, ntout_);
    case INTEGRATOR_QF: return repmat(q(), 1, ntout_);
    case INTEGRATOR_ZF: return repmat(z(), 1, ntout_);
    case INTEGRATOR_RXF: return rx();
    case INTEGRATOR_RQF: return rq();
    case INTEGRATOR_RZF: return rz();
    case INTEGRATOR_NUM_OUT: break;
    }
    return Sparsity();
  }

  int Integrator::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<IntegratorMemory*>(mem);

    // Reset statistics
    for (auto&& s : m->fstats) s.second.reset();
    m->fstats.at(name_).tic();

    // Read inputs
    const double* x0 = arg[INTEGRATOR_X0];
    const double* z0 = arg[INTEGRATOR_Z0];
    const double* p = arg[INTEGRATOR_P];
    const double* rx0 = arg[INTEGRATOR_RX0];
    const double* rz0 = arg[INTEGRATOR_RZ0];
    const double* rp = arg[INTEGRATOR_RP];
    arg += INTEGRATOR_NUM_IN;

    // Read outputs
    double* x = res[INTEGRATOR_XF];
    double* z = res[INTEGRATOR_ZF];
    double* q = res[INTEGRATOR_QF];
    double* rx = res[INTEGRATOR_RXF];
    double* rz = res[INTEGRATOR_RZF];
    double* rq = res[INTEGRATOR_RQF];
    res += INTEGRATOR_NUM_OUT;

    // Setup memory object
    setup(m, arg, res, iw, w);

    // Reset solver, take time to t0
    reset(m, grid_.front(), x0, z0, p);

    // Integrate forward
    for (casadi_int k=0; k<grid_.size(); ++k) {
      // Skip t0?
      if (k==0 && !output_t0_) continue;

      // Integrate forward
      advance(m, grid_[k], x, z, q);
      if (x) x += nx_;
      if (z) z += nz_;
      if (q) q += nq_;
    }

    // If backwards integration is needed
    if (nrx_>0) {
      // Integrate backward
      resetB(m, grid_.back(), rx0, rz0, rp);

      // Proceed to t0
      retreat(m, grid_.front(), rx, rz, rq);
    }

    // Finalize/print statistics
    m->fstats.at(name_).toc();
    if (print_stats_) print_stats(m);
    if (print_time_)  print_fstats(m);
    return 0;
  }

  Options Integrator::options_
  = {{&OracleFunction::options_},
     {{"expand",
       {OT_BOOL,
        "Replace MX with SX expressions in problem formulation [false]"}},
      {"print_stats",
       {OT_BOOL,
        "Print out statistics after integration"}},
      {"t0",
       {OT_DOUBLE,
        "Beginning of the time horizon"}},
      {"tf",
       {OT_DOUBLE,
        "End of the time horizon"}},
      {"grid",
       {OT_DOUBLEVECTOR,
        "Time grid"}},
      {"augmented_options",
       {OT_DICT,
        "Options to be passed down to the augmented integrator, if one is constructed."}},
      {"output_t0",
       {OT_BOOL,
        "Output the state at the initial time"}}
     }
  };

  void Integrator::init(const Dict& opts) {
    // Default (temporary) options
    double t0=0, tf=1;
    bool expand = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="expand") {
        expand = op.second;
      } else if (op.first=="output_t0") {
        output_t0_ = op.second;
      } else if (op.first=="print_stats") {
        print_stats_ = op.second;
      } else if (op.first=="grid") {
        grid_ = op.second;
      } else if (op.first=="augmented_options") {
        augmented_options_ = op.second;
      } else if (op.first=="t0") {
        t0 = op.second;
      } else if (op.first=="tf") {
        tf = op.second;
      }
    }

    // Replace MX oracle with SX oracle?
    if (expand) this->expand();

    // Store a copy of the options, for creating augmented integrators
    opts_ = opts;

    // If grid unset, default to [t0, tf]
    if (grid_.empty()) {
      grid_ = {t0, tf};
    }

    ngrid_ = grid_.size();
    ntout_ = output_t0_ ? ngrid_ : ngrid_-1;

    // Call the base class method
    OracleFunction::init(opts);

    // For sparsity pattern propagation
    alloc(oracle_);

    // Error if sparse input
    casadi_assert(x().is_dense(), "Sparse DAE not supported");
    casadi_assert(z().is_dense(), "Sparse DAE not supported");
    casadi_assert(p().is_dense(), "Sparse DAE not supported");
    casadi_assert(rx().is_dense(), "Sparse DAE not supported");
    casadi_assert(rz().is_dense(), "Sparse DAE not supported");
    casadi_assert(rp().is_dense(), "Sparse DAE not supported");

    // Get dimensions (excluding sensitivity equations)
    nx1_ = x().size1();
    nz1_ = z().size1();
    nq1_ = q().size1();
    np1_  = p().size1();
    nrx1_ = rx().size1();
    nrz1_ = rz().size1();
    nrp1_ = rp().size1();
    nrq1_ = rq().size1();

    // Get dimensions (including sensitivity equations)
    nx_ = x().nnz();
    nz_ = z().nnz();
    nq_ = q().nnz();
    np_  = p().nnz();
    nrx_ = rx().nnz();
    nrz_ = rz().nnz();
    nrp_ = rp().nnz();
    nrq_ = rq().nnz();

    // Number of sensitivities
    ns_ = x().size2()-1;

    // Get the sparsities of the forward and reverse DAE
    sp_jac_dae_ = sp_jac_dae();
    casadi_assert(!sp_jac_dae_.is_singular(),
                          "Jacobian of the forward problem is structurally rank-deficient. "
                          "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<"
                          + str(nx_+nz_));
    if (nrx_>0) {
      sp_jac_rdae_ = sp_jac_rdae();
      casadi_assert(!sp_jac_rdae_.is_singular(),
                            "Jacobian of the backward problem is structurally rank-deficient. "
                            "sprank(J)=" + str(sprank(sp_jac_rdae_)) + "<"
                            + str(nrx_+nrz_));
    }

    // Consistency check

    // Allocate sufficiently large work vectors
    alloc_w(nx_+nz_);
    alloc_w(nrx_+nrz_);
    alloc_w(nx_ + nz_ + nrx_ + nrz_, true);
  }

  int Integrator::init_mem(void* mem) const {
    if (OracleFunction::init_mem(mem)) return 1;

    auto m = static_cast<IntegratorMemory*>(mem);
    m->add_stat(name_);
    return 0;
  }

  template<typename MatType>
  std::map<string, MatType> Integrator::aug_fwd(casadi_int nfwd) const {
    if (verbose_) casadi_message(name_ + "::aug_fwd");

    // Get input expressions
    vector<MatType> arg = MatType::get_input(oracle_);
    vector<MatType> aug_x, aug_z, aug_p, aug_rx, aug_rz, aug_rp;
    MatType aug_t = arg.at(DE_T);
    aug_x.push_back(vec(arg.at(DE_X)));
    aug_z.push_back(vec(arg.at(DE_Z)));
    aug_p.push_back(vec(arg.at(DE_P)));
    aug_rx.push_back(vec(arg.at(DE_RX)));
    aug_rz.push_back(vec(arg.at(DE_RZ)));
    aug_rp.push_back(vec(arg.at(DE_RP)));

    // Get output expressions
    vector<MatType> res = oracle_(arg);
    vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad;
    aug_ode.push_back(vec(res.at(DE_ODE)));
    aug_alg.push_back(vec(res.at(DE_ALG)));
    aug_quad.push_back(vec(res.at(DE_QUAD)));
    aug_rode.push_back(vec(res.at(DE_RODE)));
    aug_ralg.push_back(vec(res.at(DE_RALG)));
    aug_rquad.push_back(vec(res.at(DE_RQUAD)));

    // Zero of time dimension
    MatType zero_t = MatType::zeros(t());

    // Forward directional derivatives
    vector<vector<MatType>> seed(nfwd, vector<MatType>(DE_NUM_IN));
    for (casadi_int d=0; d<nfwd; ++d) {
      seed[d][DE_T] = zero_t;
      string pref = "aug" + str(d) + "_";
      aug_x.push_back(vec(seed[d][DE_X] = MatType::sym(pref + "x", x())));
      aug_z.push_back(vec(seed[d][DE_Z] = MatType::sym(pref + "z", z())));
      aug_p.push_back(vec(seed[d][DE_P] = MatType::sym(pref + "p", p())));
      aug_rx.push_back(vec(seed[d][DE_RX] = MatType::sym(pref + "rx", rx())));
      aug_rz.push_back(vec(seed[d][DE_RZ] = MatType::sym(pref + "rz", rz())));
      aug_rp.push_back(vec(seed[d][DE_RP] = MatType::sym(pref + "rp", rp())));
    }

    // Calculate directional derivatives
    vector<vector<MatType>> sens;
    oracle_->call_forward(arg, res, seed, sens, true, false);

    // Collect sensitivity equations
    casadi_assert_dev(sens.size()==nfwd);
    for (casadi_int d=0; d<nfwd; ++d) {
      casadi_assert_dev(sens[d].size()==DE_NUM_OUT);
      aug_ode.push_back(vec(project(sens[d][DE_ODE], x())));
      aug_alg.push_back(vec(project(sens[d][DE_ALG], z())));
      aug_quad.push_back(vec(project(sens[d][DE_QUAD], q())));
      aug_rode.push_back(vec(project(sens[d][DE_RODE], rx())));
      aug_ralg.push_back(vec(project(sens[d][DE_RALG], rz())));
      aug_rquad.push_back(vec(project(sens[d][DE_RQUAD], rq())));
    }

    // Construct return object
    std::map<string, MatType> ret;
    ret["t"] = aug_t;
    ret["x"] = horzcat(aug_x);
    ret["z"] = horzcat(aug_z);
    ret["p"] = horzcat(aug_p);
    ret["ode"] = horzcat(aug_ode);
    ret["alg"] = horzcat(aug_alg);
    ret["quad"] = horzcat(aug_quad);
    ret["rx"] = horzcat(aug_rx);
    ret["rz"] = horzcat(aug_rz);
    ret["rp"] = horzcat(aug_rp);
    ret["rode"] = horzcat(aug_rode);
    ret["ralg"] = horzcat(aug_ralg);
    ret["rquad"] = horzcat(aug_rquad);
    return ret;
  }

  template<typename MatType>
  std::map<string, MatType> Integrator::aug_adj(casadi_int nadj) const {
    if (verbose_) casadi_message(name_ + "::aug_adj");

    // Get input expressions
    vector<MatType> arg = MatType::get_input(oracle_);
    vector<MatType> aug_x, aug_z, aug_p, aug_rx, aug_rz, aug_rp;
    MatType aug_t = arg.at(DE_T);
    aug_x.push_back(vec(arg.at(DE_X)));
    aug_z.push_back(vec(arg.at(DE_Z)));
    aug_p.push_back(vec(arg.at(DE_P)));
    aug_rx.push_back(vec(arg.at(DE_RX)));
    aug_rz.push_back(vec(arg.at(DE_RZ)));
    aug_rp.push_back(vec(arg.at(DE_RP)));

    // Get output expressions
    vector<MatType> res = oracle_(arg);
    vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad;
    aug_ode.push_back(vec(res.at(DE_ODE)));
    aug_alg.push_back(vec(res.at(DE_ALG)));
    aug_quad.push_back(vec(res.at(DE_QUAD)));
    aug_rode.push_back(vec(res.at(DE_RODE)));
    aug_ralg.push_back(vec(res.at(DE_RALG)));
    aug_rquad.push_back(vec(res.at(DE_RQUAD)));

    // Zero of time dimension
    MatType zero_t = MatType::zeros(t());

    // Reverse mode directional derivatives
    vector<vector<MatType>> seed(nadj, vector<MatType>(DE_NUM_OUT));
    for (casadi_int d=0; d<nadj; ++d) {
      string pref = "aug" + str(d) + "_";
      aug_rx.push_back(vec(seed[d][DE_ODE] = MatType::sym(pref + "ode", x())));
      aug_rz.push_back(vec(seed[d][DE_ALG] = MatType::sym(pref + "alg", z())));
      aug_rp.push_back(vec(seed[d][DE_QUAD] = MatType::sym(pref + "quad", q())));
      aug_x.push_back(vec(seed[d][DE_RODE] = MatType::sym(pref + "rode", rx())));
      aug_z.push_back(vec(seed[d][DE_RALG] = MatType::sym(pref + "ralg", rz())));
      aug_p.push_back(vec(seed[d][DE_RQUAD] = MatType::sym(pref + "rquad", rq())));
    }

    // Calculate directional derivatives
    vector<vector<MatType>> sens;
    oracle_->call_reverse(arg, res, seed, sens, true, false);

    // Collect sensitivity equations
    casadi_assert_dev(sens.size()==nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      casadi_assert_dev(sens[d].size()==DE_NUM_IN);
      aug_rode.push_back(vec(project(sens[d][DE_X], x())));
      aug_ralg.push_back(vec(project(sens[d][DE_Z], z())));
      aug_rquad.push_back(vec(project(sens[d][DE_P], p())));
      aug_ode.push_back(vec(project(sens[d][DE_RX], rx())));
      aug_alg.push_back(vec(project(sens[d][DE_RZ], rz())));
      aug_quad.push_back(vec(project(sens[d][DE_RP], rp())));
    }

    // Construct return object
    std::map<string, MatType> ret;
    ret["t"] = aug_t;
    ret["x"] = vertcat(aug_x);
    ret["z"] = vertcat(aug_z);
    ret["p"] = vertcat(aug_p);
    ret["ode"] = vertcat(aug_ode);
    ret["alg"] = vertcat(aug_alg);
    ret["quad"] = vertcat(aug_quad);
    ret["rx"] = vertcat(aug_rx);
    ret["rz"] = vertcat(aug_rz);
    ret["rp"] = vertcat(aug_rp);
    ret["rode"] = vertcat(aug_rode);
    ret["ralg"] = vertcat(aug_ralg);
    ret["rquad"] = vertcat(aug_rquad);

    // Make sure that forward problem does not depend on backward states
    Function f("f", {ret["t"], ret["x"], ret["z"], ret["p"]},
                    {ret["ode"], ret["alg"], ret["quad"]});
    if (f.has_free()) {
      // Replace dependencies of rx, rz and rp with zeros
      f = Function("f", {ret["t"], ret["x"], ret["z"], ret["p"],
                         ret["rx"], ret["rz"], ret["rp"]},
                        {ret["ode"], ret["alg"], ret["quad"]});
      vector<MatType> v = {ret["t"], ret["x"], ret["z"], ret["p"], 0, 0, 0};
      v = f(v);
      ret["ode"] = v.at(0);
      ret["alg"] = v.at(1);
      ret["quad"] = v.at(2);
    }

    return ret;
  }

  int Integrator::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::sp_forward");

    // Work vectors
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;
    bvec_t *tmp_rx = w; w += nrx_;
    bvec_t *tmp_rz = w; w += nrz_;

    // Propagate forward
    const bvec_t** arg1 = arg+n_in_;
    fill_n(arg1, static_cast<size_t>(DE_NUM_IN), nullptr);
    arg1[DE_X] = arg[INTEGRATOR_X0];
    arg1[DE_P] = arg[INTEGRATOR_P];
    bvec_t** res1 = res+n_out_;
    fill_n(res1, static_cast<size_t>(DE_NUM_OUT), nullptr);
    res1[DE_ODE] = tmp_x;
    res1[DE_ALG] = tmp_z;
    oracle_(arg1, res1, iw, w, 0);
    if (arg[INTEGRATOR_X0]) {
      const bvec_t *tmp = arg[INTEGRATOR_X0];
      for (casadi_int i=0; i<nx_; ++i) tmp_x[i] |= *tmp++;
    }

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    copy_n(tmp_x, nx_+nz_, w);
    fill_n(tmp_x, nx_+nz_, 0);
    sp_jac_dae_.spsolve(tmp_x, w, false);

    // Get xf and zf
    if (res[INTEGRATOR_XF]) copy_n(tmp_x, nx_, res[INTEGRATOR_XF]);
    if (res[INTEGRATOR_ZF]) copy_n(tmp_z, nz_, res[INTEGRATOR_ZF]);

    // Propagate to quadratures
    if (nq_>0 && res[INTEGRATOR_QF]) {
      arg1[DE_X] = tmp_x;
      arg1[DE_Z] = tmp_z;
      res1[DE_ODE] = res1[DE_ALG] = nullptr;
      res1[DE_QUAD] = res[INTEGRATOR_QF];
      if (oracle_(arg1, res1, iw, w, 0)) return 1;
    }

    if (nrx_>0) {
      // Propagate through g
      fill_n(arg1, static_cast<size_t>(DE_NUM_IN), nullptr);
      arg1[DE_X] = tmp_x;
      arg1[DE_P] = arg[INTEGRATOR_P];
      arg1[DE_Z] = tmp_z;
      arg1[DE_RX] = arg[INTEGRATOR_X0];
      arg1[DE_RX] = arg[INTEGRATOR_RX0];
      arg1[DE_RP] = arg[INTEGRATOR_RP];
      fill_n(res1, static_cast<size_t>(DE_NUM_OUT), nullptr);
      res1[DE_RODE] = tmp_rx;
      res1[DE_RALG] = tmp_rz;
      oracle_(arg1, res1, iw, w, 0);
      if (arg[INTEGRATOR_RX0]) {
        const bvec_t *tmp = arg[INTEGRATOR_RX0];
        for (casadi_int i=0; i<nrx_; ++i) tmp_rx[i] |= *tmp++;
      }

      // "Solve" in order to resolve interdependencies (cf. Rootfinder)
      copy_n(tmp_rx, nrx_+nrz_, w);
      fill_n(tmp_rx, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(tmp_rx, w, false);

      // Get rxf and rzf
      if (res[INTEGRATOR_RXF]) copy_n(tmp_rx, nrx_, res[INTEGRATOR_RXF]);
      if (res[INTEGRATOR_RZF]) copy_n(tmp_rz, nrz_, res[INTEGRATOR_RZF]);

      // Propagate to quadratures
      if (nrq_>0 && res[INTEGRATOR_RQF]) {
        arg1[DE_RX] = tmp_rx;
        arg1[DE_RZ] = tmp_rz;
        res1[DE_RODE] = res1[DE_RALG] = nullptr;
        res1[DE_RQUAD] = res[INTEGRATOR_RQF];
        if (oracle_(arg1, res1, iw, w, 0)) return 1;
      }
    }
    return 0;
  }

  int Integrator::sp_reverse(bvec_t** arg, bvec_t** res,
      casadi_int* iw, bvec_t* w, void* mem) const {
    if (verbose_) casadi_message(name_ + "::sp_reverse");

    // Work vectors
    bvec_t** arg1 = arg+n_in_;
    bvec_t** res1 = res+n_out_;
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;

    // Shorthands
    bvec_t* x0 = arg[INTEGRATOR_X0];
    bvec_t* p = arg[INTEGRATOR_P];
    bvec_t* xf = res[INTEGRATOR_XF];
    bvec_t* zf = res[INTEGRATOR_ZF];
    bvec_t* qf = res[INTEGRATOR_QF];

    // Propagate from outputs to state vectors
    if (xf) {
      copy_n(xf, nx_, tmp_x);
      fill_n(xf, nx_, 0);
    } else {
      fill_n(tmp_x, nx_, 0);
    }
    if (zf) {
      copy_n(zf, nz_, tmp_z);
      fill_n(zf, nz_, 0);
    } else {
      fill_n(tmp_z, nz_, 0);
    }

    if (nrx_>0) {
      // Work vectors
      bvec_t *tmp_rx = w; w += nrx_;
      bvec_t *tmp_rz = w; w += nrz_;

      // Shorthands
      bvec_t* rx0 = arg[INTEGRATOR_RX0];
      bvec_t* rp = arg[INTEGRATOR_RP];
      bvec_t* rxf = res[INTEGRATOR_RXF];
      bvec_t* rzf = res[INTEGRATOR_RZF];
      bvec_t* rqf = res[INTEGRATOR_RQF];

      // Propagate from outputs to state vectors
      if (rxf) {
        copy_n(rxf, nrx_, tmp_rx);
        fill_n(rxf, nrx_, 0);
      } else {
        fill_n(tmp_rx, nrx_, 0);
      }
      if (rzf) {
        copy_n(rzf, nrz_, tmp_rz);
        fill_n(rzf, nrz_, 0);
      } else {
        fill_n(tmp_rz, nrz_, 0);
      }

      // Get dependencies from backward quadratures
      fill_n(res1, static_cast<size_t>(DE_NUM_OUT), nullptr);
      fill_n(arg1, static_cast<size_t>(DE_NUM_IN), nullptr);
      res1[DE_RQUAD] = rqf;
      arg1[DE_X] = tmp_x;
      arg1[DE_Z] = tmp_z;
      arg1[DE_P] = p;
      arg1[DE_RX] = tmp_rx;
      arg1[DE_RZ] = tmp_rz;
      arg1[DE_RP] = rp;
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;

      // Propagate interdependencies
      fill_n(w, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(w, tmp_rx, true);
      copy_n(w, nrx_+nrz_, tmp_rx);

      // Direct dependency rx0 -> rxf
      if (rx0) for (casadi_int i=0; i<nrx_; ++i) rx0[i] |= tmp_rx[i];

      // Indirect dependency via g
      res1[DE_RODE] = tmp_rx;
      res1[DE_RALG] = tmp_rz;
      res1[DE_RQUAD] = nullptr;
      arg1[DE_RX] = rx0;
      arg1[DE_RZ] = nullptr; // arg[INTEGRATOR_RZ0] is a guess, no dependency
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    }

    // Get dependencies from forward quadratures
    fill_n(res1, static_cast<size_t>(DE_NUM_OUT), nullptr);
    fill_n(arg1, static_cast<size_t>(DE_NUM_IN), nullptr);
    res1[DE_QUAD] = qf;
    arg1[DE_X] = tmp_x;
    arg1[DE_Z] = tmp_z;
    arg1[DE_P] = p;
    if (qf && nq_>0) {
      if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    }

    // Propagate interdependencies
    fill_n(w, nx_+nz_, 0);
    sp_jac_dae_.spsolve(w, tmp_x, true);
    copy_n(w, nx_+nz_, tmp_x);

    // Direct dependency x0 -> xf
    if (x0) for (casadi_int i=0; i<nx_; ++i) x0[i] |= tmp_x[i];

    // Indirect dependency through f
    res1[DE_ODE] = tmp_x;
    res1[DE_ALG] = tmp_z;
    res1[DE_QUAD] = nullptr;
    arg1[DE_X] = x0;
    arg1[DE_Z] = nullptr; // arg[INTEGRATOR_Z0] is a guess, no dependency
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
    return 0;
  }

  Function Integrator::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    if (verbose_) casadi_message(name_ + "::get_forward");

    // Integrator options
    Dict aug_opts = getDerivativeOptions(true);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Create integrator for augmented DAE
    Function aug_dae;
    string aug_prefix = "fsens" + str(nfwd) + "_";
    string dae_name = aug_prefix + oracle_.name();
    Dict dae_opts = {{"derivative_of", oracle_}};
    if (oracle_.is_a("SXFunction")) {
      aug_dae = map2oracle(dae_name, aug_fwd<SX>(nfwd));
    } else {
      aug_dae = map2oracle(dae_name, aug_fwd<MX>(nfwd));
    }
    aug_opts["derivative_of"] = self();
    Function aug_int = integrator(aug_prefix + name_, plugin_name(),
      aug_dae, aug_opts);

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN*(1+nfwd) + INTEGRATOR_NUM_OUT);

    // Augmented state
    vector<MX> x0_aug, p_aug, z0_aug, rx0_aug, rp_aug, rz0_aug;

    // Add nondifferentiated inputs and forward seeds
    for (casadi_int dir=-1; dir<nfwd; ++dir) {
      // Suffix
      string suff;
      if (dir>=0) suff = "_" + str(dir);

      // Augmented problem
      vector<MX> din(INTEGRATOR_NUM_IN);
      x0_aug.push_back(vec(din[INTEGRATOR_X0] = MX::sym("x0" + suff, x())));
      p_aug.push_back(vec(din[INTEGRATOR_P] = MX::sym("p" + suff, p())));
      z0_aug.push_back(vec(din[INTEGRATOR_Z0] = MX::sym("z0" + suff, z())));
      rx0_aug.push_back(vec(din[INTEGRATOR_RX0] = MX::sym("rx0" + suff, rx())));
      rp_aug.push_back(vec(din[INTEGRATOR_RP] = MX::sym("rp" + suff, rp())));
      rz0_aug.push_back(vec(din[INTEGRATOR_RZ0] = MX::sym("rz0" + suff, rz())));
      ret_in.insert(ret_in.end(), din.begin(), din.end());

      // Dummy outputs
      if (dir==-1) {
        vector<MX> dout(INTEGRATOR_NUM_OUT);
        dout[INTEGRATOR_XF]  = MX::sym("xf_dummy", Sparsity(size_out(INTEGRATOR_XF)));
        dout[INTEGRATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
        dout[INTEGRATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
        dout[INTEGRATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
        dout[INTEGRATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
        dout[INTEGRATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
        ret_in.insert(ret_in.end(), dout.begin(), dout.end());
      }
    }

    // Call the integrator
    vector<MX> integrator_in(INTEGRATOR_NUM_IN);
    integrator_in[INTEGRATOR_X0] = horzcat(x0_aug);
    integrator_in[INTEGRATOR_P] = horzcat(p_aug);
    integrator_in[INTEGRATOR_Z0] = horzcat(z0_aug);
    integrator_in[INTEGRATOR_RX0] = horzcat(rx0_aug);
    integrator_in[INTEGRATOR_RP] = horzcat(rp_aug);
    integrator_in[INTEGRATOR_RZ0] = horzcat(rz0_aug);
    vector<MX> integrator_out = aug_int(integrator_in);
    for (auto&& e : integrator_out) {
      // Workaround
      if (e.size2()!=1+nfwd) e = reshape(e, -1, 1+nfwd);
    }

    // Augmented results
    vector<casadi_int> offset = range(1+nfwd+1);
    vector<MX> xf_aug = horzsplit(integrator_out[INTEGRATOR_XF], offset);
    vector<MX> qf_aug = horzsplit(integrator_out[INTEGRATOR_QF], offset);
    vector<MX> zf_aug = horzsplit(integrator_out[INTEGRATOR_ZF], offset);
    vector<MX> rxf_aug = horzsplit(integrator_out[INTEGRATOR_RXF], offset);
    vector<MX> rqf_aug = horzsplit(integrator_out[INTEGRATOR_RQF], offset);
    vector<MX> rzf_aug = horzsplit(integrator_out[INTEGRATOR_RZF], offset);

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(INTEGRATOR_NUM_OUT*nfwd);

    // Collect the forward sensitivities
    vector<MX> dd(INTEGRATOR_NUM_IN);
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      dd[INTEGRATOR_XF]  = reshape(xf_aug.at(dir+1), x().size());
      dd[INTEGRATOR_QF]  = reshape(qf_aug.at(dir+1), q().size());
      dd[INTEGRATOR_ZF]  = reshape(zf_aug.at(dir+1), z().size());
      dd[INTEGRATOR_RXF] = reshape(rxf_aug.at(dir+1), rx().size());
      dd[INTEGRATOR_RQF] = reshape(rqf_aug.at(dir+1), rq().size());
      dd[INTEGRATOR_RZF] = reshape(rzf_aug.at(dir+1), rz().size());
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }

    // Concatenate forward seeds
    vector<MX> v(nfwd);
    auto r_it = ret_in.begin() + n_in_ + n_out_;
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = *(r_it + d*n_in_);
      *r_it++ = horzcat(v);
    }
    ret_in.resize(n_in_ + n_out_ + n_in_);

    // Concatenate forward sensitivites
    r_it = ret_out.begin();
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nfwd; ++d) v[d] = *(r_it + d*n_out_);
      *r_it++ = horzcat(v);
    }
    ret_out.resize(n_out_);

    // Create derivative function and return
    return Function(name, ret_in, ret_out, inames, onames, opts);
  }

  Function Integrator::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    if (verbose_) casadi_message(name_ + "::get_reverse");

    // Integrator options
    Dict aug_opts = getDerivativeOptions(false);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Create integrator for augmented DAE
    Function aug_dae;
    string aug_prefix = "asens" + str(nadj) + "_";
    string dae_name = aug_prefix + oracle_.name();
    if (oracle_.is_a("SXFunction")) {
      aug_dae = map2oracle(dae_name, aug_adj<SX>(nadj));
    } else {
      aug_dae = map2oracle(dae_name, aug_adj<MX>(nadj));
    }
    aug_opts["derivative_of"] = self();
    Function aug_int = integrator(aug_prefix + name_, plugin_name(),
      aug_dae, aug_opts);

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN + INTEGRATOR_NUM_OUT*(1+nadj));

    // Augmented state
    vector<MX> x0_aug, p_aug, z0_aug, rx0_aug, rp_aug, rz0_aug;

    // Inputs or forward/adjoint seeds in one direction
    vector<MX> dd(INTEGRATOR_NUM_IN);
    x0_aug.push_back(vec(dd[INTEGRATOR_X0] = MX::sym("x0", x())));
    p_aug.push_back(vec(dd[INTEGRATOR_P] = MX::sym("p", p())));
    z0_aug.push_back(vec(dd[INTEGRATOR_Z0] = MX::sym("r0", z())));
    rx0_aug.push_back(vec(dd[INTEGRATOR_RX0] = MX::sym("rx0", rx())));
    rp_aug.push_back(vec(dd[INTEGRATOR_RP] = MX::sym("rp", rp())));
    rz0_aug.push_back(vec(dd[INTEGRATOR_RZ0] = MX::sym("rz0", rz())));
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add dummy inputs (outputs of the nondifferentiated funciton)
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    dd[INTEGRATOR_XF]  = MX::sym("xf_dummy", Sparsity(x().size()));
    dd[INTEGRATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
    dd[INTEGRATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
    dd[INTEGRATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
    dd[INTEGRATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
    dd[INTEGRATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add adjoint seeds
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    for (casadi_int dir=0; dir<nadj; ++dir) {
      // Suffix
      string suff;
      if (dir>=0) suff = "_" + str(dir);

      // Augmented problem
      rx0_aug.push_back(vec(dd[INTEGRATOR_XF] = MX::sym("xf" + suff, x())));
      rp_aug.push_back(vec(dd[INTEGRATOR_QF] = MX::sym("qf" + suff, q())));
      rz0_aug.push_back(vec(dd[INTEGRATOR_ZF] = MX::sym("zf" + suff, z())));
      x0_aug.push_back(vec(dd[INTEGRATOR_RXF] = MX::sym("rxf" + suff, rx())));
      p_aug.push_back(vec(dd[INTEGRATOR_RQF] = MX::sym("rqf" + suff, rq())));
      z0_aug.push_back(vec(dd[INTEGRATOR_RZF] = MX::sym("rzf" + suff, rz())));
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());
    }

    // Call the integrator
    vector<MX> integrator_in(INTEGRATOR_NUM_IN);
    integrator_in[INTEGRATOR_X0] = vertcat(x0_aug);
    integrator_in[INTEGRATOR_P] = vertcat(p_aug);
    integrator_in[INTEGRATOR_Z0] = vertcat(z0_aug);
    integrator_in[INTEGRATOR_RX0] = vertcat(rx0_aug);
    integrator_in[INTEGRATOR_RP] = vertcat(rp_aug);
    integrator_in[INTEGRATOR_RZ0] = vertcat(rz0_aug);
    vector<MX> integrator_out = aug_int(integrator_in);

    // Get offset in the splitted problem
    vector<casadi_int> off_x = {0, x().numel()};
    vector<casadi_int> off_z = {0, z().numel()};
    vector<casadi_int> off_q = {0, q().numel()};
    vector<casadi_int> off_p = {0, p().numel()};
    vector<casadi_int> off_rx = {0, rx().numel()};
    vector<casadi_int> off_rz = {0, rz().numel()};
    vector<casadi_int> off_rq = {0, rq().numel()};
    vector<casadi_int> off_rp = {0, rp().numel()};
    for (casadi_int dir=0; dir<nadj; ++dir) {
      off_x.push_back(off_x.back() + rx().numel());
      off_z.push_back(off_z.back() + rz().numel());
      off_q.push_back(off_q.back() + rp().numel());
      off_p.push_back(off_p.back() + rq().numel());
      off_rx.push_back(off_rx.back() + x().numel());
      off_rz.push_back(off_rz.back() + z().numel());
      off_rq.push_back(off_rq.back() + p().numel());
      off_rp.push_back(off_rp.back() + q().numel());
    }

    // Augmented results
    vector<MX> xf_aug = vertsplit(integrator_out[INTEGRATOR_XF], off_x);
    vector<MX> qf_aug = vertsplit(integrator_out[INTEGRATOR_QF], off_q);
    vector<MX> zf_aug = vertsplit(integrator_out[INTEGRATOR_ZF], off_z);
    vector<MX> rxf_aug = vertsplit(integrator_out[INTEGRATOR_RXF], off_rx);
    vector<MX> rqf_aug = vertsplit(integrator_out[INTEGRATOR_RQF], off_rq);
    vector<MX> rzf_aug = vertsplit(integrator_out[INTEGRATOR_RZF], off_rz);

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(INTEGRATOR_NUM_IN*nadj);

    // Collect the adjoint sensitivities
    dd.resize(INTEGRATOR_NUM_IN);
    fill(dd.begin(), dd.end(), MX());
    for (casadi_int dir=0; dir<nadj; ++dir) {
      dd[INTEGRATOR_X0]  = reshape(rxf_aug.at(dir+1), x().size());
      dd[INTEGRATOR_P]   = reshape(rqf_aug.at(dir+1), p().size());
      dd[INTEGRATOR_Z0]  = reshape(rzf_aug.at(dir+1), z().size());
      dd[INTEGRATOR_RX0] = reshape(xf_aug.at(dir+1), rx().size());
      dd[INTEGRATOR_RP]  = reshape(qf_aug.at(dir+1), rp().size());
      dd[INTEGRATOR_RZ0] = reshape(zf_aug.at(dir+1), rz().size());
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }

    // Concatenate forward seeds
    vector<MX> v(nadj);
    auto r_it = ret_in.begin() + n_in_ + n_out_;
    for (casadi_int i=0; i<n_out_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = *(r_it + d*n_out_);
      *r_it++ = horzcat(v);
    }
    ret_in.resize(n_in_ + n_out_ + n_out_);

    // Concatenate forward sensitivites
    r_it = ret_out.begin();
    for (casadi_int i=0; i<n_in_; ++i) {
      for (casadi_int d=0; d<nadj; ++d) v[d] = *(r_it + d*n_in_);
      *r_it++ = horzcat(v);
    }
    ret_out.resize(n_in_);

    // Create derivative function and return
    return Function(name, ret_in, ret_out, inames, onames, opts);
  }

  Dict Integrator::getDerivativeOptions(bool fwd) const {
    // Copy all options
    return opts_;
  }

  Sparsity Integrator::sp_jac_dae() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = oracle_.sparsity_jac(DE_X, DE_ODE);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nx_);

    // Quick return if no algebraic variables
    if (nz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = oracle_.sparsity_jac(DE_Z, DE_ODE);
    Sparsity jac_alg_x = oracle_.sparsity_jac(DE_X, DE_ALG);
    Sparsity jac_alg_z = oracle_.sparsity_jac(DE_Z, DE_ALG);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  Sparsity Integrator::sp_jac_rdae() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = oracle_.sparsity_jac(DE_RX, DE_RODE);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nrx_);

    // Quick return if no algebraic variables
    if (nrz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = oracle_.sparsity_jac(DE_RZ, DE_RODE);
    Sparsity jac_alg_x = oracle_.sparsity_jac(DE_RX, DE_RALG);
    Sparsity jac_alg_z = oracle_.sparsity_jac(DE_RZ, DE_RALG);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  std::map<std::string, Integrator::Plugin> Integrator::solvers_;

  const std::string Integrator::infix_ = "integrator";

  void Integrator::setStopTime(IntegratorMemory* mem, double tf) const {
    casadi_error("setStopTime not defined for class " + class_name());
  }

  FixedStepIntegrator::FixedStepIntegrator(const std::string& name, const Function& dae)
    : Integrator(name, dae) {

    // Default options
    nk_ = 20;
  }

  FixedStepIntegrator::~FixedStepIntegrator() {
    clear_mem();
  }

  Options FixedStepIntegrator::options_
  = {{&Integrator::options_},
     {{"number_of_finite_elements",
       {OT_INT,
        "Number of finite elements"}}
     }
  };

  void FixedStepIntegrator::init(const Dict& opts) {
    // Call the base class init
    Integrator::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="number_of_finite_elements") {
        nk_ = op.second;
      }
    }

    // Number of finite elements and time steps
    casadi_assert_dev(nk_>0);
    h_ = static_cast<double>(grid_.back() - grid_.front())/static_cast<double>(nk_);

    // Setup discrete time dynamics
    setupFG();

    // Get discrete time dimensions
    nZ_ = F_.nnz_in(DAE_Z);
    nRZ_ =  G_.is_null() ? 0 : G_.nnz_in(RDAE_RZ);
  }

  int FixedStepIntegrator::init_mem(void* mem) const {
    if (Integrator::init_mem(mem)) return 1;
    auto m = static_cast<FixedStepMemory*>(mem);

    // Discrete time algebraic variable
    m->Z.resize(F_.nnz_in(DAE_Z));
    if (!G_.is_null()) m->RZ.resize(G_.nnz_in(RDAE_RZ));

    // Allocate tape if backward states are present
    if (nrx_>0) {
      m->x_tape.resize(nk_+1, vector<double>(nx_));
      m->Z_tape.resize(nk_, vector<double>(nZ_));
    }

    // Allocate state
    m->x.resize(nx_);
    m->z.resize(nz_);
    m->p.resize(np_);
    m->q.resize(nq_);
    m->rx.resize(nrx_);
    m->rz.resize(nrz_);
    m->rp.resize(nrp_);
    m->rq.resize(nrq_);
    m->x_prev.resize(nx_);
    m->Z_prev.resize(nZ_);
    m->q_prev.resize(nq_);
    m->rx_prev.resize(nrx_);
    m->RZ_prev.resize(nRZ_);
    m->rq_prev.resize(nrq_);
    return 0;
  }

  void FixedStepIntegrator::advance(IntegratorMemory* mem, double t,
                                    double* x, double* z, double* q) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Get discrete time sought
    casadi_int k_out = static_cast<casadi_int>(std::ceil((t - grid_.front())/h_));
    k_out = std::min(k_out, nk_); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert_dev(k_out>=0);

    // Explicit discrete time dynamics
    const Function& F = getExplicit();

    // Discrete dynamics function inputs ...
    fill_n(m->arg, F.n_in(), nullptr);
    m->arg[DAE_T] = &m->t;
    m->arg[DAE_X] = get_ptr(m->x_prev);
    m->arg[DAE_Z] = get_ptr(m->Z_prev);
    m->arg[DAE_P] = get_ptr(m->p);

    // ... and outputs
    fill_n(m->res, F.n_out(), nullptr);
    m->res[DAE_ODE] = get_ptr(m->x);
    m->res[DAE_ALG] = get_ptr(m->Z);
    m->res[DAE_QUAD] = get_ptr(m->q);

    // Take time steps until end time has been reached
    while (m->k<k_out) {
      // Update the previous step
      casadi_copy(get_ptr(m->x), nx_, get_ptr(m->x_prev));
      casadi_copy(get_ptr(m->Z), nZ_, get_ptr(m->Z_prev));
      casadi_copy(get_ptr(m->q), nq_, get_ptr(m->q_prev));

      // Take step
      F(m->arg, m->res, m->iw, m->w);
      casadi_axpy(nq_, 1., get_ptr(m->q_prev), get_ptr(m->q));

      // Tape
      if (nrx_>0) {
        casadi_copy(get_ptr(m->x), nx_, get_ptr(m->x_tape.at(m->k+1)));
        casadi_copy(get_ptr(m->Z), m->Z.size(), get_ptr(m->Z_tape.at(m->k)));
      }

      // Advance time
      m->k++;
      m->t = static_cast<double>(grid_.front()) + static_cast<double>(m->k)*h_;
    }

    // Return to user TODO(@jaeandersson): interpolate
    casadi_copy(get_ptr(m->x), nx_, x);
    casadi_copy(get_ptr(m->Z)+m->Z.size()-nz_, nz_, z);
    casadi_copy(get_ptr(m->q), nq_, q);
  }

  void FixedStepIntegrator::retreat(IntegratorMemory* mem, double t,
                                    double* rx, double* rz, double* rq) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Get discrete time sought
    casadi_int k_out = static_cast<casadi_int>(std::floor((t - grid_.front())/h_));
    //  make sure that rounding errors does not result in k_out>nk_
    k_out = std::max(k_out, casadi_int(0));
    casadi_assert_dev(k_out<=nk_);

    // Explicit discrete time dynamics
    const Function& G = getExplicitB();

    // Discrete dynamics function inputs ...
    fill_n(m->arg, G.n_in(), nullptr);
    m->arg[RDAE_T] = &m->t;
    m->arg[RDAE_P] = get_ptr(m->p);
    m->arg[RDAE_RX] = get_ptr(m->rx_prev);
    m->arg[RDAE_RZ] = get_ptr(m->RZ_prev);
    m->arg[RDAE_RP] = get_ptr(m->rp);

    // ... and outputs
    fill_n(m->res, G.n_out(), nullptr);
    m->res[RDAE_ODE] = get_ptr(m->rx);
    m->res[RDAE_ALG] = get_ptr(m->RZ);
    m->res[RDAE_QUAD] = get_ptr(m->rq);

    // Take time steps until end time has been reached
    while (m->k>k_out) {
      // Advance time
      m->k--;
      m->t = static_cast<double>(grid_.front()) + static_cast<double>(m->k)*h_;

      // Update the previous step
      casadi_copy(get_ptr(m->rx), nrx_, get_ptr(m->rx_prev));
      casadi_copy(get_ptr(m->RZ), nRZ_, get_ptr(m->RZ_prev));
      casadi_copy(get_ptr(m->rq), nrq_, get_ptr(m->rq_prev));

      // Take step
      m->arg[RDAE_X] = get_ptr(m->x_tape.at(m->k));
      m->arg[RDAE_Z] = get_ptr(m->Z_tape.at(m->k));
      G(m->arg, m->res, m->iw, m->w);
      casadi_axpy(nrq_, 1., get_ptr(m->rq_prev), get_ptr(m->rq));
    }

    // Return to user TODO(@jaeandersson): interpolate
    casadi_copy(get_ptr(m->rx), nrx_, rx);
    casadi_copy(get_ptr(m->RZ)+m->RZ.size()-nrz_, nrz_, rz);
    casadi_copy(get_ptr(m->rq), nrq_, rq);
  }

  void FixedStepIntegrator::
  reset(IntegratorMemory* mem, double t,
        const double* x, const double* z, const double* p) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(p, np_, get_ptr(m->p));

    // Update the state
    casadi_copy(x, nx_, get_ptr(m->x));
    casadi_copy(z, nz_, get_ptr(m->z));

    // Reset summation states
    casadi_fill(get_ptr(m->q), nq_, 0.);

    // Bring discrete time to the beginning
    m->k = 0;

    // Get consistent initial conditions
    casadi_fill(get_ptr(m->Z), m->Z.size(), numeric_limits<double>::quiet_NaN());

    // Add the first element in the tape
    if (nrx_>0) {
      casadi_copy(x, nx_, get_ptr(m->x_tape.at(0)));
    }
  }

  void FixedStepIntegrator::resetB(IntegratorMemory* mem, double t, const double* rx,
                                   const double* rz, const double* rp) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Update time
    m->t = t;

    // Set parameters
    casadi_copy(rp, nrp_, get_ptr(m->rp));

    // Update the state
    casadi_copy(rx, nrx_, get_ptr(m->rx));
    casadi_copy(rz, nrz_, get_ptr(m->rz));

    // Reset summation states
    casadi_fill(get_ptr(m->rq), nrq_, 0.);

    // Bring discrete time to the end
    m->k = nk_;

    // Get consistent initial conditions
    casadi_fill(get_ptr(m->RZ), m->RZ.size(), numeric_limits<double>::quiet_NaN());
  }

  ImplicitFixedStepIntegrator::
  ImplicitFixedStepIntegrator(const std::string& name, const Function& dae)
    : FixedStepIntegrator(name, dae) {
  }

  ImplicitFixedStepIntegrator::~ImplicitFixedStepIntegrator() {
  }

  Options ImplicitFixedStepIntegrator::options_
  = {{&FixedStepIntegrator::options_},
     {{"rootfinder",
       {OT_STRING,
        "An implicit function solver"}},
      {"rootfinder_options",
       {OT_DICT,
        "Options to be passed to the NLP Solver"}}
     }
  };

  void ImplicitFixedStepIntegrator::init(const Dict& opts) {
    // Call the base class init
    FixedStepIntegrator::init(opts);

    // Default (temporary) options
    std::string implicit_function_name = "newton";
    Dict rootfinder_options;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="rootfinder") {
        implicit_function_name = op.second.to_string();
      } else if (op.first=="rootfinder_options") {
        rootfinder_options = op.second;
      }
    }

    // Complete rootfinder dictionary
    rootfinder_options["implicit_input"] = DAE_Z;
    rootfinder_options["implicit_output"] = DAE_ALG;

    // Allocate a solver
    rootfinder_ = rootfinder(name_ + "_rootfinder", implicit_function_name,
                                  F_, rootfinder_options);
    alloc(rootfinder_);

    // Allocate a root-finding solver for the backward problem
    if (nRZ_>0) {
      // Options
      Dict backward_rootfinder_options = rootfinder_options;
      backward_rootfinder_options["implicit_input"] = RDAE_RZ;
      backward_rootfinder_options["implicit_output"] = RDAE_ALG;
      string backward_implicit_function_name = implicit_function_name;

      // Allocate a Newton solver
      backward_rootfinder_ =
        rootfinder(name_+ "_backward_rootfinder",
                   backward_implicit_function_name,
                   G_, backward_rootfinder_options);
      alloc(backward_rootfinder_);
    }
  }

  template<typename XType>
  Function Integrator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d, const Dict& opts) {
    std::vector<XType> de_in(DE_NUM_IN), de_out(DE_NUM_OUT);

    for (auto&& i : d) {
      if (i.first=="t") {
        de_in[DE_T]=i.second;
      } else if (i.first=="x") {
        de_in[DE_X]=i.second;
      } else if (i.first=="z") {
        de_in[DE_Z]=i.second;
      } else if (i.first=="p") {
        de_in[DE_P]=i.second;
      } else if (i.first=="rx") {
        de_in[DE_RX]=i.second;
      } else if (i.first=="rz") {
        de_in[DE_RZ]=i.second;
      } else if (i.first=="rp") {
        de_in[DE_RP]=i.second;
      } else if (i.first=="ode") {
        de_out[DE_ODE]=i.second;
      } else if (i.first=="alg") {
        de_out[DE_ALG]=i.second;
      } else if (i.first=="quad") {
        de_out[DE_QUAD]=i.second;
      } else if (i.first=="rode") {
        de_out[DE_RODE]=i.second;
      } else if (i.first=="ralg") {
        de_out[DE_RALG]=i.second;
      } else if (i.first=="rquad") {
        de_out[DE_RQUAD]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }

    // Make sure x and ode exist
    casadi_assert(!de_in[DE_X].is_empty(), "Ill-posed ODE - no state");

    // Number of right-hand-sides
    casadi_int nrhs = de_in[DE_X].size2();

    // Make sure consistent number of right-hand-sides
    for (bool b : {true, false}) {
      for (auto&& e : b ? de_in : de_out) {
        // Skip time
        if (&e == &de_in[DE_T]) continue;
        // Number of rows
        casadi_int nr = e.size1();
        // Make sure no change in number of elements
        casadi_assert(e.numel()==nr*nrhs, "Inconsistent number of rhs");
        e = reshape(e, nr, nrhs);
      }
    }

    // Consistent sparsity for x
    casadi_assert(de_in[DE_X].size()==de_out[DE_ODE].size(),
      "Dimension mismatch for 'ode'");
    de_out[DE_ODE] = project(de_out[DE_ODE], de_in[DE_X].sparsity());

    // Consistent sparsity for z
    casadi_assert(de_in[DE_Z].size()==de_out[DE_ALG].size(),
      "Dimension mismatch for 'alg'");
    de_out[DE_ALG] = project(de_out[DE_ALG], de_in[DE_Z].sparsity());

    // Consistent sparsity for rx
    casadi_assert(de_in[DE_RX].size()==de_out[DE_RODE].size(),
      "Dimension mismatch for 'rode'");
    de_out[DE_RODE] = project(de_out[DE_RODE], de_in[DE_RX].sparsity());

    // Consistent sparsity for rz
    casadi_assert(de_in[DE_RZ].size()==de_out[DE_RALG].size(),
      "Dimension mismatch for 'ralg'");
    de_out[DE_RALG] = project(de_out[DE_RALG], de_in[DE_RZ].sparsity());

    // Construct
    return Function(name, de_in, de_out, DE_INPUTS, DE_OUTPUTS, opts);
  }

} // namespace casadi
