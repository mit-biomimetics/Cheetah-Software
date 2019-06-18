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


#include "dple_impl.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_dple(const string& name) {
    return Dple::has_plugin(name);
  }

  void load_dple(const string& name) {
    Dple::load_plugin(name);
  }

  string doc_dple(const string& name) {
    return Dple::getPlugin(name).doc;
  }

  MX dplesol(const MX& A, const MX& V, const std::string& solver, const Dict& opts) {
    SpDict sp;
    sp["a"] = A.sparsity();
    sp["v"] = V.sparsity();
    Function f = dplesol("dplesol", solver, sp, opts);
    MXDict f_in;
    f_in["a"] = A;
    f_in["v"] = V;
    MXDict f_out = f(f_in);
    return f_out["p"];
  }

  CASADI_EXPORT MXVector dplesol(const MXVector& A, const MXVector& V, const std::string& solver,
    const Dict& opts) {
      casadi_assert(A.size()==V.size(),
        "dplesol: sizes of A vector (" + str(A.size()) + ") and V vector "
        "(" + str(V.size()) + ") must match.");
      std::vector<MX> Adense, Vdense;

      for (casadi_int i=0;i<A.size();++i) {
        Adense.push_back(densify(A[i]));
        Vdense.push_back(densify(V[i]));
      }

      MX ret = dplesol(diagcat(Adense), diagcat(Vdense), solver, opts);
      return diagsplit(ret, ret.size1()/A.size());
  }

  CASADI_EXPORT DMVector dplesol(const DMVector& A, const DMVector& V, const std::string& solver,
    const Dict& opts) {
      casadi_assert(A.size()==V.size(),
        "dplesol: sizes of A vector (" + str(A.size()) + ") and V vector "
        "(" + str(V.size()) + ") must match.");
      std::vector<DM> Adense, Vdense;

      for (casadi_int i=0;i<A.size();++i) {
        Adense.push_back(densify(A[i]));
        Vdense.push_back(densify(V[i]));
      }

      DM Afull = diagcat(Adense);
      DM Vfull = diagcat(Vdense);

      SpDict sp;
      sp["a"] = Afull.sparsity();
      sp["v"] = Vfull.sparsity();
      Function f = dplesol("dplesol", solver, sp, opts);
      DMDict f_in;
      f_in["a"] = Afull;
      f_in["v"] = Vfull;
      DMDict f_out = f(f_in);
      return diagsplit(f_out["p"], f_out["p"].size1()/A.size());
  }

  Function dplesol(const string& name, const string& solver,
                const SpDict& st, const Dict& opts) {
    return Function::create(Dple::instantiate(name, solver, st), opts);
  }

  vector<string> dple_in() {
    vector<string> ret(dple_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=dple_in(i);
    return ret;
  }

  vector<string> dple_out() {
    vector<string> ret(dple_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=dple_out(i);
    return ret;
  }

  string dple_in(casadi_int ind) {
    switch (static_cast<DpleInput>(ind)) {
    case DPLE_A:      return "a";
    case DPLE_V:      return "v";
    case DPLE_NUM_IN: break;
    }
    return string();
  }

  string dple_out(casadi_int ind) {
    switch (static_cast<DpleOutput>(ind)) {
      case DPLE_P:      return "p";
      case DPLE_NUM_OUT: break;
    }
    return string();
  }

  casadi_int dple_n_in() {
    return DPLE_NUM_IN;
  }

  casadi_int dple_n_out() {
    return DPLE_NUM_OUT;
  }

  // Constructor
  Dple::Dple(const std::string& name, const SpDict &st)
    : FunctionInternal(name) {
    for (auto i=st.begin(); i!=st.end(); ++i) {
      if (i->first=="a") {
        A_ = i->second;
      } else if (i->first=="v") {
        V_ = i->second;
      } else {
        casadi_error("Unrecognized field in Dple structure: " + str(i->first));
      }
    }

  }

  Sparsity Dple::get_sparsity_in(casadi_int i) {
    switch (static_cast<DpleInput>(i)) {
      case DPLE_A:
        return A_;
      case DPLE_V:
        return V_;
      case DPLE_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Dple::get_sparsity_out(casadi_int i) {
    switch (static_cast<DpleOutput>(i)) {
      case DPLE_P:
        return V_;
      case DPLE_NUM_OUT: break;
    }
    return Sparsity();
  }

  Options Dple::options_
  = {{&FunctionInternal::options_},
     {{"const_dim",
       {OT_BOOL,
        "Assume constant dimension of P"}},
      {"pos_def",
        {OT_BOOL,
         "Assume P positive definite"}},
      {"error_unstable",
        {OT_BOOL,
        "Throw an exception when it is detected that Product(A_i, i=N..1)"
        "has eigenvalues greater than 1-eps_unstable"}},
      {"eps_unstable",
        {OT_DOUBLE,
        "A margin for unstability detection"}}
     }
  };

  void Dple::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    // Default options
    const_dim_ = true;
    pos_def_ = false;
    error_unstable_ = false;
    eps_unstable_ = 1e-4;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="const_dim") {
        const_dim_ = op.second;
      } else if  (op.first=="pos_def") {
        pos_def_ = op.second;
      } else if  (op.first=="error_unstable") {
        error_unstable_ = op.second;
      } else if  (op.first=="eps_unstable") {
        eps_unstable_ = op.second;
      }
    }

    casadi_assert_dev(V_.size2() % V_.size1() == 0);
    nrhs_ = V_.size2() / V_.size1();
    casadi_assert_dev(nrhs_>=1);

    std::vector<Sparsity> Vs = horzsplit(V_, V_.size1());
    Sparsity Vref = Vs[0];
    casadi_assert(Vref.is_symmetric(),
      "V must be symmetric but got " + Vref.dim() + ".");

    for (auto&& s : Vs)
      casadi_assert_dev(s==Vref);

    casadi_assert(const_dim_, "Not implemented");

    casadi_int blocksize = Vref.colind()[1];
    K_ = Vref.size1()/blocksize;
    Sparsity block = Sparsity::dense(blocksize, blocksize);

    std::vector<Sparsity> blocks(K_, block);
    casadi_assert(Vref==diagcat(blocks), "Structure not recognised.");
    casadi_assert(A_==Vref, "Structure not recognised.");


  }

  Function Dple::get_forward(casadi_int nfwd, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const {
    // Symbolic A
    MX A = MX::sym("A", A_);
    Function Vdotf;
    {
      MX P = MX::sym("P", A_);
      MX Adot = MX::sym("P", A_);
      MX Vdot = MX::sym("P", A_);

      MX temp = mtimes(std::vector<MX>{Adot, P, A.T()}) +
                mtimes(std::vector<MX>{A, P, Adot.T()}) + Vdot;
      Vdotf = Function("PAVbar", {A, P, Adot, Vdot},
                { (temp+temp.T())/2});
    }

    MX P = MX::sym("P", V_);
    MX Adot = MX::sym("Adot", repmat(A_, 1, nfwd));
    MX Vdot = MX::sym("Vdot", repmat(V_, 1, nfwd));
    MX Qdot = Vdotf.map("map", "serial", nrhs_, {0, 2}, std::vector<casadi_int>{})
         .map("map", "serial", nfwd, {0, 1}, std::vector<casadi_int>{})({A, P, Adot, Vdot})[0];
    MX Pdot = dplesol(A, Qdot, plugin_name(), opts);
    MX V = MX::sym("V", Sparsity(size_in(DPLE_V))); // We dont need V
    return Function(name, {A, V, P, Adot, Vdot}, {Pdot}, inames, onames);

  }

  Function Dple::get_reverse(casadi_int nadj, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const {

    // Symbolic A
    MX A = MX::sym("A", A_);

    // Helper function to reverse, reverse-tranpose,
    // and reverse-symmetrize one block-diagonal matrix
    casadi_int n = A_.size1()/K_;
    std::vector<MX> ret = diagsplit(A, n);
    std::reverse(ret.begin(), ret.end());
    std::vector<MX> retT;
    std::vector<MX> retS;
    for (auto & e : ret) retT.push_back(e.T());
    for (auto & e : ret) retS.push_back((e+e.T())/2);
    Function revS = Function("revS", {A}, {diagcat(retS)});
    Function revT = Function("revT", {A}, {diagcat(retT)});
    Function rev  = Function("rev", {A}, {diagcat(ret)});

    // Function to compute the formula for Abar
    Function Abarf;
    {
      MX P = MX::sym("P", A_);
      MX Vbar_rev = MX::sym("Vbar", A_);
      MX A_rev = MX::sym("A", A_);

      Abarf = Function("PAVbar", {P, A_rev, Vbar_rev},
                {2*revT(mtimes(std::vector<MX>{rev(P)[0], A_rev, Vbar_rev}))[0]});
    }

    // original output
    MX P = MX::sym("P", V_);

    // Symbolic reverse seed for P
    MX Pbar = MX::sym("Pbar", repmat(V_, 1, nadj));
    // Symmetrize the seed
    MX Pbar_rev = revS.map(nrhs_).map(nadj)(Pbar)[0];

    // Reverse A for new dple
    MX A_rev = revT(A)[0];

    // Solver a dple with nrhs*nadj right-hand sides
    MX Vbar_rev = dplesol(A_rev, Pbar_rev, plugin_name(), opts);

    // Undo the reversal for Vbar
    MX Vbar = rev.map(nrhs_).map(nadj)(Vbar_rev)[0];

    MX Abar = Abarf.map("map", "serial", nrhs_, std::vector<casadi_int>{1}, {0}).
                    map("map", "serial", nadj, {0, 1}, std::vector<casadi_int>{})(
                      {P, A_rev, Vbar_rev})[0];

    MX V = MX::sym("V", Sparsity(size_in(DPLE_V))); // We dont need V
    return Function(name, {A, V, P, Pbar}, {Abar, Vbar}, inames, onames);
  }

  Dple::~Dple() {
  }

  std::map<std::string, Dple::Plugin> Dple::solvers_;

  const std::string Dple::infix_ = "dple";

} // namespace casadi
