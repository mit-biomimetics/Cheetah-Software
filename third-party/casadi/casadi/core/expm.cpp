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


#include "expm_impl.hpp"
#include "matrix.hpp"
#include "sparsity_interface.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_expm(const string& name) {
    return Expm::has_plugin(name);
  }

  void load_expm(const string& name) {
    Expm::load_plugin(name);
  }

  string doc_expm(const string& name) {
    return Expm::getPlugin(name).doc;
  }

  Function expmsol(const string& name, const string& solver,
                const Sparsity& A, const Dict& opts) {
    return Function::create(Expm::instantiate(name, solver, A), opts);
  }

  casadi_int expm_n_in() {
    return 2;
  }

  casadi_int expm_n_out() {
    return 1;
  }

  // Constructor
  Expm::Expm(const std::string& name, const Sparsity &A)
    : FunctionInternal(name), A_(Sparsity::dense(A.size1(), A.size2())) {

    casadi_assert_dev(A.is_square());

  }

  Sparsity Expm::get_sparsity_in(casadi_int i) {
    switch (i) {
      case 0:
        return A_;
      case 1:
        return Sparsity::dense(1, 1);
      default: break;
    }
    return Sparsity();
  }

  Sparsity Expm::get_sparsity_out(casadi_int i) {
    switch (i) {
      case 0:
        return A_;
      default: break;
    }
    return Sparsity();
  }

  Options Expm::options_
  = {{&FunctionInternal::options_},
     {{"const_A",
       {OT_BOOL,
        "Assume A is constant. Default: false."}}
     }
  };

  void Expm::init(const Dict& opts) {
    // Call the init method of the base class
    FunctionInternal::init(opts);

    const_A_ = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="const_A") {
        const_A_ = op.second;
      }
    }

  }

  Function Expm::get_forward(casadi_int nfwd, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const {
    MX A = MX::sym("A", A_);
    MX t = MX::sym("t");
    MX Y = MX::sym("Y", A_);
    MX Adot = MX::sym("Adot", A_);
    MX tdot = MX::sym("tdot");

    MX Ydot = mtimes(A, Y)*tdot;

    if (!const_A_) {
      DM N = DM::zeros(A_.size());

      MX extended = MX::blockcat({{A, Adot}, {N, A}});
      MX R = expm(extended*t);

      Ydot += R(Slice(0, A_.size1()), Slice(A_.size1(), 2*A_.size1()));
    }

    Function ret = Function(name, {A, t, Y, Adot, tdot}, {Ydot});

    return ret.map(name, "serial", nfwd,
      std::vector<casadi_int>{0, 1, 2}, std::vector<casadi_int>{});

  }

  Function Expm::get_reverse(casadi_int nadj, const std::string& name,
                               const std::vector<std::string>& inames,
                               const std::vector<std::string>& onames,
                               const Dict& opts) const {
    MX A = MX::sym("A", A_);
    MX t = MX::sym("t");
    MX Y = MX::sym("Y", A_);
    MX Ybar = MX::sym("Ybar", A_);

    MX tbar = sum2(sum1(Ybar*mtimes(A, Y)));
    MX Abar;
    if (const_A_) {
      Abar = MX(Sparsity(A_.size()));
    } else {
      DM N = DM::zeros(A_.size());

      MX At = A.T();
      MX extended = MX::blockcat({{At, Ybar}, {N, At}});
      MX R = expm(extended*t);

      Abar = R(Slice(0, A_.size1()), Slice(A_.size1(), 2*A_.size1()));
    }
    Function ret = Function(name, {A, t, Y, Ybar}, {Abar, tbar});

    return ret.map(name, "serial", nadj,
      std::vector<casadi_int>{0, 1, 2}, std::vector<casadi_int>{});
  }

  Sparsity Expm::getJacSparsity(casadi_int iind, casadi_int oind, bool symmetric) const {
    if (const_A_ && iind==0) {
      return Sparsity(nnz_out(oind), nnz_in(iind));
    }
    return Sparsity::dense(nnz_out(oind), nnz_in(iind));
  }

  Expm::~Expm() {
  }

  std::map<std::string, Expm::Plugin> Expm::solvers_;

  const std::string Expm::infix_ = "expm";

} // namespace casadi
