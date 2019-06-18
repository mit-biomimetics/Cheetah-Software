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


#include "ampl_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include <cstdio>
#include <cstdlib>
#include <fstream>
using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_AMPL_EXPORT
  casadi_register_nlpsol_ampl(Nlpsol::Plugin* plugin) {
    plugin->creator = AmplInterface::creator;
    plugin->name = "ampl";
    plugin->doc = AmplInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &AmplInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_AMPL_EXPORT casadi_load_nlpsol_ampl() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_ampl);
  }

  AmplInterface::AmplInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  AmplInterface::~AmplInterface() {
    clear_mem();
  }

  Options AmplInterface::options_
  = {{&Nlpsol::options_},
     {{"solver",
       {OT_STRING,
        "AMPL solver binary"}}
     }
  };

  void AmplInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Set default options
    solver_ = "ipopt";

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="solver") {
        solver_ = op.first;
      }
    }

    // Extract the expressions
    casadi_assert(oracle().is_a("SXFunction"),
                  "Only SX supported currently.");
    vector<SX> xp = oracle().sx_in();
    vector<SX> fg = oracle()(xp);

    // Get x, p, f and g
    SX x = xp.at(NL_X);
    SX p = xp.at(NL_P);
    SX f = fg.at(NL_F);
    SX g = fg.at(NL_G);
    casadi_assert(p.is_empty(), "'p' currently not supported");

    // Names of the variables, constraints
    vector<string> x_name, g_name;
    for (casadi_int i=0; i<nx_; ++i) x_name.push_back("x[" + str(i) + "]");
    for (casadi_int i=0; i<ng_; ++i) g_name.push_back("g[" + str(i) + "]");
    casadi_int max_x_name = x_name.back().size();
    casadi_int max_g_name = g_name.empty() ? 0 : g_name.back().size();

    // Calculate the Jacobian, gradient
    Sparsity jac_g = SX::jacobian(g, x).sparsity();
    Sparsity jac_f = SX::jacobian(f, x).sparsity();

    // Extract the shared subexpressions
    vector<SX> ex = {f, g}, v, vdef;
    shared(ex, v, vdef);
    f = ex[0];
    g = ex[1];

    // Header
    nl_init_ << "g3 1 1 0\n";
    // Type of constraints
    nl_init_ << nx_ << " " // number of variables
       << ng_ << " " // number of constraints
       << 1 << " " // number of objectives
       << 0 << " " // number of ranges
       << 0 << " " // ?
       << 0 << "\n"; // number of logical constraints

    // Nonlinearity - assume all nonlinear for now TODO: Detect
    nl_init_ << ng_ << " "  // nonlinear constraints
       << 1 << "\n"; // nonlinear objectives

    // Network constraints
    nl_init_ << 0 << " " // nonlinear
       << 0 << "\n"; // linear

    // Nonlinear variables
    nl_init_ << nx_ << " " // in constraints
       << nx_ << " " // in objectives
       << nx_ << "\n"; // in both

    // Linear network ..
    nl_init_ << 0 << " " // .. variables ..
       << 0 << " " // .. arith ..
       << 0 << " " // .. functions ..
       << 0 << "\n"; // .. flags

    // Discrete variables
    nl_init_ << 0 << " " // binary
       << 0 << " " // integer
       << 0 << " " // nonlinear in both
       << 0 << " " // nonlinear in constraints
       << 0 << "\n"; // nonlinear in objective

    // Nonzeros in the Jacobian, gradients
    nl_init_ << jac_g.nnz() << " " // nnz in Jacobian
       << jac_f.nnz() << "\n"; // nnz in gradients

    // Maximum name length
    nl_init_ << max_x_name << " " // constraints
       << max_g_name << "\n"; // variables

    // Shared subexpressions
    nl_init_ << v.size() << " " // both
       << 0 << " " // constraints
       << 0 << " " // objective
       << 0 << " " // c1 - constaint, but linear?
       << 0 << "\n"; // o1 - objective, but linear?

    // Create a function which evaluates f and g
    Function F("F", {vertcat(v), x}, {vertcat(vdef), f, g},
              {"v", "x"}, {"vdef", "f", "g"});

    // Iterate over the algoritm
    vector<string> work(F.sz_w());

    // Loop over the algorithm
    for (casadi_int k=0; k<F.n_instructions(); ++k) {
      // Get the atomic operation
      casadi_int op = F.instruction_id(k);
      // Get the operation indices
      std::vector<casadi_int> o = F.instruction_output(k);
      casadi_int o0=-1, o1=-1, i0=-1, i1=-1;
      if (o.size()>0) o0 = o[0];
      if (o.size()>1) o1 = o[1];
      std::vector<casadi_int> i = F.instruction_input(k);
      if (i.size()>0) i0 = i[0];
      if (i.size()>1) i1 = i[1];
      switch (op) {
        case OP_CONST:
        work[o0] = "n" + str(F.instruction_constant(k)) + "\n";
        break;
        case OP_INPUT:
        work[o0] = "v" + str(i0*v.size() + i1) + "\n";
        break;
        case OP_OUTPUT:
        if (o0==0) {
          // Common subexpression
          nl_init_ << "V" << (x.nnz()+o1) << " 0 0\n" << work[i0];
        } else if (o0==1) {
          // Nonlinear objective term
          nl_init_ << "O" << o1 << " 0\n" << work[i0];
        } else {
          // Nonlinear constraint term
          nl_init_ << "C" << o1 << "\n" << work[i0];
        }
        break;
        case OP_ADD: work[o0] = "o0\n" + work[i0] + work[i1]; break;
        case OP_SUB: work[o0] = "o1\n" + work[i0] + work[i1]; break;
        case OP_MUL: work[o0] = "o2\n" + work[i0] + work[i1]; break;
        case OP_DIV: work[o0] = "o3\n" + work[i0] + work[i1]; break;
        case OP_SQ: work[o0] = "o5\n" + work[i0] + "n2\n"; break;
        case OP_POW: work[o0] = "o5\n" + work[i0] + work[i1]; break;
        case OP_FLOOR: work[o0] = "o13\n" + work[i0]; break;
        case OP_CEIL: work[o0] = "o14\n" + work[i0]; break;
        case OP_FABS: work[o0] = "o15\n" + work[i0]; break;
        case OP_NEG: work[o0] = "o16\n" + work[i0]; break;
        case OP_TANH: work[o0] = "o37\n" + work[i0]; break;
        case OP_TAN: work[o0] = "o38\n" + work[i0]; break;
        case OP_SQRT: work[o0] = "o39\n" + work[i0]; break;
        case OP_SINH: work[o0] = "o40\n" + work[i0]; break;
        case OP_SIN: work[o0] = "o41\n" + work[i0]; break;
        case OP_LOG: work[o0] = "o43\n" + work[i0]; break;
        case OP_EXP: work[o0] = "o44\n" + work[i0]; break;
        case OP_COSH: work[o0] = "o45\n" + work[i0]; break;
        case OP_COS: work[o0] = "o46\n" + work[i0]; break;
        case OP_ATANH: work[o0] = "o47\n" + work[i0]; break;
        case OP_ATAN2: work[o0] = "o48\n" + work[i0] + work[i1]; break;
        case OP_ATAN: work[o0] = "o49\n" + work[i0]; break;
        case OP_ASINH: work[o0] = "o50\n" + work[i0]; break;
        case OP_ASIN: work[o0] = "o51\n" + work[i0]; break;
        case OP_ACOSH: work[o0] = "o52\n" + work[i0]; break;
        case OP_ACOS: work[o0] = "o53\n" + work[i0]; break;
        default:
        if (casadi_math<double>::ndeps(op)==1) {
          casadi_error(casadi_math<double>::print(op, "x") + " not supported");
        } else {
          casadi_error(casadi_math<double>::print(op, "x", "y") + " not supported");
        }
      }
    }

    // k segments, cumulative column count in jac_g
    const casadi_int *colind = jac_g.colind(), *row = jac_g.row();
    nl_init_ << "k" << (nx_-1) << "\n";
    for (casadi_int i=1; i<nx_; ++i) nl_init_ << colind[i] << "\n";

    // J segments, rows in jac_g
    Sparsity sp = jac_g.T();
    colind = sp.colind(), row = sp.row();
    for (casadi_int i=0; i<ng_; ++i) {
      nl_init_ << "J" << i << " " << (colind[i+1]-colind[i]) << "\n";
      for (casadi_int k=colind[i]; k<colind[i+1]; ++k) {
        casadi_int r=row[k];
        nl_init_ << r << " " << 0 << "\n"; // no linear term
      }
    }

    // G segments, rows in jac_f
    sp = jac_f.T();
    colind = sp.colind(), row = sp.row();
    nl_init_ << "G" << 0 << " " << (colind[0+1]-colind[0]) << "\n";
    for (casadi_int k=colind[0]; k<colind[0+1]; ++k) {
      casadi_int r=row[k];
      nl_init_ << r << " " << 0 << "\n"; // no linear term
    }
  }

  int AmplInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    //auto m = static_cast<AmplInterfaceMemory*>(mem);

    return 0;
  }

  void AmplInterface::set_work(void* mem, const double**& arg, double**& res,
                                   casadi_int*& iw, double*& w) const {
    //auto m = static_cast<AmplInterfaceMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

  }

  int AmplInterface::solve(void* mem) const {
    auto m = static_cast<AmplInterfaceMemory*>(mem);

    // Create .nl file and add preamble
    std::string nlname = temporary_file("casadi_ampl_tmp", ".nl");
    ofstream nl(nlname, ofstream::out);
    if (verbose_) casadi_message("Opened " + nlname);
    nl << nl_init_.str();

    // Primal intial guess
    nl << "x" << nx_ << "\n";
    for (casadi_int i=0; i<nx_; ++i) {
      nl << i << " " << m->x[i] << "\n";
    }


    // Add constraint bounds
    nl << "r\n";
    for (casadi_int i=0; i<ng_; ++i) {
      double lbg = m->lbg ? m->lbg[i] : 0;
      double ubg = m->ubg ? m->ubg[i] : 0;
      if (isinf(lbg)) {
        if (isinf(ubg)) { // no constraint
          nl << "3\n";
        } else { // only upper
          nl << "1 " << ubg <<  "\n";
        }
      } else {
        if (isinf(ubg)) { // only lower
          nl << "2 " << lbg << "\n";
        } else if (ubg==lbg) { // equality
          nl << "4 " << lbg << "\n";
        } else { // range
          nl << "0 " << lbg << " " << ubg << "\n";
        }
      }
    }

    // Add variable bounds
    nl << "b\n";
    for (casadi_int i=0; i<nx_; ++i) {
      double lbx = m->lbx ? m->lbx[i] : 0;
      double ubx = m->ubx ? m->ubx[i] : 0;
      if (isinf(lbx)) {
        if (isinf(ubx)) { // no constraint
          nl << "3\n";
        } else { // only upper
          nl << "1 " << ubx <<  "\n";
        }
      } else {
        if (isinf(ubx)) { // only lower
          nl << "2 " << lbx << "\n";
        } else if (ubx==lbx) { // equality
          nl << "4 " << lbx << "\n";
        } else { // range
          nl << "0 " << lbx << " " << ubx << "\n";
        }
      }
    }

    // Close file
    nl.close();

    // Temporary name for the .sol file
    std::string solname = temporary_file("casadi_ampl_tmp", ".sol");

    // Temporary name for the solver output
    std::string outname = temporary_file("casadi_ampl_tmp", ".out");
    // Call executable
    string system_cmd = solver_ + " -o" + solname + " " + nlname + " > " + outname;
    int ret = system(system_cmd.c_str());
    casadi_assert_dev(ret==0);

    // Delete the nl file
    if (remove(nlname.c_str())!=0) {
      casadi_warning("Failed to remove " + nlname);
    }
    if (verbose_) casadi_message("Removed " + nlname);

    // Open .out file and dump to screen
    ifstream out(outname, ifstream::in);
    casadi_assert(out.is_open(), "Failed to open " + outname);
    string line;
    while (!out.eof()) {
      getline(out, line);
      uout() << line << "\n";
    }

    // Close and delete .out file
    out.close();
    if (remove(outname.c_str())!=0) {
      casadi_warning("Failed to remove " + outname);
    }
    if (verbose_) casadi_message("Removed " + outname);

    // Open .sol file
    ifstream sol(solname, ifstream::in);
    casadi_assert(sol.is_open(), "Failed to open " + solname);
    if (verbose_) casadi_message("Opened " + solname);

    // Get all the lines
    vector<string> sol_lines;
    while (!sol.eof()) {
      getline(sol, line);
      sol_lines.push_back(line);
    }

    // Get the primal solution
    for (casadi_int i=0; i<nx_; ++i) {
      istringstream s(sol_lines.at(sol_lines.size()-nx_+i-1));
      s >> m->x[i];
    }

    // Get the dual solution
    for (casadi_int i=0; i<ng_; ++i) {
      istringstream s(sol_lines.at(sol_lines.size()-ng_-nx_+i-1));
      s >> m->lam_g[i];
      m->lam_g[i] *= -1;
    }

    // Close and delete .sol file
    sol.close();
    if (remove(solname.c_str())!=0) {
      casadi_warning("Failed to remove " + solname);
    }
    if (verbose_) casadi_message("Removed " + solname);

    return 0;
  }


} // namespace casadi
