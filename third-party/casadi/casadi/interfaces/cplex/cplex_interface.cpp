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

#include "cplex_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_CPLEX_EXPORT
  casadi_register_conic_cplex(Conic::Plugin* plugin) {
    plugin->creator = CplexInterface::creator;
    plugin->name = "cplex";
    plugin->doc = CplexInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &CplexInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_CPLEX_EXPORT casadi_load_conic_cplex() {
    Conic::registerPlugin(casadi_register_conic_cplex);
  }

  CplexInterface::CplexInterface(const std::string& name,
                                 const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  Options CplexInterface::options_
  = {{&Conic::options_},
     {{"cplex",
       {OT_DICT,
        "Options to be passed to CPLEX"}},
      {"qp_method",
       {OT_INT,
        "Determines which CPLEX algorithm to use."}},
      {"dump_to_file",
       {OT_BOOL,
        "Dumps QP to file in CPLEX format."}},
      {"dump_filename",
       {OT_STRING,
        "The filename to dump to."}},
      {"tol",
       {OT_DOUBLE,
        "Tolerance of solver"}},
      {"dep_check",
       {OT_INT,
        "Detect redundant constraints."}},
      {"warm_start",
       {OT_BOOL,
        "Use warm start with simplex methods (affects only the simplex methods)."}}
     }
  };

  void CplexInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Default options
    qp_method_ = 0;
    dump_to_file_ = false;
    dump_filename_ = "qp.dat";
    tol_ = 1e-6;
    dep_check_ = 0;
    warm_start_ = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="cplex") {
        opts_ = op.second;
      } else if (op.first=="qp_method") {
        qp_method_ = op.second;
      } else if (op.first=="dump_to_file") {
        dump_to_file_ = op.second;
      } else if (op.first=="dump_filename") {
        dump_filename_ = op.second.to_string();
      } else if (op.first=="tol") {
        tol_ = op.second;
      } else if (op.first=="dep_check") {
        dep_check_ = op.second;
      } else if (op.first=="warm_start") {
        warm_start_ = op.second;
      }
    }

    // Are we solving a mixed-integer problem?
    mip_ = !discrete_.empty()
      && find(discrete_.begin(), discrete_.end(), true)!=discrete_.end();

    // Type of variable
    if (mip_) {
      ctype_.resize(nx_);
      for (casadi_int i=0; i<nx_; ++i) {
        ctype_[i] = discrete_[i] ? 'I' : 'C';
      }
    }

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnz_in(CONIC_H), true); // H
    alloc_w(nnz_in(CONIC_A), true); // A
    alloc_w(nx_, true); // x
    alloc_w(nx_, true); // lam_x
    alloc_w(na_, true); // lam_a
  }

  int CplexInterface::init_mem(void* mem) const {
    if (!mem) return 1;
    auto m = static_cast<CplexMemory*>(mem);

    // Start CPLEX
    int status;
    casadi_assert_dev(m->env==nullptr);
    m->env = CPXXopenCPLEX(&status);
    if (m->env==nullptr) {
      char errmsg[CPXMESSAGEBUFSIZE];
      CPXXgeterrorstring(m->env, status, errmsg);
      casadi_error(string("Cannot initialize CPLEX environment: ") + errmsg);
    }

    // Set parameters to their default values
    if (CPXXsetdefaults(m->env)) {
      casadi_error("CPXXsetdefaults failed");
    }

    // Enable output by default
    if (CPXXsetintparam(m->env, CPX_PARAM_SCRIND, CPX_ON)) {
      casadi_error("Failure setting CPX_PARAM_SCRIND");
    }

    // Optimality tolerance
    if (CPXXsetdblparam(m->env, CPX_PARAM_EPOPT, tol_)) {
      casadi_error("Failure setting CPX_PARAM_EPOPT");
    }

    // Feasibility tolerance
    if (CPXXsetdblparam(m->env, CPX_PARAM_EPRHS, tol_)) {
      casadi_error("Failure setting CPX_PARAM_EPRHS");
    }

    // We start with barrier if crossover was chosen.
    if (CPXXsetintparam(m->env, CPX_PARAM_QPMETHOD, qp_method_ == 7 ? 4 : qp_method_)) {
      casadi_error("Failure setting CPX_PARAM_QPMETHOD");
    }

    // Setting dependency check option
    if (CPXXsetintparam(m->env, CPX_PARAM_DEPIND, dep_check_)) {
      casadi_error("Failure setting CPX_PARAM_DEPIND");
    }

    // Setting crossover algorithm
    if (qp_method_ == 7) {
      if (CPXXsetintparam(m->env, CPX_PARAM_BARCROSSALG, 1)) {
        casadi_error("Failure setting CPX_PARAM_BARCROSSALG");
      }
    }

    // Set parameters
    for (auto&& op : opts_) {
      // Get parameter index
      int whichparam;
      if (CPXXgetparamnum(m->env, op.first.c_str(), &whichparam)) {
        casadi_error("No such CPLEX parameter: " + op.first);
      }

      // Get type of parameter
      int paramtype;
      if (CPXXgetparamtype(m->env, whichparam, &paramtype)) {
        casadi_error("CPXXgetparamtype failed");
      }

      // Pass to CPLEX
      switch (paramtype) {
      case CPX_PARAMTYPE_NONE:
        casadi_error("CPX_PARAMTYPE_NONE unsupported");
        break;
      case CPX_PARAMTYPE_INT:
        status = CPXXsetintparam(m->env, whichparam, op.second);
        break;
      case CPX_PARAMTYPE_DOUBLE:
        status = CPXXsetdblparam(m->env, whichparam, op.second);
        break;
      case CPX_PARAMTYPE_STRING:
        status = CPXXsetstrparam(m->env, whichparam,
                                static_cast<string>(op.second).c_str());
        break;
      case CPX_PARAMTYPE_LONG:
        status = CPXXsetlongparam(m->env, whichparam,
                                 static_cast<CPXLONG>(static_cast<casadi_int>(op.second)));
        break;
        default:
          casadi_error("Unknown CPLEX parameter type (" + str(paramtype) + ") for " + op.first);
      }
      // Error handling
      if (status) {
        casadi_error("Failure setting option " + op.first);
      }
    }

    // Doing allocation of CPLEX data
    // Objective is to be minimized
    m->objsen = CPX_MIN;

    // Allocation of data
    // Type of constraint
    m->sense.resize(na_);
    // Right-hand side of constraints
    m->rhs.resize(na_);
    // Range value for lower AND  upper bounded constraints
    m->rngval.resize(na_);
    // Basis for primal variables
    m->cstat.resize(nx_);
    m->rstat.resize(na_);

    // Matrix A, count the number of elements per column
    m->matcnt.resize(A_.size2());
    transform(A_.colind()+1, A_.colind() + A_.size2()+1, A_.colind(), m->matcnt.begin(),
              std::minus<casadi_int>());

    // Matrix H, count the number of elements per column
    m->qmatcnt.resize(H_.size2());
    transform(H_.colind()+1, H_.colind() + H_.size2()+1, H_.colind(), m->qmatcnt.begin(),
              std::minus<casadi_int>());

    // Create problem object
    casadi_assert_dev(m->lp==nullptr);
    m->lp = CPXXcreateprob(m->env, &status, "QP from CasADi");
    casadi_assert(m->lp!=nullptr, "CPXXcreateprob failed");

    m->a_colind.resize(A_.size2()+1);
    m->a_row.resize(A_.nnz());
    m->h_colind.resize(H_.size2()+1);
    m->h_row.resize(H_.nnz());

    copy_vector(A_.colind(), m->a_colind);
    copy_vector(A_.row(), m->a_row);
    copy_vector(H_.colind(), m->h_colind);
    copy_vector(H_.row(), m->h_row);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }


  inline std::string return_status_string(int status) {
    switch (status) {
    case CPX_STAT_OPTIMAL:
      return "Optimal solution found";
    case CPX_STAT_UNBOUNDED:
      return "Model is unbounded";
    case CPX_STAT_INForUNBD:
      return "Model is infeasible or unbounded";
    case CPX_STAT_OPTIMAL_INFEAS:
      return "Optimal solution is available but with infeasibilities";
    case CPX_STAT_NUM_BEST:
      return "Solution available, but not proved optimal due to numeric difficulties";
    case CPX_STAT_FIRSTORDER:
      return "Solution satisfies first-order optimality conditions, "
             "but is not necessarily globally optimal";
    default:
      return "unknown";
    }
  }

  int CplexInterface::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<CplexMemory*>(mem);

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();

    m->fstats.at("preprocessing").tic();

    // Problem has not been solved at this point
    m->success = false;
    m->return_status = -1;

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Get inputs
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lbx=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lbx);
    double* ubx=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ubx);
    double* lba=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lba);
    double* uba=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, uba);
    double* H=w; w += nnz_in(CONIC_H);
    casadi_copy(arg[CONIC_H], nnz_in(CONIC_H), H);
    double* A=w; w += nnz_in(CONIC_A);
    casadi_copy(arg[CONIC_A], nnz_in(CONIC_A), A);
    double* x=w; w += nx_;
    casadi_copy(arg[CONIC_X0], nx_, x);
    double* lam_x=w; w += nx_;
    casadi_copy(arg[CONIC_LAM_X0], nx_, lam_x);

    // Temporaries
    double* lam_a=w; w += na_;

    // We change method in crossover
    if (m->is_warm && qp_method_ == 7) {
      (void)CPXXsetintparam(m->env, CPX_PARAM_QPMETHOD, 1);
    }

    for (casadi_int i = 0; i < na_; ++i) {
      // CPX_INFBOUND

      // Equality
      if (uba[i] - lba[i] < 1e-20) {
        m->sense[i] = 'E';
        m->rhs[i] = lba[i];
        m->rngval[i] = 0.;
      } else if (lba[i] < -CPX_INFBOUND) {
        // Ineq - no lower bound
        m->sense[i] = 'L';
        m->rhs[i] = uba[i];
        m->rngval[i] = 0.;
      } else if (uba[i] > CPX_INFBOUND) {
        // Ineq - no upper bound
        m->sense[i] = 'G';
        m->rhs[i] = lba[i];
        m->rngval[i] = 0.;
      } else { // Inew both upper and lower bounds
        m->sense[i] = 'R';
        m->rhs[i] = lba[i];
        m->rngval[i] = uba[i] - lba[i];
      }
    }

    // Copying objective, constraints, and bounds.
    const CPXNNZ* matbeg = get_ptr(m->a_colind);
    const CPXDIM* matind = get_ptr(m->a_row);

    const double* matval = A;
    const double* obj = g;
    const double* lb = lbx;
    const double* ub = ubx;
    if (CPXXcopylp(m->env, m->lp, nx_, na_, m->objsen, obj, get_ptr(m->rhs), get_ptr(m->sense),
                  matbeg, get_ptr(m->matcnt), matind, matval, lb, ub, get_ptr(m->rngval))) {
      casadi_error("CPXXcopylp failed");
    }

    // Preparing coefficient matrix Q
    const CPXNNZ* qmatbeg = get_ptr(m->h_colind);
    const CPXDIM* qmatind = get_ptr(m->h_row);
    const double* qmatval = H;
    if (CPXXcopyquad(m->env, m->lp, qmatbeg, get_ptr(m->qmatcnt), qmatind, qmatval)) {
    }

    if (dump_to_file_) {
      CPXXwriteprob(m->env, m->lp, dump_filename_.c_str(), "LP");
      casadi_error("CPXXwriteprob failed");
    }

    // Warm-starting if possible
    if (qp_method_ != 0 && qp_method_ != 4 && m->is_warm) {
      // TODO(Joel): Initialize slacks and dual variables of bound constraints
      if (CPXXcopystart(m->env, m->lp, get_ptr(m->cstat), get_ptr(m->rstat), x,
           nullptr, nullptr, lam_x)) {
        casadi_error("CPXXcopystart failed");
      }
    } else {
      if (CPXXcopystart(m->env, m->lp, nullptr, nullptr, x,
           nullptr, nullptr, lam_x)) {
        casadi_error("CPXXcopystart failed");
      }
    }

    // Solution
    double f;
    std::vector<double> slack(na_);
    int solstat;

    if (mip_) {
      // Pass type of variables
      if (CPXXcopyctype(m->env, m->lp, &ctype_[0])) {
        casadi_error("CPXXcopyctype failed");
      }

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();
      // Optimize
      if (CPXXmipopt(m->env, m->lp)) {
        casadi_error("CPXXmipopt failed");
      }
      m->fstats.at("solver").toc();
      m->fstats.at("postprocessing").tic();

      // Get objective value
      if (CPXXgetobjval(m->env, m->lp, &f)) {
        casadi_error("CPXXgetobjval failed");
      }

      // Get primal solution
      casadi_int cur_numcols = CPXXgetnumcols(m->env, m->lp);
      if (CPXXgetx(m->env, m->lp, x, 0, cur_numcols-1)) {
        casadi_error("CPXXgetx failed");
      }

      // Get slacks
      casadi_int cur_numrows = CPXXgetnumrows(m->env, m->lp);
      if (CPXXgetslack(m->env, m->lp, get_ptr(slack), 0, cur_numrows-1)) {
        casadi_error("CPXXgetslack failed");
      }

      // Not a number as dual variables (not calculated with MIQP algorithm)
      casadi_fill(lam_a, na_, nan);
      casadi_fill(lam_x, nx_, nan);

    } else {
      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();
      // Optimize
      if (CPXXqpopt(m->env, m->lp)) {
        casadi_error("CPXXqpopt failed");
      }
      m->fstats.at("solver").toc();
      m->fstats.at("postprocessing").tic();

      casadi_int problem_type = CPXXgetprobtype(m->env, m->lp);

      // The solver switched to a MIQP
      if (problem_type == CPXPROB_MIQP) {

          // Get objective value
          if (CPXXgetobjval(m->env, m->lp, &f)) {
            casadi_error("CPXXgetobjval failed");
          }

          // Get primal solution
          casadi_int cur_numcols = CPXXgetnumcols(m->env, m->lp);
          if (CPXXgetx(m->env, m->lp, x, 0, cur_numcols-1)) {
            casadi_error("CPXXgetx failed");
          }

          // Not a number as dual variables (not calculated with MIQP algorithm)
          casadi_fill(lam_a, na_, nan);
          casadi_fill(lam_x, nx_, nan);
          if (verbose_) casadi_message("CPLEX does not compute dual variables for nonconvex QPs");

      } else {

          // Retrieving solution
          if (CPXXsolution(m->env, m->lp, &solstat, &f, x, lam_a, get_ptr(slack), lam_x)) {
            casadi_error("CPXXsolution failed");
          }
      }
    }

    // Retrieving the basis
    if (qp_method_ != 0 && qp_method_ != 4) {
      (void)CPXXgetbase(m->env, m->lp, get_ptr(m->cstat), get_ptr(m->rstat));
    }

    // Flip the sign of the multipliers
    casadi_scal(na_, -1., lam_a);
    casadi_scal(nx_, -1., lam_x);

    m->return_status = CPXXgetstat(m->env, m->lp);
    m->success = m->return_status==CPX_STAT_OPTIMAL || m->return_status==CPX_STAT_FIRSTORDER;

    if (verbose_) casadi_message("CPLEX return status: " + return_status_string(m->return_status));

    // Next time we warm start
    if (warm_start_) {
      m->is_warm = true;
    }

    // Get the outputs
    if (res[CONIC_COST]) *res[CONIC_COST] = f;
    casadi_copy(lam_a, na_, res[CONIC_LAM_A]);
    casadi_copy(lam_x, nx_, res[CONIC_LAM_X]);
    casadi_copy(x, nx_, res[CONIC_X]);

    m->fstats.at("postprocessing").toc();

    // Show statistics
    if (print_time_)  print_fstats(static_cast<ConicMemory*>(mem));
    return 0;
  }

  CplexInterface::~CplexInterface() {
    clear_mem();
  }

  Dict CplexInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<CplexMemory*>(mem);
    stats["return_status"] = return_status_string(m->return_status);
    stats["success"] = m->success;
    return stats;
  }


  CplexMemory::CplexMemory() {
    // Setting warm-start flag
    this->is_warm = false;

    // Set pointer to zero to avoid deleting a nonexisting instance
    this->env = nullptr;
    this->lp = nullptr;
  }

  CplexMemory::~CplexMemory() {
    // Return flag
    casadi_int status;

    // Only free if Cplex problem if it has been allocated
    if (this->lp!=nullptr) {
      status = CPXXfreeprob(this->env, &this->lp);
      if (status!=0) {
        uerr() << "CPXXfreeprob failed, error code " << status << ".\n";
      }
      this->lp = nullptr;
    }

    // Closing down license
    if (this->env!=nullptr) {
      status = CPXXcloseCPLEX(&this->env);
      if (status!=0) {
        uerr() << "CPXXcloseCPLEX failed, error code " << status << ".\n";
      }
      this->env = nullptr;
    }
  }



} // end namespace casadi
