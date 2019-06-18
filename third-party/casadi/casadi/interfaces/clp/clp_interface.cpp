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

#include "clp_interface.hpp"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_CLP_EXPORT
  casadi_register_conic_clp(Conic::Plugin* plugin) {
    plugin->creator = ClpInterface::creator;
    plugin->name = "clp";
    plugin->doc = ClpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ClpInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_CLP_EXPORT casadi_load_conic_clp() {
    Conic::registerPlugin(casadi_register_conic_clp);
  }

  ClpInterface::ClpInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }


  inline std::string return_status_string(int status) {
    switch (status) {
    case 0:
      return "optimal";
    case 1:
      return "primal infeasible";
    case 2:
      return "dual infeasible";
    case 3:
      return "stopped on iterations or time";
    case 4:
      return "stopped due to errors";
    case 5:
      return "stopped by event handler";
    default:
      return "unknown";
    }
  }


  inline std::string return_secondary_status_string(int status) {
    switch (status) {
    case 0:
      return "none";
    case 1:
      return "primal infeasible because dual limit reached OR (probably primal"
         " infeasible but can't prove it  - main status was 4)";
    case 2:
      return "scaled problem optimal - unscaled problem has primal infeasibilities";
    case 3:
      return "scaled problem optimal - unscaled problem has dual infeasibilities";
    case 4:
      return "scaled problem optimal - unscaled problem has primal and dual infeasibilities";
    case 5:
      return "giving up in primal with flagged variables";
    case 6:
      return "failed due to empty problem check";
    case 7:
      return "postSolve says not optimal";
    case 8:
      return "failed due to bad element check";
    case 9:
      return "status was 3 and stopped on time";
    case 10:
      return "status was 3 but stopped as primal feasibles";
    case ClpEventHandler::Event::endOfIteration:
      return "endOfIteration";
    case ClpEventHandler::Event::endOfFactorization:
      return "endOfFactorization";
    case ClpEventHandler::Event::endOfValuesPass:
      return "endOfValuesPass";
    case ClpEventHandler::Event::node:
      return "node";
    case ClpEventHandler::Event::treeStatus:
      return "treeStatus";
    case ClpEventHandler::Event::solution:
      return "solution";
    case ClpEventHandler::Event::theta:
      return "theta";
    case ClpEventHandler::Event::pivotRow:
      return "pivotRow";
    case ClpEventHandler::Event::presolveStart:
      return "presolveStart";
    case ClpEventHandler::Event::presolveSize:
      return "presolveSize";
    case ClpEventHandler::Event::presolveInfeasible:
      return "presolveInfeasible";
    case ClpEventHandler::Event::presolveBeforeSolve:
      return "presolveBeforeSolve";
    case ClpEventHandler::Event::presolveAfterFirstSolve:
      return "presolveAfterFirstSolve";
    case ClpEventHandler::Event::presolveAfterSolve:
      return "presolveAfterSolve";
    case ClpEventHandler::Event::presolveEnd:
      return "presolveEnd";
    case ClpEventHandler::Event::goodFactorization:
      return "goodFactorization";
    case ClpEventHandler::Event::complicatedPivotIn:
      return "complicatedPivotIn";
    case ClpEventHandler::Event::noCandidateInPrimal:
      return "noCandidateInPrimal";
    case ClpEventHandler::Event::looksEndInPrimal:
      return "looksEndInPrimal";
    case ClpEventHandler::Event::endInPrimal:
      return "endInPrimal";
    case ClpEventHandler::Event::beforeStatusOfProblemInPrimal:
      return "beforeStatusOfProblemInPrimal";
    case ClpEventHandler::Event::startOfStatusOfProblemInPrimal:
      return "startOfStatusOfProblemInPrimal";
    case ClpEventHandler::Event::complicatedPivotOut:
      return "complicatedPivotOut";
    case ClpEventHandler::Event::noCandidateInDual:
      return "noCandidateInDual";
    case ClpEventHandler::Event::looksEndInDual:
      return "looksEndInDual";
    case ClpEventHandler::Event::endInDual:
      return "endInDual";
    case ClpEventHandler::Event::beforeStatusOfProblemInDual:
      return "beforeStatusOfProblemInDual";
    case ClpEventHandler::Event::startOfStatusOfProblemInDual:
      return "startOfStatusOfProblemInDual";
    case ClpEventHandler::Event::startOfIterationInDual:
      return "startOfIterationInDual";
    case ClpEventHandler::Event::updateDualsInDual:
      return "updateDualsInDual";
    case ClpEventHandler::Event::endOfCreateRim:
      return "endOfCreateRim";
    case ClpEventHandler::Event::slightlyInfeasible:
      return "slightlyInfeasible";
    case ClpEventHandler::Event::modifyMatrixInMiniPresolve:
      return "modifyMatrixInMiniPresolve";
    case ClpEventHandler::Event::moreMiniPresolve:
      return "moreMiniPresolve";
    case ClpEventHandler::Event::modifyMatrixInMiniPostsolve:
      return "modifyMatrixInMiniPostsolve";
    case ClpEventHandler::Event::startOfCrossover:
      return "startOfCrossover";
    case ClpEventHandler::Event::noTheta:
      return "noTheta";
    default:
      return "unknown";
    }
  }

  class CasadiHandler : public CoinMessageHandler {
    public:
      virtual int print() ;
  };

  int CasadiHandler::print() {
    uout() << messageBuffer() << std::endl;
    return 0;
  }

  void ClpInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Default options
    casadi_assert(H_.nnz()==0, "Not an LP");

    // Read options
    //for (auto&& op : opts) {
    //}

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnz_in(CONIC_H), true); // H
    alloc_w(nnz_in(CONIC_A), true); // A
  }

  int ClpInterface::init_mem(void* mem) const {
    if (!mem) return 1;
    auto m = static_cast<ClpMemory*>(mem);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();

    m->colind.resize(A_.size2()+1);
    m->row.resize(A_.nnz());

    // Problem has not been solved at this point
    m->success = false;
    return 0;
  }

  int ClpInterface::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<ClpMemory*>(mem);

    // Problem has not been solved at this point
    m->success = false;
    m->return_status = -1;
    m->secondary_return_status = -1;

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();
    m->fstats.at("preprocessing").tic();

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

    copy_vector(A_.colind(), m->colind);
    copy_vector(A_.row(), m->row);

    // Create model
    ClpSimplex model;
    model.loadProblem(A_.size2(), A_.size1(), get_ptr(m->colind), get_ptr(m->row), A,
                      lbx, ubx, g, lba, uba, nullptr);

    CasadiHandler ch;
    model.passInMessageHandler(&ch);

    m->fstats.at("preprocessing").toc();
    m->fstats.at("solver").tic();

    // Solve the problem using the primal simplex algorithm
    model.primal();

    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();

    // Primal solution
    double* x = model.primalColumnSolution();
    casadi_copy(x, nx_, res[CONIC_X]);

    // Dual solution (x)
    double* minus_lam_x = model.dualColumnSolution();
    if (res[CONIC_LAM_X]) {
      casadi_copy(minus_lam_x, nx_, res[CONIC_LAM_X]);
      casadi_scal(nx_, -1., res[CONIC_LAM_X]);
    }

    // Dual solution (A)
    double* minus_lam_a = model.dualRowSolution();
    if (res[CONIC_LAM_A]) {
      casadi_copy(minus_lam_a, na_, res[CONIC_LAM_A]);
      casadi_scal(na_, -1., res[CONIC_LAM_A]);
    }

    // Optimal cost
    double f = model.rawObjectiveValue();
    if (res[CONIC_COST]) *res[CONIC_COST] = f;

    m->fstats.at("postprocessing").toc();

    // Show statistics
    if (print_time_)  print_fstats(static_cast<ConicMemory*>(mem));

    m->return_status = model.status();
    m->success = m->return_status==0;
    m->secondary_return_status = model.secondaryStatus();

    if (verbose_) casadi_message("CLP return status: " + return_status_string(m->return_status));
    if (verbose_) casadi_message(
      "CLP secondary return status: " + return_secondary_status_string(m->secondary_return_status));

    return 0;
  }

  ClpInterface::~ClpInterface() {
    clear_mem();
  }

  ClpMemory::ClpMemory() {
  }

  ClpMemory::~ClpMemory() {
  }


  Dict ClpInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<ClpMemory*>(mem);
    stats["return_status"] = return_status_string(m->return_status);
    stats["secondary_return_status"] = return_secondary_status_string(m->secondary_return_status);
    stats["success"] = m->success;
    return stats;
  }

} // end namespace casadi
