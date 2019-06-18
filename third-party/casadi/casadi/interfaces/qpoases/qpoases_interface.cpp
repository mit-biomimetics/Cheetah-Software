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


#include "qpoases_interface.hpp"

// Bug in qpOASES?
#define ALLOW_QPROBLEMB true
#define ALLOW_ALL_OPTIONS

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QPOASES_EXPORT
  casadi_register_conic_qpoases(Conic::Plugin* plugin) {
    plugin->creator = QpoasesInterface::creator;
    plugin->name = "qpoases";
    plugin->doc = QpoasesInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &QpoasesInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QPOASES_EXPORT casadi_load_conic_qpoases() {
    Conic::registerPlugin(casadi_register_conic_qpoases);
  }

  QpoasesInterface::QpoasesInterface(const std::string& name,
                                     const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
    // Redirect output to CasADi
    static bool first_call = true;
    if (first_call) {
      qpOASES::setPrintf(qpoases_printf);
      first_call = false;
    }
  }

  void QpoasesInterface::qpoases_printf(const char* s) {
    uout() << s;
  }

  QpoasesInterface::~QpoasesInterface() {
    clear_mem();
  }

  Options QpoasesInterface::options_
  = {{&Conic::options_},
     {{"sparse",
       {OT_BOOL,
        "Formulate the QP using sparse matrices. [false]"}},
      {"schur",
       {OT_BOOL,
        "Use Schur Complement Approach [false]"}},
      {"hessian_type",
       {OT_STRING,
        "Type of Hessian - see qpOASES documentation "
        "[UNKNOWN|posdef|semidef|indef|zero|identity]]"}},
      {"max_schur",
       {OT_INT,
        "Maximal number of Schur updates [75]"}},
      {"linsol_plugin",
       {OT_STRING,
        "Linear solver plugin"}},
      {"nWSR",
       {OT_INT,
        "The maximum number of working set recalculations to be performed during "
        "the initial homotopy. Default is 5(nx + nc)"}},
      {"CPUtime",
       {OT_DOUBLE,
        "The maximum allowed CPU time in seconds for the whole initialisation"
        " (and the actually required one on output). Disabled if unset."}},
      {"printLevel",
       {OT_STRING,
        "Defines the amount of text output during QP solution, see Section 5.7"}},
      {"enableRamping",
       {OT_BOOL,
        "Enables ramping."}},
      {"enableFarBounds",
       {OT_BOOL,
        "Enables the use of  far bounds."}},
      {"enableFlippingBounds",
       {OT_BOOL,
        "Enables the use of  flipping bounds."}},
      {"enableRegularisation",
       {OT_BOOL,
        "Enables automatic  Hessian regularisation."}},
      {"enableFullLITests",
       {OT_BOOL,
        "Enables condition-hardened  (but more expensive) LI test."}},
      {"enableNZCTests",
       {OT_BOOL,
        "Enables nonzero curvature  tests."}},
      {"enableDriftCorrection",
       {OT_INT,
        "Specifies the frequency of drift corrections: 0: turns them off."}},
      {"enableCholeskyRefactorisation",
       {OT_INT,
        "Specifies the frequency of a full re-factorisation of projected "
        "Hessian matrix: 0: turns them off,  1: uses them at each iteration etc."}},
      {"enableEqualities",
       {OT_BOOL,
        "Specifies whether equalities should be treated  as always active "
        "(True) or not (False)"}},
      {"terminationTolerance",
       {OT_DOUBLE,
        "Relative termination tolerance to stop homotopy."}},
      {"boundTolerance",
       {OT_DOUBLE,
        "If upper and lower bounds differ less than this tolerance, they are regarded "
        "equal, i.e. as  equality constraint."}},
      {"boundRelaxation",
       {OT_DOUBLE,
        "Initial relaxation of bounds to start homotopy  and initial value for far bounds."}},
      {"epsNum",
       {OT_DOUBLE,
        "Numerator tolerance for ratio tests."}},
      {"epsDen",
       {OT_DOUBLE,
        "Denominator tolerance for ratio tests."}},
      {"maxPrimalJump",
       {OT_DOUBLE,
        "Maximum allowed jump in primal variables in  nonzero curvature tests."}},
      {"maxDualJump",
       {OT_DOUBLE,
        "Maximum allowed jump in dual variables in  linear independence tests."}},
      {"initialRamping",
       {OT_DOUBLE,
        "Start value for ramping strategy."}},
      {"finalRamping",
       {OT_DOUBLE,
        "Final value for ramping strategy."}},
      {"initialFarBounds",
       {OT_DOUBLE,
        "Initial size for far bounds."}},
      {"growFarBounds",
       {OT_DOUBLE,
        "Factor to grow far bounds."}},
      {"initialStatusBounds",
       {OT_STRING,
        "Initial status of bounds at first iteration."}},
      {"epsFlipping",
       {OT_DOUBLE,
        "Tolerance of squared Cholesky diagonal factor  which triggers flipping bound."}},
      {"numRegularisationSteps",
       {OT_INT,
        "Maximum number of successive regularisation steps."}},
      {"epsRegularisation",
       {OT_DOUBLE,
        "Scaling factor of identity matrix used for  Hessian regularisation."}},
      {"numRefinementSteps",
       {OT_INT,
        "Maximum number of iterative refinement steps."}},
      {"epsIterRef",
       {OT_DOUBLE,
        "Early termination tolerance for iterative  refinement."}},
      {"epsLITests",
       {OT_DOUBLE,
        "Tolerance for linear independence tests."}},
      {"epsNZCTests",
       {OT_DOUBLE,
        "Tolerance for nonzero curvature tests."}},
      {"enableInertiaCorrection",
       {OT_BOOL,
        "Should working set be repaired when negative curvature is discovered during hotstart."}}
     }
  };

  void QpoasesInterface::init(const Dict& opts) {
    Conic::init(opts);

    // Default options
    sparse_ = false;
    schur_ = false;
    hess_ = qpOASES::HessianType::HST_UNKNOWN;
    max_schur_ = 75;
    max_nWSR_ = 5 *(nx_ + na_);
    max_cputime_ = -1;
    ops_.setToDefault();
    linsol_plugin_ = "ma27";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="sparse") {
        sparse_ = op.second;
      } else if (op.first=="schur") {
        schur_=  op.second;
      } else if (op.first=="hessian_type") {
        string h = op.second;
        if (h=="unknown") {
          hess_ = qpOASES::HessianType::HST_UNKNOWN;
        } else if (h=="posdef") {
          hess_ = qpOASES::HessianType::HST_POSDEF;
        } else if (h=="semidef") {
          hess_ = qpOASES::HessianType::HST_SEMIDEF;
        } else if (h=="indef") {
          hess_ = qpOASES::HessianType::HST_INDEF;
        } else if (h=="zero") {
          hess_ = qpOASES::HessianType::HST_ZERO;
        } else if (h=="identity") {
          hess_ = qpOASES::HessianType::HST_IDENTITY;
        } else {
          casadi_error("Unknown Hessian type \"" + h + "\"");
        }
      } else if (op.first=="max_schur") {
        max_schur_ = op.second;
      } else if (op.first=="linsol_plugin") {
        linsol_plugin_ = string(op.second);
      } else if (op.first=="nWSR") {
        max_nWSR_ = op.second;
      } else if (op.first=="CPUtime") {
        max_cputime_ = op.second;
      } else if (op.first=="printLevel") {
        ops_.printLevel = to_PrintLevel(op.second);
      } else if (op.first=="enableRamping") {
        ops_.enableRamping = to_BooleanType(op.second);
      } else if (op.first=="enableFarBounds") {
        ops_.enableFarBounds = to_BooleanType(op.second);
      } else if (op.first=="enableFlippingBounds") {
        ops_.enableFlippingBounds = to_BooleanType(op.second);
      } else if (op.first=="enableRegularisation") {
        ops_.enableRegularisation = to_BooleanType(op.second);
      } else if (op.first=="enableFullLITests") {
        ops_.enableFullLITests = to_BooleanType(op.second);
      } else if (op.first=="enableNZCTests") {
        ops_.enableNZCTests = to_BooleanType(op.second);
      } else if (op.first=="enableDriftCorrection") {
        ops_.enableRegularisation = to_BooleanType(op.second);
      } else if (op.first=="enableCholeskyRefactorisation") {
        ops_.enableCholeskyRefactorisation = op.second;
      } else if (op.first=="enableEqualities") {
        ops_.enableEqualities = to_BooleanType(op.second);
      } else if (op.first=="terminationTolerance") {
        ops_.terminationTolerance = op.second;
      } else if (op.first=="boundTolerance") {
        ops_.boundTolerance = op.second;
      } else if (op.first=="boundRelaxation") {
        ops_.boundRelaxation = op.second;
      } else if (op.first=="epsNum") {
        ops_.epsNum = op.second;
      } else if (op.first=="epsDen") {
        ops_.epsDen = op.second;
      } else if (op.first=="maxPrimalJump") {
        ops_.maxPrimalJump = op.second;
      } else if (op.first=="maxDualJump") {
        ops_.maxDualJump = op.second;
      } else if (op.first=="initialRamping") {
        ops_.initialRamping = op.second;
      } else if (op.first=="finalRamping") {
        ops_.finalRamping = op.second;
      } else if (op.first=="initialFarBounds") {
        ops_.initialFarBounds = op.second;
      } else if (op.first=="growFarBounds") {
        ops_.growFarBounds = op.second;
      } else if (op.first=="initialStatusBounds") {
        ops_.initialStatusBounds = to_SubjectToStatus(op.second);
      } else if (op.first=="epsFlipping") {
        ops_.epsFlipping = op.second;
      } else if (op.first=="numRegularisationSteps") {
        ops_.numRegularisationSteps = op.second;
      } else if (op.first=="epsRegularisation") {
        ops_.epsRegularisation = op.second;
      } else if (op.first=="numRefinementSteps") {
        ops_.numRefinementSteps = op.second;
      } else if (op.first=="epsIterRef") {
        ops_.epsIterRef = op.second;
      } else if (op.first=="epsLITests") {
        ops_.epsLITests = op.second;
      } else if (op.first=="epsNZCTests") {
        ops_.epsNZCTests = op.second;
      } else if (op.first=="enableInertiaCorrection") {
        ops_.enableInertiaCorrection = to_BooleanType(op.second);
      }
    }

    // Allocate work vectors
    if (sparse_) {
      alloc_w(nnz_in(CONIC_H), true); // h
      alloc_w(nnz_in(CONIC_A), true); // a
    } else {
      alloc_w(nx_*nx_, true); // h
      alloc_w(nx_*na_, true); // a
    }
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nx_+na_, true); // dual
  }

  int QpoasesInterface::init_mem(void* mem) const {
    auto m = static_cast<QpoasesMemory*>(mem);
    m->called_once = false;

    // Linear solver, if any
    m->linsol_plugin = linsol_plugin_;

    // Create qpOASES instance
    if (m->qp) delete m->qp;
    if (schur_) {
      m->sqp = new qpOASES::SQProblemSchur(nx_, na_, hess_, max_schur_,
        mem, qpoases_init, qpoases_sfact, qpoases_nfact, qpoases_solve);
    } else if (na_==0) {
      m->qp = new qpOASES::QProblemB(nx_, hess_);
    } else {
      m->sqp = new qpOASES::SQProblem(nx_, na_, hess_);
    }

    // Pass to qpOASES
    m->qp->setOptions(ops_);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    m->h_row.resize(H_.nnz());
    m->h_colind.resize(H_.size2()+1);
    m->a_row.resize(A_.nnz());
    m->a_colind.resize(A_.size2()+1);

    return 0;
  }

  int QpoasesInterface::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QpoasesMemory*>(mem);

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();

    m->fstats.at("preprocessing").tic();

    // Problem has not been solved at this point
    m->success = false;
    m->return_status = -1;

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Maxiumum number of working set changes
    int nWSR = max_nWSR_;
    double cputime = max_cputime_;
    double *cputime_ptr = cputime<=0 ? nullptr : &cputime;

    // Get the arguments to call qpOASES with
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lb=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lb);
    double* ub=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ub);
    double* lbA=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lbA);
    double* ubA=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, ubA);

    // Return flag
    casadi_int flag;

    // Sparse or dense mode?
    if (sparse_) {
      // Get quadratic term
      copy_vector(H_.colind(), m->h_colind);
      copy_vector(H_.row(), m->h_row);
      double* h=w; w += H_.nnz();
      casadi_copy(arg[CONIC_H], H_.nnz(), h);
      if (m->h) delete m->h;
      m->h = new qpOASES::SymSparseMat(H_.size1(), H_.size2(),
        get_ptr(m->h_row), get_ptr(m->h_colind), h);
      m->h->createDiagInfo();

      // Get linear term
      copy_vector(A_.colind(), m->a_colind);
      copy_vector(A_.row(), m->a_row);
      double* a=w; w += A_.nnz();
      casadi_copy(arg[CONIC_A], A_.nnz(), a);
      if (m->a) delete m->a;
      m->a = new qpOASES::SparseMatrix(A_.size1(), A_.size2(),
        get_ptr(m->a_row), get_ptr(m->a_colind), a);

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();

      // Solve sparse
      if (m->called_once) {
        flag = m->sqp->hotstart(m->h, g, m->a, lb, ub, lbA, ubA, nWSR, cputime_ptr);
      } else {
        flag = m->sqp->init(m->h, g, m->a, lb, ub, lbA, ubA, nWSR, cputime_ptr);
      }
      m->fstats.at("solver").toc();

    } else {
      // Get quadratic term
      double* h=w; w += nx_*nx_;
      casadi_densify(arg[CONIC_H], H_, h, false);

      // Get linear term
      double* a = w; w += nx_*na_;
      casadi_densify(arg[CONIC_A], A_, a, true);

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();
      // Solve dense
      if (na_==0) {
        if (m->called_once) {
          // Broken?
          //flag = m->qp->hotstart(g, lb, ub, nWSR, cputime_ptr);
          m->qp->reset();
          flag = m->qp->init(h, g, lb, ub, nWSR, cputime_ptr);
        } else {
          flag = m->qp->init(h, g, lb, ub, nWSR, cputime_ptr);
        }
      } else {
        if (m->called_once) {
          flag = m->sqp->hotstart(h, g, a, lb, ub, lbA, ubA, nWSR, cputime_ptr);
        } else {
          flag = m->sqp->init(h, g, a, lb, ub, lbA, ubA, nWSR, cputime_ptr);
        }
      }
      m->fstats.at("solver").toc();
    }

    // Solver is "warm" now
    m->called_once = true;

    m->fstats.at("postprocessing").tic();

    m->return_status = flag;
    m->success = flag==qpOASES::SUCCESSFUL_RETURN;

    if (verbose_) casadi_message("qpOASES return status: " + getErrorMessage(m->return_status));

    if (flag!=qpOASES::SUCCESSFUL_RETURN && flag!=qpOASES::RET_MAX_NWSR_REACHED) {
      casadi_error("qpOASES failed: " + getErrorMessage(flag));
    }

    // Get optimal cost
    if (res[CONIC_COST]) *res[CONIC_COST] = m->qp->getObjVal();

    // Get the primal solution
    if (res[CONIC_X]) m->qp->getPrimalSolution(res[CONIC_X]);

    // Get the dual solution
    if (res[CONIC_LAM_X] || res[CONIC_LAM_A]) {
      double* dual=w; w += nx_+na_;
      m->qp->getDualSolution(dual);
      casadi_scal(nx_+na_, -1., dual);
      casadi_copy(dual, nx_, res[CONIC_LAM_X]);
      casadi_copy(dual+nx_, na_, res[CONIC_LAM_A]);
    }

    m->fstats.at("postprocessing").toc();

    // Show statistics
    if (print_time_)  print_fstats(static_cast<ConicMemory*>(mem));
    return 0;
  }

  std::string QpoasesInterface::getErrorMessage(casadi_int flag) {
    switch (flag) {
    case qpOASES::SUCCESSFUL_RETURN:
      return "Successful return.";
    case qpOASES::RET_DIV_BY_ZERO:
      return "Division by zero.";
    case qpOASES::RET_INDEX_OUT_OF_BOUNDS:
      return "Index out of bounds.";
    case qpOASES::RET_INVALID_ARGUMENTS:
      return "At least one of the arguments is invalid.";
    case qpOASES::RET_ERROR_UNDEFINED:
      return "Error number undefined.";
    case qpOASES::RET_WARNING_UNDEFINED:
      return "Warning number undefined.";
    case qpOASES::RET_INFO_UNDEFINED:
      return "Info number undefined.";
    case qpOASES::RET_EWI_UNDEFINED:
      return "Error/warning/info number undefined.";
    case qpOASES::RET_AVAILABLE_WITH_LINUX_ONLY:
      return "This function is available under Linux only.";
    case qpOASES::RET_UNKNOWN_BUG:
      return "The error occured is not yet known.";
    case qpOASES::RET_PRINTLEVEL_CHANGED:
      return "Print level changed.";
    case qpOASES::RET_NOT_YET_IMPLEMENTED:
      return "Requested function is not yet implemented in this version of qpOASES.";
      // Indexlist
    case qpOASES::RET_INDEXLIST_MUST_BE_REORDERD:
      return "Index list has to be reordered.";
    case qpOASES::RET_INDEXLIST_EXCEEDS_MAX_LENGTH:
      return "Index list exceeds its maximal physical length.";
    case qpOASES::RET_INDEXLIST_CORRUPTED:
      return "Index list corrupted.";
    case qpOASES::RET_INDEXLIST_OUTOFBOUNDS:
      return "Physical index is out of bounds.";
    case qpOASES::RET_INDEXLIST_ADD_FAILED:
      return "Adding indices from another index set failed.";
    case qpOASES::RET_INDEXLIST_INTERSECT_FAILED:
      return "Intersection with another index set failed.";
      // SubjectTo / Bounds / Constraints
    case qpOASES::RET_INDEX_ALREADY_OF_DESIRED_STATUS:
      return "Index is already of desired status.";
    case qpOASES::RET_ADDINDEX_FAILED:
      return "Adding index to index set failed.";
    case qpOASES::RET_REMOVEINDEX_FAILED:
      return "Removing index from index set failed.";
    case qpOASES::RET_SWAPINDEX_FAILED:
      return "Cannot swap between different indexsets.";
    case qpOASES::RET_NOTHING_TO_DO:
      return "Nothing to do.";
    case qpOASES::RET_SETUP_BOUND_FAILED:
      return "Setting up bound index failed.";
    case qpOASES::RET_SETUP_CONSTRAINT_FAILED:
      return "Setting up constraint index failed.";
    case qpOASES::RET_MOVING_BOUND_FAILED:
      return "Moving bound between index sets failed.";
    case qpOASES::RET_MOVING_CONSTRAINT_FAILED:
      return "Moving constraint between index sets failed.";
    case qpOASES::RET_SHIFTING_FAILED:
      return "Shifting of bounds/constraints failed.";
    case qpOASES::RET_ROTATING_FAILED:
      return "Rotating of bounds/constraints failed.";
      // QProblem
    case qpOASES::RET_QPOBJECT_NOT_SETUP:
      return "The QP object has not been setup correctly, use another constructor.";
    case qpOASES::RET_QP_ALREADY_INITIALISED:
      return "QProblem has already been initialized.";
    case qpOASES::RET_NO_INIT_WITH_STANDARD_SOLVER:
      return "Initialisation via extern QP solver is not yet implemented.";
    case qpOASES::RET_RESET_FAILED:
      return "Reset failed.";
    case qpOASES::RET_INIT_FAILED:
      return "Initialisation failed.";
    case qpOASES::RET_INIT_FAILED_TQ:
      return "Initialisation failed due to TQ factorisation.";
    case qpOASES::RET_INIT_FAILED_CHOLESKY:
      return "Initialisation failed due to Cholesky decomposition.";
    case qpOASES::RET_INIT_FAILED_HOTSTART:
      return "Initialisation failed! QP could not be solved!";
    case qpOASES::RET_INIT_FAILED_INFEASIBILITY:
      return "Initial QP could not be solved due to infeasibility!";
    case qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS:
      return "Initial QP could not be solved due to unboundedness!";
    case qpOASES::RET_INIT_SUCCESSFUL:
      return "Initialisation done.";
    case qpOASES::RET_OBTAINING_WORKINGSET_FAILED:
      return "Failed to obtain working set for auxiliary QP.";
    case qpOASES::RET_SETUP_WORKINGSET_FAILED:
      return "Failed to setup working set for auxiliary QP.";
    case qpOASES::RET_SETUP_AUXILIARYQP_FAILED:
      return "Failed to setup auxiliary QP for initialized homotopy.";
    case qpOASES::RET_NO_EXTERN_SOLVER:
      return "No extern QP solver available.";
    case qpOASES::RET_QP_UNBOUNDED:
      return "QP is unbounded.";
    case qpOASES::RET_QP_INFEASIBLE:
      return "QP is infeasible.";
    case qpOASES::RET_QP_NOT_SOLVED:
      return "Problems occured while solving QP with standard solver.";
    case qpOASES::RET_QP_SOLVED:
      return "QP successfully solved.";
    case qpOASES::RET_UNABLE_TO_SOLVE_QP:
      return "Problems occured while solving QP.";
    case qpOASES::RET_INITIALISATION_STARTED:
      return "Starting problem initialisation.";
    case qpOASES::RET_HOTSTART_FAILED:
      return "Unable to perform homotopy due to internal error.";
    case qpOASES::RET_HOTSTART_FAILED_TO_INIT:
      return "Unable to initialise problem.";
    case qpOASES::RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED:
      return "Unable to perform homotopy as previous QP is not solved.";
    case qpOASES::RET_ITERATION_STARTED:
      return "Iteration...";
    case qpOASES::RET_SHIFT_DETERMINATION_FAILED:
      return "Determination of shift of the QP data failed.";
    case qpOASES::RET_STEPDIRECTION_DETERMINATION_FAILED:
      return "Determination of step direction failed.";
    case qpOASES::RET_STEPLENGTH_DETERMINATION_FAILED:
      return "Determination of step direction failed.";
    case qpOASES::RET_OPTIMAL_SOLUTION_FOUND:
      return "Optimal solution of neighbouring QP found.";
    case qpOASES::RET_HOMOTOPY_STEP_FAILED:
      return "Unable to perform homotopy step.";
    case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY:
      return "Premature homotopy termination because QP is infeasible.";
    case qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
      return "Premature homotopy termination because QP is unbounded.";
    case qpOASES::RET_WORKINGSET_UPDATE_FAILED:
      return "Unable to update working sets according to initial guesses.";
    case qpOASES::RET_MAX_NWSR_REACHED:
      return "Maximum number of working set recalculations performed.";
    case qpOASES::RET_CONSTRAINTS_NOT_SPECIFIED:
      return "Problem does comprise constraints! "
        "You also have to specify new constraints' bounds.";
    case qpOASES::RET_INVALID_FACTORISATION_FLAG:
      return "Invalid factorisation flag.";
    case qpOASES::RET_UNABLE_TO_SAVE_QPDATA:
      return "Unable to save QP data.";
    case qpOASES::RET_STEPDIRECTION_FAILED_TQ:
      return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY:
      return "Abnormal termination due to Cholesky factorisation.";
    case qpOASES::RET_CYCLING_DETECTED:
      return "Cycling detected.";
    case qpOASES::RET_CYCLING_NOT_RESOLVED:
      return "Cycling cannot be resolved, QP probably infeasible.";
    case qpOASES::RET_CYCLING_RESOLVED:
      return "Cycling probably resolved.";
    case qpOASES::RET_STEPSIZE:
      return "For displaying performed stepsize.";
    case qpOASES::RET_STEPSIZE_NONPOSITIVE:
      return "For displaying non-positive stepsize.";
    case qpOASES::RET_SETUPSUBJECTTOTYPE_FAILED:
      return "Setup of SubjectToTypes failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED:
      return "Addition of constraint to working set failed.";
    case qpOASES::RET_ADDCONSTRAINT_FAILED_INFEASIBILITY:
      return "Addition of constraint to working set failed (due to QP infeasibility).";
    case qpOASES::RET_ADDBOUND_FAILED:
      return "Addition of bound to working set failed.";
    case qpOASES::RET_ADDBOUND_FAILED_INFEASIBILITY:
      return "Addition of bound to working set failed (due to QP infeasibility).";
    case qpOASES::RET_REMOVECONSTRAINT_FAILED:
      return "Removal of constraint from working set failed.";
    case qpOASES::RET_REMOVEBOUND_FAILED:
      return "Removal of bound from working set failed.";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET:
      return "Removing from active set...";
    case qpOASES::RET_ADD_TO_ACTIVESET:
      return "Adding to active set...";
    case qpOASES::RET_REMOVE_FROM_ACTIVESET_FAILED:
      return "Removing from active set failed.";
    case qpOASES::RET_ADD_TO_ACTIVESET_FAILED:
      return "Adding to active set failed.";
    case qpOASES::RET_CONSTRAINT_ALREADY_ACTIVE:
      return "Constraint is already active.";
    case qpOASES::RET_ALL_CONSTRAINTS_ACTIVE:
      return "All constraints are active, no further constraint can be added.";
    case qpOASES::RET_LINEARLY_DEPENDENT:
      return "New bound/constraint is linearly dependent.";
    case qpOASES::RET_LINEARLY_INDEPENDENT:
      return "New bound/constraint is linearly independent.";
    case qpOASES::RET_LI_RESOLVED:
      return "Linear indepence of active contraint matrix successfully resolved.";
    case qpOASES::RET_ENSURELI_FAILED:
      return "Failed to ensure linear indepence of active contraint matrix.";
    case qpOASES::RET_ENSURELI_FAILED_TQ:
      return "Abnormal termination due to TQ factorisation.";
    case qpOASES::RET_ENSURELI_FAILED_NOINDEX:
      return "QP is infeasible.";
    case qpOASES::RET_ENSURELI_FAILED_CYCLING:
      return "QP is infeasible.";
    case qpOASES::RET_BOUND_ALREADY_ACTIVE:
      return "Bound is already active.";
    case qpOASES::RET_ALL_BOUNDS_ACTIVE:
      return "All bounds are active, no further bound can be added.";
    case qpOASES::RET_CONSTRAINT_NOT_ACTIVE:
      return "Constraint is not active.";
    case qpOASES::RET_BOUND_NOT_ACTIVE:
      return "Bound is not active.";
    case qpOASES::RET_HESSIAN_NOT_SPD:
      return "Projected Hessian matrix not positive definite.";
    case qpOASES::RET_HESSIAN_INDEFINITE:
      return "Hessian matrix is indefinite.";
    case qpOASES::RET_MATRIX_SHIFT_FAILED:
      return "Unable to update matrices or to transform vectors.";
    case qpOASES::RET_MATRIX_FACTORISATION_FAILED:
      return "Unable to calculate new matrix factorisations.";
    case qpOASES::RET_PRINT_ITERATION_FAILED:
      return "Unable to print information on current iteration.";
    case qpOASES::RET_NO_GLOBAL_MESSAGE_OUTPUTFILE:
      return "No global message output file initialized.";
    case qpOASES::RET_DISABLECONSTRAINTS_FAILED:
      return "Unable to disbable constraints.";
    case qpOASES::RET_ENABLECONSTRAINTS_FAILED:
      return "Unable to enbable constraints.";
    case qpOASES::RET_ALREADY_ENABLED:
      return "Bound or constraint is already enabled.";
    case qpOASES::RET_ALREADY_DISABLED:
      return "Bound or constraint is already disabled.";
    case qpOASES::RET_NO_HESSIAN_SPECIFIED:
      return "No Hessian matrix has been specified.";
    case qpOASES::RET_USING_REGULARISATION:
      return "Using regularisation as Hessian matrix is not positive definite.";
    case qpOASES::RET_EPS_MUST_BE_POSITVE:
      return "Eps for regularisation must be sufficiently positive.";
    case qpOASES::RET_REGSTEPS_MUST_BE_POSITVE:
      return "Maximum number of regularisation steps must be non-negative.";
    case qpOASES::RET_HESSIAN_ALREADY_REGULARISED:
      return "Hessian has been already regularised.";
    case qpOASES::RET_CANNOT_REGULARISE_IDENTITY:
      return "Identity Hessian matrix cannot be regularised.";
    case qpOASES::RET_NO_REGSTEP_NWSR:
      return "No additional regularisation step could be performed due to limits.";
    case qpOASES::RET_FEWER_REGSTEPS_NWSR:
      return "Fewer additional regularisation steps have been performed due to limits.";
    case qpOASES::RET_CHOLESKY_OF_ZERO_HESSIAN:
      return "Cholesky decomposition of (unregularised) zero Hessian matrix.";
    case qpOASES::RET_CONSTRAINTS_ARE_NOT_SCALED:
      return "When defining __MANY_CONSTRAINTS__, l1 norm of each "
        "constraint must be not greater than one.";
    case qpOASES::RET_ERROR_IN_CONSTRAINTPRODUCT:
      return "Error in user-defined constraint product function.";
      // SQProblem
    case qpOASES::RET_UPDATEMATRICES_FAILED:
      return "Unable to update QP matrices.";
    case qpOASES::RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED:
      return "Unable to update matrices as previous QP is not solved.";
      // Utils
    case qpOASES::RET_UNABLE_TO_OPEN_FILE:
      return "Unable to open file.";
    case qpOASES::RET_UNABLE_TO_WRITE_FILE:
      return "Unable to write into file.";
    case qpOASES::RET_UNABLE_TO_READ_FILE:
      return "Unable to read from file.";
    case qpOASES::RET_FILEDATA_INCONSISTENT:
      return "File contains inconsistent data.";
      // SolutionAnalysis
    case qpOASES::RET_UNABLE_TO_ANALYSE_QPROBLEM:
      return "Unable to analyse (S)QProblem(B) object";
      // Benchmark
    case qpOASES::RET_NWSR_SET_TO_ONE:
      return "Maximum number of working set changes was set to 1.";
    case qpOASES::RET_BENCHMARK_ABORTED:
      return "Benchmark aborted.";
    case qpOASES::RET_UNABLE_TO_READ_BENCHMARK:
      return "Unable to read benchmark data.";
    case qpOASES::RET_INITIAL_QP_SOLVED:
      return "Initial QP solved.";
    case qpOASES::RET_QP_SOLUTION_STARTED:
      return "Solving QP...";
    case qpOASES::RET_BENCHMARK_SUCCESSFUL:
      return "Benchmark terminated successfully.";
    }

    // Default error message
    stringstream ss;
    ss << "Unknown error flag: " << flag << ". Consult qpOASES documentation.";
    return ss.str();
  }

  bool QpoasesInterface::from_BooleanType(qpOASES::BooleanType b) {
    switch (b) {
    case qpOASES::BT_TRUE:              return true;
    case qpOASES::BT_FALSE:             return false;
    }
    casadi_error("not_implemented");
  }

  qpOASES::BooleanType QpoasesInterface::to_BooleanType(bool b) {
    return b ? qpOASES::BT_TRUE : qpOASES::BT_FALSE;
  }

  std::string QpoasesInterface::from_SubjectToStatus(qpOASES::SubjectToStatus b) {
    switch (b) {
    case qpOASES::ST_INACTIVE:          return "inactive";
    case qpOASES::ST_LOWER:             return "lower";
    case qpOASES::ST_UPPER:             return "upper";
    case qpOASES::ST_INFEASIBLE_LOWER:  return "infeasible_lower";
    case qpOASES::ST_INFEASIBLE_UPPER:  return "infeasible_upper";
    case qpOASES::ST_UNDEFINED:         return "undefined";
    }
    casadi_error("not_implemented");
  }

  qpOASES::SubjectToStatus QpoasesInterface::to_SubjectToStatus(std::string b) {
    if (b == "inactive") {
      return qpOASES::ST_INACTIVE;
    } else if (b == "lower") {
      return qpOASES::ST_LOWER;
    } else if (b == "infeasible_lower") {
      return qpOASES::ST_INFEASIBLE_LOWER;
    } else if (b == "infeasible_upper") {
      return qpOASES::ST_INFEASIBLE_UPPER;
    } else if (b == "undefined") {
      return qpOASES::ST_UNDEFINED;
    } else {
      casadi_error("No such qpOASES::SubjectToStatus: " + b);
    }
  }

  std::string QpoasesInterface::from_PrintLevel(qpOASES::PrintLevel b) {
    switch (b) {
    case qpOASES::PL_TABULAR:           return "tabular";
    case qpOASES::PL_NONE:              return "none";
    case qpOASES::PL_LOW:               return "low";
    case qpOASES::PL_MEDIUM:            return "medium";
    case qpOASES::PL_HIGH:              return "high";
    case qpOASES::PL_DEBUG_ITER:        return "debug_iter";
    }
    casadi_error("not_implemented");
  }

  qpOASES::PrintLevel QpoasesInterface::to_PrintLevel(std::string b) {
    if (b == "tabular") {
      return qpOASES::PL_TABULAR;
    } else if (b == "none") {
      return qpOASES::PL_NONE;
    } else if (b == "low") {
      return qpOASES::PL_LOW;
    } else if (b == "medium") {
      return qpOASES::PL_MEDIUM;
    } else if (b == "high") {
      return qpOASES::PL_HIGH;
    } else if (b == "debug_iter") {
      return qpOASES::PL_DEBUG_ITER;
    } else {
      casadi_error("No such qpOASES::PrintLevel: " + b);
    }
  }

  QpoasesMemory::QpoasesMemory() {
    this->qp = nullptr;
    this->h = nullptr;
    this->a = nullptr;
  }

  QpoasesMemory::~QpoasesMemory() {
    if (this->qp) delete this->qp;
    if (this->h) delete this->h;
    if (this->a) delete this->a;
  }

  int QpoasesInterface::
  qpoases_init(void* mem, int dim, int nnz, const int* row, const int* col) {
    casadi_assert_dev(mem!=nullptr);
    QpoasesMemory* m = static_cast<QpoasesMemory*>(mem);

    // Get sparsity pattern in sparse triplet format
    m->row.clear();
    m->col.clear();
    m->nz_map.clear();
    for (casadi_int k=0; k<nnz; ++k) {
      // Add upper(?) triangular part (and diagonal)
      m->row.push_back(row[k]-1);
      m->col.push_back(col[k]-1);
      m->nz_map.push_back(k);
      // Add lower(?) triangular part
      if (row[k]!=col[k]) {
        m->row.push_back(col[k]-1);
        m->col.push_back(row[k]-1);
        m->nz_map.push_back(k);
      }
    }

    // Create sparsity pattern: TODO(@jaeandersson) No memory allocation
    Sparsity sp = Sparsity::triplet(dim, dim, m->row, m->col, m->lin_map, false);
    for (casadi_int& e : m->lin_map) e = m->nz_map[e];

    // Allocate memory for nonzeros
    m->nz.resize(sp.nnz());

    // Create linear solver
    m->linsol = Linsol("linsol", m->linsol_plugin, sp);

    return 0;
  }

  int QpoasesInterface::qpoases_sfact(void* mem, const double* vals) {
    casadi_assert_dev(mem!=nullptr);
    QpoasesMemory* m = static_cast<QpoasesMemory*>(mem);

    // Get nonzero elements (entire elements)
    for (int i=0; i<m->nz.size(); ++i) m->nz[i] = vals[m->lin_map[i]];

    // Pass to linear solver
    m->linsol.sfact(get_ptr(m->nz));

    return 0;
  }

  int QpoasesInterface::
  qpoases_nfact(void* mem, const double* vals, int* neig, int* rank) {
    casadi_assert_dev(mem!=nullptr);
    QpoasesMemory* m = static_cast<QpoasesMemory*>(mem);

    // Get nonzero elements (entire elements)
    for (casadi_int i=0; i<m->nz.size(); ++i) m->nz[i] = vals[m->lin_map[i]];

    // Pass to linear solver
    m->linsol.nfact(get_ptr(m->nz));

    // Number of negative eigenvalues
    if (neig) *neig = m->linsol.neig(get_ptr(m->nz));

    // Rank of the matrix
    if (rank) *rank = m->linsol.rank(get_ptr(m->nz));

    return 0;
  }

  int QpoasesInterface::qpoases_solve(void* mem, int nrhs, double* rhs) {
    casadi_assert_dev(mem!=nullptr);
    QpoasesMemory* m = static_cast<QpoasesMemory*>(mem);

    // Pass to linear solver
    m->linsol.solve(get_ptr(m->nz), rhs, nrhs);

    return 0;
  }

  Dict QpoasesInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpoasesMemory*>(mem);
    stats["return_status"] = getErrorMessage(m->return_status);
    stats["success"] = m->success;
    return stats;
  }

} // namespace casadi
