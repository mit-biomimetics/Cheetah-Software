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


#include "worhp_interface.hpp"

#include "casadi/core/casadi_misc.hpp"
#include <ctime>
#include <cstring>

using namespace std;

namespace casadi {

  extern "C"
  int CASADI_NLPSOL_WORHP_EXPORT
  casadi_register_nlpsol_worhp(Nlpsol::Plugin* plugin) {
    plugin->creator = WorhpInterface::creator;
    plugin->name = "worhp";
    plugin->doc = WorhpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &WorhpInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_WORHP_EXPORT casadi_load_nlpsol_worhp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_worhp);
  }

  WorhpInterface::WorhpInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  WorhpInterface::~WorhpInterface() {
    clear_mem();
  }

  Options WorhpInterface::options_
  = {{&Nlpsol::options_},
     {{"worhp",
       {OT_DICT,
        "Options to be passed to WORHP"}}
     }
  };

  void WorhpInterface::init(const Dict& opts) {

    // Call the init method of the base class
    Nlpsol::init(opts);

    if (CheckWorhpVersion(WORHP_MAJOR, WORHP_MINOR, WORHP_PATCH)) {
      casadi_warning("Worhp incompatibility. Interface was compiled for Worhp " +
        str(WORHP_MAJOR) + "." + str(WORHP_MINOR) + "." + std::string(WORHP_PATCH));
    }

    // Default options
    Dict worhp_opts;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="worhp") {
        worhp_opts = op.second;
      }
    }

    // Sort Worhp options
    casadi_int nopts = WorhpGetParamCount();
    for (auto&& op : worhp_opts) {
      if (op.first.compare("qp")==0) {
        qp_opts_ = op.second;
        continue;
      }

      // Get corresponding index using a linear search
      casadi_int ind;
      for (ind=1; ind<=nopts; ++ind) {
        // Get name in WORHP
        const char* name = WorhpGetParamName(ind);
        // Break if matching name
        if (op.first.compare(name)==0) break;
      }
      if (ind>nopts) casadi_error("No such Worhp option: " + op.first);

      // Add to the corresponding list
      switch (WorhpGetParamType(ind)) {
      case WORHP_BOOL_T:
        bool_opts_[op.first] = op.second;
        break;
      case WORHP_DOUBLE_T:
        double_opts_[op.first] = op.second;
        break;
      case WORHP_INT_T:
        int_opts_[op.first] = op.second;
        break;
      default:
        casadi_error("Cannot handle WORHP option \"" + op.first + "\": Unknown type " +
          str(WorhpGetParamType(ind)) + ".");
        break;
      }
    }

    // Setup NLP functions
    f_fcn_ = create_function("nlp_f", {"x", "p"}, {"f"});
    g_fcn_ = create_function("nlp_g", {"x", "p"}, {"g"});
    grad_f_fcn_ = create_function("nlp_grad_f", {"x", "p"}, {"f", "grad:f:x"});
    jac_g_fcn_ = create_function("nlp_jac_g", {"x", "p"}, {"g", "jac:g:x"});
    jacg_sp_ = jac_g_fcn_.sparsity_out(1);
    hess_l_fcn_ = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                  {"transpose:hess:gamma:x:x"},
                                  {{"gamma", {"f", "g"}}});
    hesslag_sp_ = hess_l_fcn_.sparsity_out(0);

    // Temporary vectors
    alloc_w(nx_); // for fetching diagonal entries form Hessian
  }

  void worhp_disp(int mode, const char message[]) {
    if (mode & WORHP_PRINT_MESSAGE) {
      uout() << message << std::endl;
    }
    if (mode & WORHP_PRINT_WARNING) {
      uerr() << message << std::endl;
    }
    if (mode & WORHP_PRINT_ERROR) {
      uerr() << message << std::endl;
    }
  }

  int WorhpInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<WorhpMemory*>(mem);

    SetWorhpPrint(&worhp_disp);

    WorhpPreInit(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
    m->worhp_o.initialised = false;
    m->worhp_w.initialised = false;
    m->worhp_p.initialised = false;
    m->worhp_c.initialised = false;

    // Initialize parameters to default values
    int status;
    InitParams(&status, &m->worhp_p);
    casadi_assert(status==0, "Problem in Worhp InitParams. Status: " + str(status));


    // Pass boolean parameters
    for (auto&& op : bool_opts_) {
      casadi_assert(
        WorhpSetBoolParam(&m->worhp_p, op.first.c_str(), op.second),
        "Problem setting boolean Worhp parameter " + op.first);
    }

    // Pass double parameters
    for (auto&& op : double_opts_) {
      casadi_assert(
        WorhpSetDoubleParam(&m->worhp_p, op.first.c_str(), op.second),
        "Problem setting double Worhp parameter " + op.first);
    }

    // Pass integer parameters
    for (auto&& op : int_opts_) {
      casadi_assert(
        WorhpSetIntParam(&m->worhp_p, op.first.c_str(), op.second),
        "Problem setting integer Worhp parameter " + op.first);
    }

    // Pass qp parameters
    for (auto&& op : qp_opts_) {
      if (op.first=="ipBarrier") {
        m->worhp_p.qp.ipBarrier = op.second;
      } else if (op.first=="ipComTol") {
        m->worhp_p.qp.ipComTol = op.second;
      } else if (op.first=="ipFracBound") {
        m->worhp_p.qp.ipFracBound = op.second;
      } else if (op.first=="ipMinAlpha") {
        m->worhp_p.qp.ipMinAlpha = op.second;
      } else if (op.first=="ipRelaxDiv") {
        m->worhp_p.qp.ipRelaxDiv = op.second;
      } else if (op.first=="ipRelaxMax") {
        m->worhp_p.qp.ipRelaxMax = op.second;
      } else if (op.first=="ipRelaxMin") {
        m->worhp_p.qp.ipRelaxMin = op.second;
      } else if (op.first=="ipRelaxMult") {
        m->worhp_p.qp.ipRelaxMult = op.second;
      } else if (op.first=="ipResTol") {
        m->worhp_p.qp.ipResTol = op.second;
      } else if (op.first=="lsTol") {
        m->worhp_p.qp.lsTol = op.second;
      } else if (op.first=="nsnBeta") {
        m->worhp_p.qp.nsnBeta = op.second;
      } else if (op.first=="nsnKKT") {
        m->worhp_p.qp.nsnKKT = op.second;
      } else if (op.first=="nsnMinAlpha") {
        m->worhp_p.qp.nsnMinAlpha = op.second;
      } else if (op.first=="nsnSigma") {
        m->worhp_p.qp.nsnSigma = op.second;
      } else if (op.first=="ipLsMethod") {
        m->worhp_p.qp.ipLsMethod = op.second;
      } else if (op.first=="lsItMaxIter") {
        m->worhp_p.qp.lsItMaxIter = op.second;
      } else if (op.first=="lsItMethod") {
        m->worhp_p.qp.lsItMethod = op.second;
      } else if (op.first=="lsItPrecondMethod") {
        m->worhp_p.qp.lsItPrecondMethod = op.second;
      } else if (op.first=="lsRefineMaxIter") {
        m->worhp_p.qp.lsRefineMaxIter = op.second;
      } else if (op.first=="maxIter") {
        m->worhp_p.qp.maxIter = op.second;
      } else if (op.first=="method") {
        m->worhp_p.qp.method = op.second;
      } else if (op.first=="nsnLsMethod") {
        m->worhp_p.qp.nsnLsMethod = op.second;
      } else if (op.first=="printLevel") {
        m->worhp_p.qp.printLevel = op.second;
      } else if (op.first=="ipTryRelax") {
        m->worhp_p.qp.ipTryRelax = op.second;
      } else if (op.first=="lsScale") {
        m->worhp_p.qp.lsScale = op.second;
      } else if (op.first=="lsTrySimple") {
        m->worhp_p.qp.lsTrySimple = op.second;
      } else if (op.first=="nsnGradStep") {
        m->worhp_p.qp.nsnGradStep = op.second;
      } else if (op.first=="scaleIntern") {
        m->worhp_p.qp.scaleIntern = op.second;
      } else if (op.first=="strict") {
        m->worhp_p.qp.strict = op.second;
      } else {
        casadi_error("No such Worhp option: qp." + op.first);
      }
    }

    // Mark the parameters as set
    m->worhp_p.initialised = true;
    m->init_ = false;

    return 0;
  }

  void WorhpInterface::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<WorhpMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Free existing Worhp memory (except parameters)
    m->worhp_p.initialised = false; // Avoid freeing the memory for parameters
    if (m->worhp_o.initialised || m->worhp_w.initialised || m->worhp_c.initialised) {
      WorhpFree(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
    }
    m->worhp_p.initialised = true;

    // Number of (free) variables
    m->worhp_o.n = nx_;

    // Number of constraints
    m->worhp_o.m = ng_;

    /// Control data structure needs to be reset every time
    m->worhp_c.initialised = false;
    m->worhp_w.initialised = false;
    m->worhp_o.initialised = false;

    // Worhp uses the CS format internally, hence it is the preferred sparse matrix format.
    m->worhp_w.DF.nnz = nx_;
    if (m->worhp_o.m>0) {
      m->worhp_w.DG.nnz = jacg_sp_.nnz();  // Jacobian of G
    } else {
      m->worhp_w.DG.nnz = 0;
    }

    if (true /*m->worhp_w.HM.NeedStructure*/) { // not initialized
      m->worhp_w.HM.nnz = nx_ + hesslag_sp_.nnz_lower(true);
    } else {
      m->worhp_w.HM.nnz = 0;
    }

    /* Data structure initialisation. */
    WorhpInit(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
    m->init_ = true;
    if (m->worhp_c.status != FirstCall) {
      string msg = return_codes(m->worhp_c.status);
      casadi_error("Main: Initialisation failed. Status: " + msg);
    }

    if (m->worhp_w.DF.NeedStructure) {
      for (casadi_int i=0; i<nx_; ++i) {
        m->worhp_w.DF.row[i] = i + 1; // Index-1 based
      }
    }

    if (m->worhp_o.m>0 && m->worhp_w.DG.NeedStructure) {
      casadi_int nz=0;
      const casadi_int* colind = jacg_sp_.colind();
      const casadi_int* row = jacg_sp_.row();
      for (casadi_int c=0; c<nx_; ++c) {
        for (casadi_int el=colind[c]; el<colind[c+1]; ++el) {
          casadi_int r = row[el];
          m->worhp_w.DG.col[nz] = c + 1; // Index-1 based
          m->worhp_w.DG.row[nz] = r + 1;
          nz++;
        }
      }
    }

    if (m->worhp_w.HM.NeedStructure) {
      // Get the sparsity pattern of the Hessian
      const casadi_int* colind = hesslag_sp_.colind();
      const casadi_int* row = hesslag_sp_.row();

      casadi_int nz=0;

      // Strictly lower triangular part of the Hessian (note CCS -> CRS format change)
      for (casadi_int c=0; c<nx_; ++c) {
        for (casadi_int el=colind[c]; el<colind[c+1]; ++el) {
          if (row[el]>c) {
            m->worhp_w.HM.row[nz] = row[el] + 1;
            m->worhp_w.HM.col[nz] = c + 1;
            nz++;
          }
        }
      }

      // Diagonal always included
      for (casadi_int r=0; r<nx_; ++r) {
        m->worhp_w.HM.row[nz] = r + 1;
        m->worhp_w.HM.col[nz] = r + 1;
        nz++;
      }
    }
  }

  int WorhpInterface::solve(void* mem) const {
    auto m = static_cast<WorhpMemory*>(mem);

    if (m->lbg && m->ubg) {
      for (casadi_int i=0; i<ng_; ++i) {
        casadi_assert(!(m->lbg[i]==-inf && m->ubg[i] == inf),
                        "WorhpInterface::evaluate: Worhp cannot handle the case when both "
                        "LBG and UBG are infinite."
                        "You have that case at non-zero " + str(i)+ "."
                        "Reformulate your problem eliminating the corresponding constraint.");
      }
    }

    // Pass inputs to WORHP data structures
    casadi_copy(m->x, nx_, m->worhp_o.X);
    casadi_copy(m->lbx, nx_, m->worhp_o.XL);
    casadi_copy(m->ubx, nx_, m->worhp_o.XU);
    casadi_copy(m->lam_x, nx_, m->worhp_o.Lambda);
    if (m->worhp_o.m>0) {
      casadi_copy(m->lam_g, ng_, m->worhp_o.Mu);
      casadi_copy(m->lbg, ng_, m->worhp_o.GL);
      casadi_copy(m->ubg, ng_, m->worhp_o.GU);
    }

    // Replace infinite bounds with m->worhp_p.Infty
    double inf = numeric_limits<double>::infinity();
    for (casadi_int i=0; i<nx_; ++i)
      if (m->worhp_o.XL[i]==-inf) m->worhp_o.XL[i] = -m->worhp_p.Infty;
    for (casadi_int i=0; i<nx_; ++i)
      if (m->worhp_o.XU[i]== inf) m->worhp_o.XU[i] =  m->worhp_p.Infty;
    for (casadi_int i=0; i<ng_; ++i)
      if (m->worhp_o.GL[i]==-inf) m->worhp_o.GL[i] = -m->worhp_p.Infty;
    for (casadi_int i=0; i<ng_; ++i)
      if (m->worhp_o.GU[i]== inf) m->worhp_o.GU[i] =  m->worhp_p.Infty;

    if (verbose_) casadi_message("WorhpInterface::starting iteration");

    bool firstIteration = true;

    // Reverse Communication loop
    while (m->worhp_c.status < TerminateSuccess &&  m->worhp_c.status > TerminateError) {
      if (GetUserAction(&m->worhp_c, callWorhp)) {
        Worhp(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
      }


      if (GetUserAction(&m->worhp_c, iterOutput)) {

        if (!firstIteration) {
          firstIteration = true;

          if (!fcallback_.is_null()) {
            m->iter = m->worhp_w.MajorIter;
            m->iter_sqp = m->worhp_w.MinorIter;
            m->inf_pr = m->worhp_w.NormMax_CV;
            m->inf_du = m->worhp_p.ScaledKKT;
            m->alpha_pr = m->worhp_w.ArmijoAlpha;

            // Inputs
            fill_n(m->arg, fcallback_.n_in(), nullptr);
            m->arg[NLPSOL_X] = m->worhp_o.X;
            m->arg[NLPSOL_F] = &m->worhp_o.F;
            m->arg[NLPSOL_G] = m->worhp_o.G;
            m->arg[NLPSOL_LAM_P] = nullptr;
            m->arg[NLPSOL_LAM_X] = m->worhp_o.Lambda;
            m->arg[NLPSOL_LAM_G] = m->worhp_o.Mu;

            // Outputs
            fill_n(m->res, fcallback_.n_out(), nullptr);
            double ret_double;
            m->res[0] = &ret_double;

            m->fstats.at("callback_fun").tic();
            // Evaluate the callback function
            fcallback_(m->arg, m->res, m->iw, m->w, 0);
            m->fstats.at("callback_fun").toc();
            casadi_int ret = static_cast<casadi_int>(ret_double);

            if (ret) m->worhp_c.status = TerminateError;
          }
        }


        IterationOutput(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
        DoneUserAction(&m->worhp_c, iterOutput);
      }

      if (GetUserAction(&m->worhp_c, evalF)) {
        m->arg[0] = m->worhp_o.X;
        m->arg[1] = m->p;
        m->res[0] = &m->worhp_o.F;
        calc_function(m, "nlp_f");
        m->f = m->worhp_o.F; // Store cost, before scaling
        m->worhp_o.F *= m->worhp_w.ScaleObj;
        DoneUserAction(&m->worhp_c, evalF);
      }

      if (GetUserAction(&m->worhp_c, evalG)) {
        m->arg[0] = m->worhp_o.X;
        m->arg[1] = m->p;
        m->res[0] = m->worhp_o.G;
        calc_function(m, "nlp_g");
        DoneUserAction(&m->worhp_c, evalG);
      }

      if (GetUserAction(&m->worhp_c, evalDF)) {
        m->arg[0] = m->worhp_o.X;
        m->arg[1] = m->p;
        m->res[0] = nullptr;
        m->res[1] = m->worhp_w.DF.val;
        calc_function(m, "nlp_grad_f");
        casadi_scal(nx_, m->worhp_w.ScaleObj, m->worhp_w.DF.val);
        DoneUserAction(&m->worhp_c, evalDF);
      }

      if (GetUserAction(&m->worhp_c, evalDG)) {
        m->arg[0] = m->worhp_o.X;
        m->arg[1] = m->p;
        m->res[0] = nullptr;
        m->res[1] = m->worhp_w.DG.val;
        calc_function(m, "nlp_jac_g");
        DoneUserAction(&m->worhp_c, evalDG);
      }

      if (GetUserAction(&m->worhp_c, evalHM)) {
        m->arg[0] = m->worhp_o.X;
        m->arg[1] = m->p;
        m->arg[2] = &m->worhp_w.ScaleObj;
        m->arg[3] = m->worhp_o.Mu;
        m->res[0] = m->worhp_w.HM.val;
        calc_function(m, "nlp_hess_l");
        // Diagonal values
        double *dval = m->w;
        casadi_fill(dval, nx_, 0.);

        // Remove diagonal
        const casadi_int* colind = hesslag_sp_.colind();
        const casadi_int* row = hesslag_sp_.row();
        casadi_int ind=0;
        for (casadi_int c=0; c<nx_; ++c) {
          for (casadi_int el=colind[c]; el<colind[c+1]; ++el) {
            if (row[el]==c) {
              dval[c] = m->worhp_w.HM.val[el];
            } else {
              m->worhp_w.HM.val[ind++] = m->worhp_w.HM.val[el];
            }
          }
        }

        // Add diagonal entries at the end
        casadi_copy(dval, nx_, m->worhp_w.HM.val+ind);
        DoneUserAction(&m->worhp_c, evalHM);
      }

      if (GetUserAction(&m->worhp_c, fidif)) {
        WorhpFidif(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);
      }
    }

    // Copy outputs
    casadi_copy(m->worhp_o.X, nx_, m->x);
    casadi_copy(m->worhp_o.G, ng_, m->g);
    casadi_copy(m->worhp_o.Lambda, nx_, m->lam_x);
    casadi_copy(m->worhp_o.Mu, ng_, m->lam_g);

    StatusMsg(&m->worhp_o, &m->worhp_w, &m->worhp_p, &m->worhp_c);

    m->return_code = m->worhp_c.status;
    m->return_status = return_codes(m->worhp_c.status);
    m->success = m->return_code > TerminateSuccess;
    return 0;
  }

  const char* WorhpInterface::return_codes(casadi_int flag) {
    switch (flag) {
    case OptimalSolution: return "OptimalSolution";
    case OptimalSolutionConstantF: return "OptimalSolutionConstantF";
    case SearchDirectionZero: return "SearchDirectionZero";
    case SearchDirectionSmall: return "SearchDirectionSmall";
    case FritzJohn: return "FritzJohn";
    case NotDiffable: return "NotDiffable";
    case Unbounded: return "Unbounded";
    case FeasibleSolution: return "FeasibleSolution";
    case LowPassFilterOptimal: return "LowPassFilterOptimal";
    case LowPassFilterAcceptable: return "LowPassFilterAcceptable";
    case AcceptableSolution: return "AcceptableSolution";
    case AcceptablePrevious: return "AcceptablePrevious";
    case AcceptableSolutionConstantF: return "AcceptableSolutionConstantF";
    case AcceptablePreviousConstantF: return "AcceptablePreviousConstantF";
    case AcceptableSolutionSKKT: return "AcceptableSolutionSKKT";
    case AcceptableSolutionScaled: return "AcceptableSolutionScaled";
    case AcceptablePreviousScaled: return "AcceptablePreviousScaled";
    case TerminateError: return "TerminateError";
    case MaxCalls: return "MaxCalls";
    case MaxIter: return "MaxIter";
    case Timeout: return "Timeout";
    case TooBig: return "TooBig";
    case evalsNaN: return "evalsNaN";
    case DivergingPrimal: return "DivergingPrimal";
    case DivergingDual: return "DivergingDual";
    case MinimumStepsize: return "MinimumStepsize";
    case RegularizationFailed: return "RegularizationFailed";
    case InitError: return "InitError";
    case DataError: return "DataError";
    case RestartError: return "RestartError";
    case QPerror: return "QPerror";
    case LinearSolverFailed: return "LinearSolverFailed";
    case TerminatedByCheckFD: return "TerminatedByCheckFD";
    case LicenseError: return "LicenseError";
    case Debug: return "Debug";
    }
    return "Unknown WORHP return code";
  }

  WorhpMemory::WorhpMemory() {
    this->worhp_o.initialised = false;
    this->worhp_w.initialised = false;
    this->worhp_p.initialised = false;
    this->worhp_c.initialised = false;
  }

  WorhpMemory::~WorhpMemory() {
    if (this->init_) {
      if (this->worhp_p.initialised || this->worhp_o.initialised ||
          this->worhp_w.initialised || this->worhp_c.initialised) {
        WorhpFree(&this->worhp_o, &this->worhp_w, &this->worhp_p, &this->worhp_c);
      }
    }
  }

  Dict WorhpInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<WorhpMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

} // namespace casadi
