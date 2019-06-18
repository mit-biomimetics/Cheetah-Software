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


#include "kinsol_interface.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_ROOTFINDER_KINSOL_EXPORT
  casadi_register_rootfinder_kinsol(Rootfinder::Plugin* plugin) {
    plugin->creator = KinsolInterface::creator;
    plugin->name = "kinsol";
    plugin->doc = KinsolInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &KinsolInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_KINSOL_EXPORT casadi_load_rootfinder_kinsol() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_kinsol);
  }

  /** \brief Kinsol solver class
   *
   * @copydoc Rootfinder_doc
   * You can provide an initial guess by setting output(0).\n
   * A good initial guess may be needed to avoid errors like
   * "The linear solver's setup function failed in an unrecoverable manner."
   *
   The constraints option expects an integer entry for each variable u:\n

   0 then no constraint is imposed on \p ui. \n
   1 then \p ui will be constrained to be \p ui >= 0.0. \n
   -1 then \p ui will be constrained to be \p ui <= 0.0. \n
   2 then \p ui will be constrained to be \p ui > 0.0. \n
   -2 then \p ui will be constrained to be \p ui < 0.0. \n

   *
   * \see Rootfinder for more information
   *
   */
  KinsolInterface::KinsolInterface(const std::string& name, const Function& f)
    : Rootfinder(name, f) {

    u_scale_ = nullptr;
    f_scale_ = nullptr;

    // Default options
    exact_jac_ = true;
    disable_internal_warnings_ = false;
    max_iter_ = 0;
    maxl_ = 0;
    upper_bandwidth_ = -1;
    lower_bandwidth_ = -1;
    use_preconditioner_ = false;
    abstol_ = 1e-6;
  }

  KinsolInterface::~KinsolInterface() {
    if (u_scale_) N_VDestroy_Serial(u_scale_);
    if (f_scale_) N_VDestroy_Serial(f_scale_);
    clear_mem();
  }

  Options KinsolInterface::options_
  = {{&Rootfinder::options_},
     {{"max_iter",
       {OT_INT,
        "Maximum number of Newton iterations. Putting 0 sets the default value of KinSol."}},
      {"abstol",
       {OT_DOUBLE,
        "Stopping criterion tolerance"}},
      {"linear_solver_type",
       {OT_STRING,
        "dense|banded|iterative|user_defined"}},
      {"upper_bandwidth",
       {OT_INT,
        "Upper bandwidth for banded linear solvers"}},
      {"lower_bandwidth",
       {OT_INT,
        "Lower bandwidth for banded linear solvers"}},
      {"max_krylov",
       {OT_INT,
        "Maximum Krylov space dimension"}},
      {"exact_jacobian",
       {OT_BOOL,
        "Use exact Jacobian information"}},
      {"iterative_solver",
       {OT_STRING,
        "gmres|bcgstab|tfqmr"}},
      {"f_scale",
       {OT_DOUBLEVECTOR,
        "Equation scaling factors"}},
      {"u_scale",
       {OT_DOUBLEVECTOR,
        "Variable scaling factors"}},
      {"pretype",
       {OT_STRING,
        "Type of preconditioner"}},
      {"use_preconditioner",
       {OT_BOOL,
        "Precondition an iterative solver"}},
      {"strategy",
       {OT_STRING,
        "Globalization strategy"}},
      {"disable_internal_warnings",
       {OT_BOOL,
        "Disable KINSOL internal warning messages"}}
     }
  };

  void KinsolInterface::init(const Dict& opts) {
    // Initialize the base classes
    Rootfinder::init(opts);

    // Default (temporary) options
    string strategy = "none";
    vector<double> u_scale;
    vector<double> f_scale;
    string linear_solver_type = "dense";
    string iterative_solver = "gmres";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="strategy") {
        strategy = op.second.to_string();
      } else if (op.first=="exact_jacobian") {
        exact_jac_ = op.second;
      } else if (op.first=="u_scale") {
        u_scale = op.second;
      } else if (op.first=="f_scale") {
        f_scale = op.second;
      } else if (op.first=="disable_internal_warnings") {
        disable_internal_warnings_ = op.second;
      } else if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="linear_solver_type") {
        linear_solver_type = op.second.to_string();
      } else if (op.first=="max_krylov") {
        maxl_ = op.second;
      } else if (op.first=="upper_bandwidth") {
        upper_bandwidth_ = op.second;
      } else if (op.first=="lower_bandwidth") {
        lower_bandwidth_ = op.second;
      } else if (op.first=="iterative_solver") {
        iterative_solver = op.second.to_string();
      } else if (op.first=="use_preconditioner") {
        use_preconditioner_ = op.second;
      } else if (op.first=="abstol") {
        abstol_ = op.second;
      }
    }

    // Get globalization strategy
    if (strategy=="linesearch") {
      strategy_ = KIN_LINESEARCH;
    } else {
      casadi_assert_dev(strategy=="none");
      strategy_ = KIN_NONE;
    }

    // Allocate N_Vectors
    if (u_scale_) N_VDestroy_Serial(u_scale_);
    if (f_scale_) N_VDestroy_Serial(f_scale_);
    u_scale_ = N_VNew_Serial(n_);
    f_scale_ = N_VNew_Serial(n_);

    // Set scaling factors on variables
    if (!u_scale.empty()) {
      casadi_assert_dev(u_scale.size()==NV_LENGTH_S(u_scale_));
      copy(u_scale.begin(), u_scale.end(), NV_DATA_S(u_scale_));
    } else {
      N_VConst(1.0, u_scale_);
    }

    // Set scaling factors on equations
    if (!f_scale.empty()) {
      casadi_assert_dev(f_scale.size()==NV_LENGTH_S(f_scale_));
      copy(f_scale.begin(), f_scale.end(), NV_DATA_S(f_scale_));
    } else {
      N_VConst(1.0, f_scale_);
    }

    // Type of linear solver
    if (linear_solver_type=="dense") {
      linear_solver_type_ = DENSE;
      if (exact_jac_) {
        // For storing Jacobian nonzeros
        alloc_w(sp_jac_.nnz(), true);
      }
    } else if (linear_solver_type=="banded") {
      linear_solver_type_ = BANDED;
      casadi_assert_dev(upper_bandwidth_>=0);
      casadi_assert_dev(lower_bandwidth_>=0);
      if (exact_jac_) {
        // For storing Jacobian nonzeros
        alloc_w(sp_jac_.nnz(), true);
      }
    } else if (linear_solver_type=="iterative") {
      linear_solver_type_ = ITERATIVE;
      if (iterative_solver=="gmres") {
        iterative_solver_ = GMRES;
      } else if (iterative_solver=="bcgstab") {
        iterative_solver_ = BCGSTAB;
      } else if (iterative_solver=="tfqmr") {
        iterative_solver_ = TFQMR;
      } else {
        casadi_error("KINSOL: Unknown sparse solver");
      }
      if (exact_jac_) {
        get_jtimes();
      }
    } else if (linear_solver_type=="user_defined") {
      linear_solver_type_ = USER_DEFINED;

      // Form the Jacobian-times-vector function
      get_jtimes();

      // Allocate space for Jacobian
      alloc_w(sp_jac_.nnz(), true);
    } else {
      casadi_error("Unknown linear solver");
    }
  }

  void KinsolInterface::set_work(void* mem, const double**& arg, double**& res,
                        casadi_int*& iw, double*& w) const {
      Rootfinder::set_work(mem, arg, res, iw, w);
      auto m = static_cast<KinsolMemory*>(mem);
      m->jac = w; w += sp_jac_.nnz();
   }

  void KinsolInterface::get_jtimes() {
    vector<string> jtimes_in = oracle_.name_in();
    jtimes_in.push_back("fwd:" + oracle_.name_in(iin_));
    vector<string> jtimes_out = {"fwd:" + oracle_.name_out(iout_)};
    jtimes_ = oracle_.factory("jtimes", jtimes_in, jtimes_out);
    alloc(jtimes_);
  }

  int KinsolInterface::solve(void* mem) const {
    auto m = static_cast<KinsolMemory*>(mem);

    // Get the initial guess
    casadi_copy(m->iarg[iin_], nnz_in(iin_), NV_DATA_S(m->u));

    // Solve the nonlinear system of equations
    int flag = KINSol(m->mem, m->u, strategy_, u_scale_, f_scale_);
    m->success = flag>= KIN_SUCCESS;
    if (flag<KIN_SUCCESS) kinsol_error("KINSol", flag, error_on_fail_);

    // Warn if not successful return
    if (verbose_) {
      if (flag!=KIN_SUCCESS) kinsol_error("KINSol", flag, false);
    }

    // Get the solution
    casadi_copy(NV_DATA_S(m->u), nnz_out(iout_), m->ires[iout_]);

    // Evaluate auxiliary outputs
    if (n_out_>0) {
      // Evaluate f_
      copy_n(m->iarg, n_in_, m->arg);
      m->arg[iin_] = NV_DATA_S(m->u);
      copy_n(m->ires, n_out_, m->res);
      m->res[iout_] = nullptr;
      oracle_(m->arg, m->res, m->iw, m->w, 0);
    }
    return 0;
  }

  void KinsolInterface::func(KinsolMemory& m, N_Vector u, N_Vector fval) const {
    // Evaluate f_
    copy_n(m.iarg, n_in_, m.arg);
    m.arg[iin_] = NV_DATA_S(u);
    fill_n(m.res, n_out_, nullptr);
    m.res[iout_] = NV_DATA_S(fval);
    oracle_(m.arg, m.res, m.iw, m.w, 0);

    // Make sure that all entries of the linear system are valid
    double *fdata = NV_DATA_S(fval);
    for (int k=0; k<n_; ++k) {
      try {
        casadi_assert(!isnan(fdata[k]),
          "Nonzero " + str(k) + " is not-a-number");
        casadi_assert(!isinf(fdata[k]),
          "Nonzero " + str(k) + " is infinite");
      } catch(exception& ex) {
        stringstream ss;
        ss << ex.what() << endl;
        if (verbose_) {
          uout() << "u = ";
          N_VPrint_Serial(u);
        }

        throw CasadiException(ss.str());
      }
    }
  }

  int KinsolInterface::func_wrapper(N_Vector u, N_Vector fval, void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.func(*this_, u, fval);
      return 0;
    } catch(exception& e) {
      uerr() << "func failed: " << e.what() << endl;
      return 1;
    }
  }

  int KinsolInterface::djac_wrapper(long N, N_Vector u, N_Vector fu, DlsMat J,
                                 void *user_data, N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.djac(*this_, N, u, fu, J, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      uerr() << "djac failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolInterface::djac(KinsolMemory& m, long N, N_Vector u, N_Vector fu, DlsMat J,
                          N_Vector tmp1, N_Vector tmp2) const {
    // Evaluate jac_
    copy_n(m.iarg, n_in_, m.arg);
    m.arg[iin_] = NV_DATA_S(u);
    fill_n(m.res, n_out_+1, nullptr);
    m.res[0] = m.jac;
    calc_function(&m, "jac_f_z");

    // Get sparsity and non-zero elements
    const casadi_int* colind = sp_jac_.colind();
    casadi_int ncol = sp_jac_.size2();
    const casadi_int* row = sp_jac_.row();

    // Loop over columns
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        casadi_int rr = row[el];

        // Set the element
        DENSE_ELEM(J, rr, cc) = m.jac[el];
      }
    }
  }

  int KinsolInterface::bjac_wrapper(long N, long mupper, long mlower, N_Vector u,
                                            N_Vector fu, DlsMat J, void *user_data,
                                            N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.bjac(*this_, N, mupper, mlower, u, fu, J, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      uerr() << "bjac failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolInterface::bjac(KinsolMemory& m, long N, long mupper, long mlower, N_Vector u,
                          N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) const {
    // Evaluate jac_
    copy_n(m.iarg, n_in_, m.arg);
    m.arg[iin_] = NV_DATA_S(u);
    fill_n(m.res, n_out_+1, nullptr);
    m.res[0] = m.jac;
    calc_function(&m, "jac_f_z");

    // Get sparsity and non-zero elements
    const casadi_int* colind = sp_jac_.colind();
    casadi_int ncol = sp_jac_.size2();
    const casadi_int* row = sp_jac_.row();

    // Loop over cols
    for (casadi_int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (casadi_int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        casadi_int rr = row[el];

        // Set the element
        if (rr-cc>=-mupper && rr-cc<=mlower) {
          BAND_ELEM(J, rr, cc) = m.jac[el];
        }
      }
    }
  }

  int KinsolInterface::jtimes_wrapper(N_Vector v, N_Vector Jv, N_Vector u, int* new_u,
                                   void *user_data) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.jtimes(*this_, v, Jv, u, new_u);
      return 0;
    } catch(exception& e) {
      uerr() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolInterface::jtimes(KinsolMemory& m, N_Vector v, N_Vector Jv,
                            N_Vector u, int* new_u) const {
    // Evaluate f_fwd_
    copy_n(m.iarg, n_in_, m.arg);
    m.arg[iin_] = NV_DATA_S(u);
    m.arg[n_in_] = NV_DATA_S(v);
    m.res[0] = NV_DATA_S(Jv);
    jtimes_(m.arg, m.res, m.iw, m.w, 0);
  }

  int KinsolInterface::
  psetup_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                 void* user_data, N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.psetup(*this_, u, uscale, fval, fscale, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      uerr() << "psetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolInterface::
  psetup(KinsolMemory& m, N_Vector u, N_Vector uscale, N_Vector fval,
         N_Vector fscale, N_Vector tmp1, N_Vector tmp2) const {
    // Evaluate jac_
    copy_n(m.iarg, n_in_, m.arg);
    m.arg[iin_] = NV_DATA_S(u);
    fill_n(m.res, n_out_+1, nullptr);
    m.res[0] = m.jac;
    if (calc_function(&m, "jac_f_z")) casadi_error("Jacobian calculation failed");

    // Get sparsity and non-zero elements
    //const int* colind = sp_jac_.colind();
    //int ncol = sp_jac_.size2();
    //const int* row = sp_jac_.row();

    // Factorize the linear system
    if (linsol_.nfact(m.jac)) casadi_error("'nfact' failed");
  }

  int KinsolInterface::psolve_wrapper(N_Vector u, N_Vector uscale, N_Vector fval,
                                     N_Vector fscale, N_Vector v, void* user_data, N_Vector tmp) {
    try {
      casadi_assert_dev(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->self.psolve(*this_, u, uscale, fval, fscale, v, tmp);
      return 0;
    } catch(exception& e) {
      uerr() << "psolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolInterface::psolve(KinsolMemory& m, N_Vector u, N_Vector uscale, N_Vector fval,
                            N_Vector fscale, N_Vector v, N_Vector tmp) const {
    // Solve the factorized system
    if (linsol_.solve(m.jac, NV_DATA_S(v))) casadi_error("'solve' failed");
  }

  int KinsolInterface::lsetup(KINMem kin_mem) {
    try {
      auto m = to_mem(kin_mem->kin_lmem);
      auto& s = m->self;

      N_Vector u =  kin_mem->kin_uu;
      N_Vector uscale = kin_mem->kin_uscale;
      N_Vector fval = kin_mem->kin_fval;
      N_Vector fscale = kin_mem->kin_fscale;
      N_Vector tmp1 = kin_mem->kin_vtemp1;
      N_Vector tmp2 = kin_mem->kin_vtemp2;
      s.psetup(*m, u, uscale, fval, fscale, tmp1, tmp2);

      return 0;
    } catch(exception& e) {
      uerr() << "lsetup failed: " << e.what() << endl;;
      return -1;
    }
  }

  int KinsolInterface::
  lsolve(KINMem kin_mem, N_Vector x, N_Vector b, double *sJpnorm, double *sFdotJp) {
    try {
      auto m = to_mem(kin_mem->kin_lmem);
      auto& s = m->self;

      // Get vectors
      N_Vector u =  kin_mem->kin_uu;
      N_Vector uscale = kin_mem->kin_uscale;
      N_Vector fval = kin_mem->kin_fval;
      N_Vector fscale = kin_mem->kin_fscale;
      N_Vector tmp1 = kin_mem->kin_vtemp1;
      //N_Vector tmp2 = kin_mem->kin_vtemp2;

      // Solve the linear system
      N_VScale(1.0, b, x);
      s.psolve(*m, u, uscale, fval, fscale, x, tmp1);

      // Calculate residuals
      int flag = KINSpilsAtimes(kin_mem, x, b);
      if (flag) return flag;
      *sJpnorm = N_VWL2Norm(b, fscale);
      N_VProd(b, fscale, b);
      N_VProd(b, fscale, b);
      *sFdotJp = N_VDotProd(fval, b);

      return 0;
    } catch(exception& e) {
      uerr() << "lsolve failed: " << e.what() << endl;;
      return -1;
    }
  }

  void KinsolInterface::
  ehfun(int error_code, const char *module, const char *function,
                char *msg, void *eh_data) {
    try {
      auto m = to_mem(eh_data);
      auto& s = m->self;
      if (!s.disable_internal_warnings_) {
        uerr() << msg << endl;
      }
    } catch(exception& e) {
      uerr() << "ehfun failed: " << e.what() << endl;
    }
  }

  void KinsolInterface::kinsol_error(const string& module, int flag, bool fatal) const {
    // Get the error message
    const char *id, *msg;
    switch (flag) {
    case KIN_SUCCESS:
      id = "KIN_SUCCES";
      msg = "KINSol succeeded; the scaled norm of F(u) is less than fnormtol";
    break;
  case KIN_INITIAL_GUESS_OK:
    id = "KIN_INITIAL_GUESS_OK";
    msg = "The guess u = u0 satisfied the system F(u) = 0 within the tolerances specified.";
    break;
    case KIN_STEP_LT_STPTOL:
      id = "KIN_STEP_LT_STPTOL";
      msg = "KINSol stopped based on scaled step length. This "
        "means that the current iterate may be an approximate solution of the "
        "given nonlinear system, but it is also quite possible that the algorithm"
        " is 'stalled' (making insufficient progress) near an invalid solution, "
        "or that the scalar scsteptol is too large.";
      break;
    case KIN_MEM_NULL:
      id = "KIN_MEM_NULL";
      msg = "The kinsol memory block pointer was NULL.";
      break;
    case KIN_ILL_INPUT:
      id = "KIN_ILL_INPUT";
      msg = "An input parameter was invalid.";
      break;
    case KIN_NO_MALLOC:
      id = "KIN_NO_MALLOC";
      msg = "The kinsol memory was not allocated by a call to KINCreate.";
      break;
    case KIN_LINESEARCH_NONCONV:
      id = "KIN_LINESEARCH_NONCONV";
      msg = "The line search algorithm was unable to find "
        "an iterate sufficiently distinct from the current iterate, or could not"
        " find an iterate satisfying the sufficient decrease condition. Failure"
        " to satisfy the sufficient decrease condition could mean the current "
        "iterate is 'close' to an approximate solution of the given nonlinear "
        "system, the difference approximation of the matrix-vector product J(u)v"
        " is inaccurate, or the real scalar scsteptol is too large.";
      break;
    case KIN_MAXITER_REACHED:
      id = "KIN_MAXITER_REACHED";
      msg = "The maximum number of nonlinear iterations "
        "has been reached.";
      break;
    case KIN_MXNEWT_5X_EXCEEDED:
      id = "KIN_MXNEWT_5X_EXCEEDED";
      msg = "Five consecutive steps have been taken that "
        "satisfy the inequality  || D_u p ||_L2 > 0.99 mxnewtstep, where p "
        "denotes the current step and mxnewtstep is a scalar upper bound on the "
        "scaled step length. Such a failure may mean that || D_F F(u)||_L2 "
        "asymptotes from above to a positive value, or the real scalar "
        "mxnewtstep is too small.";
      break;
    case KIN_LINESEARCH_BCFAIL:
      id = "KIN_LINESEARCH_BCFAIL";
      msg = "The line search algorithm was unable to satisfy "
        "the “beta-condition” for MXNBCF +1 nonlinear iterations (not necessarily "
        "consecutive), which may indicate the algorithm is making poor progress.";
      break;
    case KIN_LINSOLV_NO_RECOVERY:
      id = "KIN_LINSOLV_NO_RECOVERY";
      msg = "The user-supplied routine psolve encountered a"
        " recoverable error, but the preconditioner is already current.";
      break;
    case KIN_LINIT_FAIL:
      id = "KIN_LINIT_FAIL";
      msg = "The linear solver initialization routine (linit) encountered an error.";
    break;
    case KIN_LSETUP_FAIL:
      id = "KIN_LSETUP_FAIL";
      msg = "The user-supplied routine pset (used to set up the "
        "preconditioner data) encountered an unrecoverable error.";
      break;
    case KIN_LSOLVE_FAIL:
      id = "KIN_LSOLVE_FAIL";
      msg = "Either the user-supplied routine psolve "
        "(used to to solve the preconditioned linear system) encountered an "
        "unrecoverable error, or the linear solver routine (lsolve) "
        "encountered an error condition.";
      break;
    case KIN_SYSFUNC_FAIL:
      id = "KIN_SYSFUNC_FAIL";
      msg = "The system function failed in an unrecoverable manner.";
      break;
    case KIN_FIRST_SYSFUNC_ERR:
      id = "KIN_FIRST_SYSFUNC_ERR";
      msg = "The system function failed recoverably at the first call.";
      break;
    case KIN_REPTD_SYSFUNC_ERR:
      id = "KIN_REPTD_SYSFUNC_ERR";
      msg = "The system function had repeated recoverable errors. "
        "No recovery is possible.";
      break;
    default:
      id = "N/A";
      msg = nullptr;
    }

    // Construct message
    stringstream ss;
    if (msg==nullptr) {
      ss << "Unknown " << (fatal? "error" : "warning") <<" (" << flag << ")"
        " from module \"" << module << "\".";
    } else {
      ss << "Module \"" << module << "\" returned flag \"" << id << "\"." << endl;
      ss << "The description of this flag is: " << endl;
      ss << "\"" << msg << "\"" << endl;
    }
    ss << "Consult KINSOL documentation for more information.";
    if (fatal) {
      casadi_error("nlpsol process failed. "
                   "Set 'error_on_fail' option to false to ignore this error. "
                   + ss.str());
    } else {
      casadi_warning(ss.str());
    }
  }

  KinsolMemory::KinsolMemory(const KinsolInterface& s) : self(s) {
    this->u = nullptr;
    this->mem = nullptr;
  }

  KinsolMemory::~KinsolMemory() {
    if (this->u) N_VDestroy_Serial(this->u);
    if (this->mem) KINFree(&this->mem);
  }

  int KinsolInterface::init_mem(void* mem) const {
    if (Rootfinder::init_mem(mem)) return 1;
    auto m = static_cast<KinsolMemory*>(mem);

    // Current solution
    m->u = N_VNew_Serial(n_);

    // Create KINSOL memory block
    m->mem = KINCreate();

    // KINSOL bugfix
    KINMem kin_mem = KINMem(m->mem);
    kin_mem->kin_inexact_ls = FALSE;

    // Set optional inputs
    int flag = KINSetUserData(m->mem, m);
    casadi_assert(flag==KIN_SUCCESS, "KINSetUserData");

    // Set error handler function
    flag = KINSetErrHandlerFn(m->mem, ehfun, m);
    casadi_assert(flag==KIN_SUCCESS, "KINSetErrHandlerFn");

    // Initialize KINSOL
    flag = KINInit(m->mem, func_wrapper, m->u);
    casadi_assert_dev(flag==KIN_SUCCESS);

    // Setting maximum number of Newton iterations
    flag = KINSetMaxNewtonStep(m->mem, max_iter_);
    casadi_assert_dev(flag==KIN_SUCCESS);

    // Set constraints
    if (!u_c_.empty()) {
      N_Vector domain  = N_VNew_Serial(n_);
      copy(u_c_.begin(), u_c_.end(), NV_DATA_S(domain));

      // Pass to KINSOL
      flag = KINSetConstraints(m->mem, domain);
      casadi_assert_dev(flag==KIN_SUCCESS);

      // Free the temporary vector
      N_VDestroy_Serial(domain);
    }

    switch (linear_solver_type_) {
    case KinsolInterface::DENSE:
      // Dense Jacobian
      flag = KINDense(m->mem, n_);
      casadi_assert(flag==KIN_SUCCESS, "KINDense");

      if (exact_jac_) {
        flag = KINDlsSetDenseJacFn(m->mem, djac_wrapper);
        casadi_assert(flag==KIN_SUCCESS, "KINDlsSetDenseJacFn");
      }
      break;
    case KinsolInterface::BANDED:
      // Banded Jacobian
      flag = KINBand(m->mem, n_, upper_bandwidth_, lower_bandwidth_);
      casadi_assert(flag==KIN_SUCCESS, "KINBand");

      if (exact_jac_) {
        flag = KINDlsSetBandJacFn(m->mem, bjac_wrapper);
        casadi_assert(flag==KIN_SUCCESS, "KINDlsBandJacFn");
      }
      break;
    case KinsolInterface::ITERATIVE:
      // Attach the sparse solver
      switch (iterative_solver_) {
      case KinsolInterface::GMRES:
        flag = KINSpgmr(m->mem, maxl_);
        casadi_assert(flag==KIN_SUCCESS, "KINSpgmr");
        break;
      case KinsolInterface::BCGSTAB:
        flag = KINSpbcg(m->mem, maxl_);
        casadi_assert(flag==KIN_SUCCESS, "KINSpbcg");
        break;
      case KinsolInterface::TFQMR:
        flag = KINSptfqmr(m->mem, maxl_);
        casadi_assert(flag==KIN_SUCCESS, "KINSptfqmr");
        break;
      }

      // Attach functions for Jacobian information
      if (exact_jac_) {
        flag = KINSpilsSetJacTimesVecFn(m->mem, jtimes_wrapper);
        casadi_assert(flag==KIN_SUCCESS, "KINSpilsSetJacTimesVecFn");
      }

      // Add a preconditioner
      if (use_preconditioner_) {
        flag = KINSpilsSetPreconditioner(m->mem, psetup_wrapper, psolve_wrapper);
        casadi_assert_dev(flag==KIN_SUCCESS);
      }
      break;
    case KinsolInterface::USER_DEFINED:
      // Set fields in the IDA memory
      KINMem kin_mem = KINMem(m->mem);
      kin_mem->kin_lmem   = m;
      kin_mem->kin_lsetup = lsetup;
      kin_mem->kin_lsolve = lsolve;
      kin_mem->kin_setupNonNull = TRUE;
      break;
    }

    // Set stop criterion
    flag = KINSetFuncNormTol(m->mem, abstol_);
    casadi_assert_dev(flag==KIN_SUCCESS);
    return 0;
  }

} // namespace casadi
