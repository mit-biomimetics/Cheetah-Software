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


#include "knitro_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_KNITRO_EXPORT
  casadi_register_nlpsol_knitro(Nlpsol::Plugin* plugin) {
    plugin->creator = KnitroInterface::creator;
    plugin->name = "knitro";
    plugin->doc = KnitroInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &KnitroInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_KNITRO_EXPORT casadi_load_nlpsol_knitro() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_knitro);
  }

  KnitroInterface::KnitroInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  KnitroInterface::~KnitroInterface() {
    clear_mem();
  }

  Options KnitroInterface::options_
  = {{&Nlpsol::options_},
     {{"knitro",
       {OT_DICT,
        "Options to be passed to KNITRO"}},
      {"detect_linear_constraints",
       {OT_BOOL,
        "Detect type of constraints"}},
      {"contype",
       {OT_INTVECTOR,
        "Type of constraint"}}
     }
  };

  void KnitroInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    bool detect_linear_constraints = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="knitro") {
        opts_ = op.second;
      } else if (op.first=="contype") {
        contype_ = op.second;
      } else if (op.first=="detect_linear_constraints") {
        detect_linear_constraints = op.second;
      }
    }

    casadi_assert(!detect_linear_constraints || contype_.empty(),
      "When specifying 'contype', set 'detect_linear_constraints' to false");

    // Type of constraints, general by default
    if (contype_.empty()) {
      contype_.resize(ng_, KTR_CONTYPE_GENERAL);
      if (detect_linear_constraints) {
        std::vector<bool> nl_g = oracle_.which_depends("x", {"g"}, 2, true);
        for (casadi_int i=0;i<ng_;++i)
          contype_[i] = nl_g[i] ? KTR_CONTYPE_GENERAL : KTR_CONTYPE_LINEAR;
      }
    }

    casadi_assert_dev(contype_.size()==ng_);

    // Setup NLP functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    Function gf_jg_fcn = create_function("nlp_gf_jg", {"x", "p"}, {"grad:f:x", "jac:g:x"});
    jacg_sp_ = gf_jg_fcn.sparsity_out(1);
    Function hess_l_fcn = create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                                  {"hess:gamma:x:x"},
                                  {{"gamma", {"f", "g"}}});
    hesslag_sp_ = hess_l_fcn.sparsity_out(0);

    // Allocate persistent memory
    alloc_w(nx_, true); // wlbx_
    alloc_w(nx_, true); // wubx_
    alloc_w(ng_, true); // wlbg_
    alloc_w(ng_, true); // wubg_
  }

  int KnitroInterface::init_mem(void* mem) const {
    return Nlpsol::init_mem(mem);
    //auto m = static_cast<KnitroMemory*>(mem);

    // Commented out since I have not found out how to change the bounds
    // Allocate KNITRO memory block
    /*  m.kc = KTR_new(); */
  }

  void KnitroInterface::set_work(void* mem, const double**& arg, double**& res,
                                 casadi_int*& iw, double*& w) const {
    auto m = static_cast<KnitroMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Copy inputs to temporary arrays
    m->wlbx = w; w += nx_;
    m->wubx = w; w += nx_;
    m->wlbg = w; w += ng_;
    m->wubg = w; w += ng_;
  }

  int casadi_KTR_puts(const char * const str, void * const userParams) {
    std::string s(str);
    uout() << s << std::flush;
    return s.size();
  }

  int KnitroInterface::solve(void* mem) const {
    auto m = static_cast<KnitroMemory*>(mem);

    // Allocate KNITRO memory block (move back to init!)
    casadi_assert_dev(m->kc==nullptr);
    m->kc = KTR_new_puts(casadi_KTR_puts, nullptr);
    casadi_assert_dev(m->kc!=nullptr);
    casadi_int status;

    // Jacobian sparsity
    vector<int> Jcol, Jrow;
    if (!jacg_sp_.is_null()) {
      assign_vector(jacg_sp_.get_col(), Jcol);
      assign_vector(jacg_sp_.get_row(), Jrow);
    }

    // Hessian sparsity
    casadi_int nnzH = hesslag_sp_.is_null() ? 0 : hesslag_sp_.nnz();
    vector<int> Hcol, Hrow;
    if (nnzH>0) {
      assign_vector(hesslag_sp_.get_col(), Hcol);
      assign_vector(hesslag_sp_.get_row(), Hrow);
      status = KTR_set_int_param_by_name(m->kc, "hessopt", KTR_HESSOPT_EXACT);
      casadi_assert(status==0, "KTR_set_int_param failed");
    } else {
      status = KTR_set_int_param_by_name(m->kc, "hessopt", KTR_HESSOPT_LBFGS);
      casadi_assert(status==0, "KTR_set_int_param failed");
    }

    // Pass user set options
    for (auto&& op : opts_) {
      int param_id;
      casadi_assert(KTR_get_param_id(m->kc, op.first.c_str(), &param_id)==0,
        "Unknown parameter '" + op.first + "'.");

      int param_type;
      casadi_assert(!KTR_get_param_type(m->kc, param_id, &param_type),
        "Error when setting option '" + op.first + "'.");

      switch (param_type) {
        case KTR_PARAMTYPE_INTEGER:
          casadi_assert(!KTR_set_int_param(m->kc, param_id, op.second),
            "Error when setting option '" + op.first + "'.");
          continue;
        case KTR_PARAMTYPE_FLOAT:
          casadi_assert(!KTR_set_double_param(m->kc, param_id, op.second),
            "Error when setting option '" + op.first + "'.");
          continue;
        case KTR_PARAMTYPE_STRING:
          {
            string str = op.second.to_string();
            casadi_assert(!KTR_set_char_param(m->kc, param_id, str.c_str()),
              "Error when setting option '" + op.first + "'.");
          }
          continue;
        default:
          casadi_error("Error when setting option '" + op.first + "'.");
      }
    }

    // "Correct" upper and lower bounds
    casadi_copy(m->lbx, nx_, m->wlbx);
    casadi_copy(m->ubx, nx_, m->wubx);
    casadi_copy(m->lbg, ng_, m->wlbg);
    casadi_copy(m->ubg, ng_, m->wubg);
    for (casadi_int i=0; i<nx_; ++i) if (isinf(m->wlbx[i])) m->wlbx[i] = -KTR_INFBOUND;
    for (casadi_int i=0; i<nx_; ++i) if (isinf(m->wubx[i])) m->wubx[i] =  KTR_INFBOUND;
    for (casadi_int i=0; i<ng_; ++i) if (isinf(m->wlbg[i])) m->wlbg[i] = -KTR_INFBOUND;
    for (casadi_int i=0; i<ng_; ++i) if (isinf(m->wubg[i])) m->wubg[i] =  KTR_INFBOUND;

    if (mi_) {
      // Convexity status of the objective function
      const casadi_int objFnType = KTR_FNTYPE_UNCERTAIN;

      // Types of variables
      vector<int> vtype;
      vtype.reserve(nx_);
      for (auto&& e : discrete_) {
        vtype.push_back(e ? KTR_VARTYPE_INTEGER : KTR_VARTYPE_CONTINUOUS);
      }

      // Convexity status of the constraint functions
      vector<int> ftype(ng_, KTR_FNTYPE_UNCERTAIN);

      // Intialize
      status =
      KTR_mip_init_problem(m->kc, nx_, KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
                           objFnType, get_ptr(vtype), m->wlbx, m->wubx,
                           ng_, get_ptr(contype_), get_ptr(ftype),
                           m->wlbg, m->wubg, Jcol.size(), get_ptr(Jcol), get_ptr(Jrow),
                           nnzH, get_ptr(Hrow), get_ptr(Hcol), m->x, nullptr);
      casadi_assert(status==0, "KTR_mip_init_problem failed");
    } else {
      status =
      KTR_init_problem(m->kc, nx_, KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
                       m->wlbx, m->wubx, ng_, get_ptr(contype_),
                       m->wlbg, m->wubg, Jcol.size(), get_ptr(Jcol), get_ptr(Jrow),
                       nnzH, get_ptr(Hrow), get_ptr(Hcol), m->x, nullptr); // initial lambda
      casadi_assert(status==0, "KTR_init_problem failed");
    }

    // Register callback functions
    status = KTR_set_func_callback(m->kc, &callback);
    casadi_assert(status==0, "KTR_set_func_callback failed");

    status = KTR_set_grad_callback(m->kc, &callback);
    casadi_assert(status==0, "KTR_set_grad_callbackfailed");

    if (nnzH>0) {
      status = KTR_set_hess_callback(m->kc, &callback);
      casadi_assert(status==0, "KTR_set_hess_callbackfailed");
    }

    // Lagrange multipliers
    vector<double> lambda(nx_+ng_);

    // Solve NLP
    double f;
    if (mi_) {
      status =
      KTR_mip_solve(m->kc, m->x, get_ptr(lambda), 0, &f,
                    nullptr, nullptr, nullptr, nullptr, nullptr, static_cast<void*>(m));

    } else {
      status =
      KTR_solve(m->kc, m->x, get_ptr(lambda), 0, &f,
                nullptr, nullptr, nullptr, nullptr, nullptr, static_cast<void*>(m));
    }
    m->return_status = return_codes(status);
    m->success = status==KTR_RC_OPTIMAL_OR_SATISFACTORY ||
                 status==KTR_RC_NEAR_OPT;

    // Output dual solution
    casadi_copy(get_ptr(lambda), ng_, m->lam_g);
    casadi_copy(get_ptr(lambda)+ng_, nx_, m->lam_x);

    // Output optimal cost
    m->f = f;

    // Calculate constraints
    if (m->g) {
      m->arg[0] = m->x;
      m->arg[1] = m->p;
      m->res[0] = nullptr;
      m->res[1] = m->g;
      calc_function(m, "nlp_fg");
    }

    // Free memory (move to destructor!)
    KTR_free(&m->kc);
    m->kc = nullptr;
    return 0;
  }

  int KnitroInterface::callback(const int evalRequestCode, const int n,
                              const int m, const int nnzJ,
                              const int nnzH, const double* const x,
                              const double* const lambda,
                              double* const obj, double* const c, double* const objGrad,
                              double* const jac, double* const hessian, double* const hessVector,
                              void *userParams) {
    try {
      // Get a pointer to the calling object
      auto m = static_cast<KnitroMemory*>(userParams);

      // Direct to the correct function
      switch (evalRequestCode) {
      case KTR_RC_EVALFC:
      m->arg[0] = x;
      m->arg[1] = m->p;
      m->res[0] = obj;
      m->res[1] = c;
      m->self.calc_function(m, "nlp_fg");
      break;
      case KTR_RC_EVALGA:
      m->arg[0] = x;
      m->arg[1] = m->p;
      m->res[0] = objGrad;
      m->res[1] = jac;
      m->self.calc_function(m, "nlp_gf_jg");
      break;
      case KTR_RC_EVALH:
        {
          double sigma = 1.;
          m->arg[0] = x;
          m->arg[1] = m->p;
          m->arg[2] = &sigma;
          m->arg[3] = lambda;
          m->res[0] = hessian;
          if (m->self.calc_function(m, "nlp_hess_l")) {
            casadi_error("calc_hess_l failed");
          }
        }
        break;
      default:
        casadi_error("KnitroInterface::callback: unknown method");
      }

      return 0;
    } catch(KeyboardInterruptException& ex) {
      return KTR_RC_USER_TERMINATION;
    } catch(exception& ex) {
      uerr() << "KnitroInterface::callback caught exception: "
                               << ex.what() << endl;
      return -1;
    }

  }

  const char* KnitroInterface::return_codes(int flag) {
    switch (flag) {
    case KTR_RC_OPTIMAL_OR_SATISFACTORY: return "KTR_RC_OPTIMAL_OR_SATISFACTORY";
    case KTR_RC_NEAR_OPT: return "KTR_RC_NEAR_OPT";
    case KTR_RC_FEAS_XTOL: return "KTR_RC_FEAS_XTOL";
    case KTR_RC_FEAS_NO_IMPROVE: return "KTR_RC_FEAS_NO_IMPROVE";
    case KTR_RC_FEAS_FTOL: return "KTR_RC_FEAS_FTOL";
    case KTR_RC_INFEASIBLE: return "KTR_RC_INFEASIBLE";
    case KTR_RC_INFEAS_XTOL: return "KTR_RC_INFEAS_XTOL";
    case KTR_RC_INFEAS_NO_IMPROVE: return "KTR_RC_INFEAS_NO_IMPROVE";
    case KTR_RC_INFEAS_MULTISTART: return "KTR_RC_INFEAS_MULTISTART";
    case KTR_RC_INFEAS_CON_BOUNDS: return "KTR_RC_INFEAS_CON_BOUNDS";
    case KTR_RC_INFEAS_VAR_BOUNDS: return "KTR_RC_INFEAS_VAR_BOUNDS";
    case KTR_RC_UNBOUNDED: return "KTR_RC_UNBOUNDED";
    case KTR_RC_ITER_LIMIT_FEAS: return "KTR_RC_ITER_LIMIT_FEAS";
    case KTR_RC_TIME_LIMIT_FEAS: return "KTR_RC_TIME_LIMIT_FEAS";
    case KTR_RC_FEVAL_LIMIT_FEAS: return "KTR_RC_FEVAL_LIMIT_FEAS";
    case KTR_RC_MIP_EXH_FEAS: return "KTR_RC_MIP_EXH_FEAS";
    case KTR_RC_MIP_TERM_FEAS: return "KTR_RC_MIP_TERM_FEAS";
    case KTR_RC_MIP_SOLVE_LIMIT_FEAS: return "KTR_RC_MIP_SOLVE_LIMIT_FEAS";
    case KTR_RC_MIP_NODE_LIMIT_FEAS: return "KTR_RC_MIP_NODE_LIMIT_FEAS";
    case KTR_RC_ITER_LIMIT_INFEAS: return "KTR_RC_ITER_LIMIT_INFEAS";
    case KTR_RC_TIME_LIMIT_INFEAS: return "KTR_RC_TIME_LIMIT_INFEAS";
    case KTR_RC_FEVAL_LIMIT_INFEAS: return "KTR_RC_FEVAL_LIMIT_INFEAS";
    case KTR_RC_MIP_EXH_INFEAS: return "KTR_RC_MIP_EXH_INFEAS";
    case KTR_RC_MIP_SOLVE_LIMIT_INFEAS: return "KTR_RC_MIP_SOLVE_LIMIT_INFEAS";
    case KTR_RC_MIP_NODE_LIMIT_INFEAS: return "KTR_RC_MIP_NODE_LIMIT_INFEAS";
    case KTR_RC_CALLBACK_ERR: return "KTR_RC_CALLBACK_ERR";
    case KTR_RC_LP_SOLVER_ERR: return "KTR_RC_LP_SOLVER_ERR";
    case KTR_RC_EVAL_ERR: return "KTR_RC_EVAL_ERR";
    case KTR_RC_OUT_OF_MEMORY: return "KTR_RC_OUT_OF_MEMORY";
    case KTR_RC_USER_TERMINATION: return "KTR_RC_USER_TERMINATION";
    case KTR_RC_OPEN_FILE_ERR: return "KTR_RC_OPEN_FILE_ERR";
    case KTR_RC_BAD_N_OR_F: return "KTR_RC_BAD_N_OR_F";
    case KTR_RC_BAD_CONSTRAINT: return "KTR_RC_BAD_CONSTRAINT";
    case KTR_RC_BAD_JACOBIAN: return "KTR_RC_BAD_JACOBIAN";
    case KTR_RC_BAD_HESSIAN: return "KTR_RC_BAD_HESSIAN";
    case KTR_RC_BAD_CON_INDEX: return "KTR_RC_BAD_CON_INDEX";
    case KTR_RC_BAD_JAC_INDEX: return "KTR_RC_BAD_JAC_INDEX";
    case KTR_RC_BAD_HESS_INDEX: return "KTR_RC_BAD_HESS_INDEX";
    case KTR_RC_BAD_CON_BOUNDS: return "KTR_RC_BAD_CON_BOUNDS";
    case KTR_RC_BAD_VAR_BOUNDS: return "KTR_RC_BAD_VAR_BOUNDS";
    case KTR_RC_ILLEGAL_CALL: return "KTR_RC_ILLEGAL_CALL";
    case KTR_RC_BAD_KCPTR: return "KTR_RC_BAD_KCPTR";
    case KTR_RC_NULL_POINTER: return "KTR_RC_NULL_POINTER";
    case KTR_RC_BAD_INIT_VALUE: return "KTR_RC_BAD_INIT_VALUE";
    case KTR_RC_NEWPOINT_HALT: return "KTR_RC_NEWPOINT_HALT";
    case KTR_RC_BAD_LICENSE: return "KTR_RC_BAD_LICENSE";
    case KTR_RC_BAD_PARAMINPUT: return "KTR_RC_BAD_PARAMINPUT";
    case KTR_RC_LINEAR_SOLVER_ERR: return "KTR_RC_LINEAR_SOLVER_ERR";
    case KTR_RC_DERIV_CHECK_FAILED: return "KTR_RC_DERIV_CHECK_FAILED";
    case KTR_RC_DERIV_CHECK_TERMINATE: return "KTR_RC_DERIV_CHECK_TERMINATE";
    case KTR_RC_INTERNAL_ERROR: return "KTR_RC_INTERNAL_ERROR";
    }
    return nullptr;
  }

  Dict KnitroInterface::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<KnitroMemory*>(mem);
    stats["return_status"] = m->return_status;

    return stats;
  }

  KnitroMemory::KnitroMemory(const KnitroInterface& self) : self(self) {
    this->kc = nullptr;
  }

  KnitroMemory::~KnitroMemory() {
    // Currently no persistent memory since KNITRO requires knowledge of nature of bounds
    // if (this->kc) {
    //   KTR_free(&this->kc);
    // }
  }

} // namespace casadi
