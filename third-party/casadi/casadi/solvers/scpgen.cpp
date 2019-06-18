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


#include "scpgen.hpp"
#include "casadi/core/core.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_SCPGEN_EXPORT
      casadi_register_nlpsol_scpgen(Nlpsol::Plugin* plugin) {
    plugin->creator = Scpgen::creator;
    plugin->name = "scpgen";
    plugin->doc = Scpgen::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Scpgen::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SCPGEN_EXPORT casadi_load_nlpsol_scpgen() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_scpgen);
  }

  Scpgen::Scpgen(const std::string& name, const Function& nlp) : Nlpsol(name, nlp) {
    casadi_warning("SCPgen is under development");
  }

  Scpgen::~Scpgen() {
    clear_mem();
  }

  Options Scpgen::options_
  = {{&Nlpsol::options_},
     {{"qpsol",
       {OT_STRING,
        "The QP solver to be used by the SQP method"}},
      {"qpsol_options",
       {OT_DICT,
        "Options to be passed to the QP solver"}},
      {"hessian_approximation",
       {OT_STRING,
        "gauss-newton|exact"}},
      {"max_iter",
       {OT_INT,
        "Maximum number of SQP iterations"}},
      {"max_iter_ls",
       {OT_INT,
        "Maximum number of linesearch iterations"}},
      {"tol_pr",
       {OT_DOUBLE,
        "Stopping criterion for primal infeasibility"}},
      {"tol_du",
       {OT_DOUBLE,
        "Stopping criterion for dual infeasability"}},
      {"tol_reg",
       {OT_DOUBLE,
        "Stopping criterion for regularization"}},
      {"tol_pr_step",
       {OT_DOUBLE,
        "Stopping criterion for the step size"}},
      {"c1",
       {OT_DOUBLE,
        "Armijo condition, coefficient of decrease in merit"}},
      {"beta",
       {OT_DOUBLE,
        "Line-search parameter, restoration factor of stepsize"}},
      {"merit_memsize",
       {OT_INT,
        "Size of memory to store history of merit function values"}},
      {"merit_start",
       {OT_DOUBLE,
        "Lower bound for the merit function parameter"}},
      {"lbfgs_memory",
       {OT_INT,
        "Size of L-BFGS memory."}},
      {"regularize",
       {OT_BOOL,
        "Automatic regularization of Lagrange Hessian."}},
      {"print_header",
       {OT_BOOL,
        "Print the header with problem statistics"}},
      {"codegen",
       {OT_BOOL,
        "C-code generation"}},
      {"reg_threshold",
       {OT_DOUBLE,
        "Threshold for the regularization."}},
      {"name_x",
       {OT_STRINGVECTOR,
        "Names of the variables."}},
      {"print_x",
       {OT_INTVECTOR,
        "Which variables to print."}}
     }
  };

  void Scpgen::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Default options
    max_iter_ = 50;
    max_iter_ls_ = 1;
    c1_ = 1e-4;
    beta_ = 0.8;
    lbfgs_memory_ = 10;
    tol_pr_ = 1e-6;
    tol_du_ = 1e-6;
    tol_reg_ = 1e-11;
    regularize_ = false;
    codegen_ = false;
    reg_threshold_ = 1e-8;
    tol_pr_step_ = 1e-6;
    merit_memsize_ = 4;
    merit_start_ = 1e-8;
    string hessian_approximation = "exact";
    string qpsol_plugin;
    Dict qpsol_options;
    print_header_ = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="max_iter_ls") {
        max_iter_ls_ = op.second;
      } else if (op.first=="c1") {
        c1_ = op.second;
      } else if (op.first=="beta") {
        beta_ = op.second;
      } else if (op.first=="lbfgs_memory") {
        lbfgs_memory_ = op.second;
      } else if (op.first=="tol_pr") {
        tol_pr_ = op.second;
      } else if (op.first=="tol_du") {
        tol_du_ = op.second;
      } else if (op.first=="tol_reg") {
        tol_reg_ = op.second;
      } else if (op.first=="regularize") {
        regularize_ = op.second;
      } else if (op.first=="codegen") {
        codegen_ = op.second;
      } else if (op.first=="reg_threshold") {
        reg_threshold_ = op.second;
      } else if (op.first=="tol_pr_step") {
        tol_pr_step_ = op.second;
      } else if (op.first=="merit_memsize") {
        merit_memsize_ = op.second;
      } else if (op.first=="merit_start") {
        merit_start_ = op.second;
      } else if (op.first=="hessian_approximation") {
        hessian_approximation = op.second.to_string();
      } else if (op.first=="name_x") {
        name_x_ = op.second;
      } else if (op.first=="print_x") {
        print_x_ = op.second;
      } else if (op.first=="qpsol") {
        qpsol_plugin = op.second.to_string();
      } else if (op.first=="qpsol_options") {
        qpsol_options = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      }
    }

    // Gauss-Newton Hessian?
    gauss_newton_ = hessian_approximation == "gauss-newton";

    // Name the components
    if (name_x_.empty()) {
      name_x_.resize(nx_);
      for (casadi_int i=0; i<nx_; ++i) {
        stringstream ss;
        ss << "x" << i;
        name_x_[i] = ss.str();
      }
    } else {
      casadi_assert_dev(name_x_.size()==nx_);
    }

    // Generate lifting functions
    Function fg = oracle();
    Function vdef_fcn, vinit_fcn;
    fg.generate_lifted(vdef_fcn, vinit_fcn);
    vinit_fcn_ = vinit_fcn;
    alloc(vinit_fcn_);

    // Extract the expressions
    vector<MX> vdef_in = vdef_fcn.mx_in();
    vector<MX> vdef_out = vdef_fcn(vdef_in);

    // Get the dimensions
    MX x = vdef_in.at(0);
    MX p = vdef_in.at(1);
    v_.resize(vdef_in.size()-2);
    for (casadi_int i=0; i<v_.size(); ++i) {
      v_[i].v = vdef_in.at(i+2);
      v_[i].v_def = vdef_out.at(i+2);
      v_[i].n = v_[i].v.nnz();
    }

    // Scalar objective function
    MX f;

    // Multipliers
    MX g_lam;

    // Definition of the lifted dual variables
    MX p_defL, gL_defL;

    if (gauss_newton_) {
      // Least square objective
      f = dot(vdef_out[0], vdef_out[0])/2;
      gL_defL = vdef_out[0];
      ngn_ = gL_defL.nnz();
    } else {
      // Scalar objective function
      f = vdef_out[0];

      // Lagrange multipliers corresponding to the definition of the dependent variables
      stringstream ss;
      casadi_int i=0;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        ss.str(string());
        ss << "lam_x" << i++;
        it->v_lam = MX::sym(ss.str(), it->v.sparsity());
      }

      // Lagrange multipliers for the nonlinear constraints
      g_lam = MX::sym("g_lam", ng_);

      if (verbose_) {
        uout() << "Allocated intermediate variables." << endl;
      }

      // Adjoint sweep to get the definitions of the lifted dual variables
      // (Equation 3.8 in Albersmeyer2010)
      vector<vector<MX> > aseed(1), asens(1);
      aseed[0].push_back(1.0);
      aseed[0].push_back(g_lam);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        aseed[0].push_back(it->v_lam);
      }
      vdef_fcn->call_reverse(vdef_in, vdef_out, aseed, asens, true, false);
      i=0;

      gL_defL = asens[0].at(i++);
      if (gL_defL.is_null()) gL_defL = MX::zeros(x.sparsity()); // Needed?

      p_defL = asens[0].at(i++);
      if (p_defL.is_null()) p_defL = MX::zeros(p.sparsity()); // Needed?

      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        it->v_defL = asens[0].at(i++);
        if (it->v_defL.is_null()) {
          it->v_defL = MX::zeros(it->v.sparsity());
        }
      }

      if (verbose_) {
        uout() << "Generated the gradient of the Lagrangian." << endl;
      }
    }

    // Residual function

    // Inputs
    vector<MX> res_fcn_in;
    casadi_int n=0;
    res_fcn_in.push_back(x);             res_x_ = n++;
    res_fcn_in.push_back(p);             res_p_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_in.push_back(it->v);        it->res_var = n++;
    }
    if (!gauss_newton_) {
      res_fcn_in.push_back(g_lam);        res_g_lam_ = n++;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        res_fcn_in.push_back(it->v_lam);  it->res_lam = n++;
      }
    }

    // Outputs
    vector<MX> res_fcn_out;
    n=0;
    res_fcn_out.push_back(f);                              res_f_ = n++;
    res_fcn_out.push_back(gL_defL);                        res_gl_ = n++;
    res_fcn_out.push_back(vdef_out[1]);                    res_g_ = n++;
    res_fcn_out.push_back(p_defL);                         res_p_d_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_out.push_back(it->v_def - it->v);             it->res_d = n++;
    }

    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        res_fcn_out.push_back(it->v_defL - it->v_lam);     it->res_lam_d = n++;
      }
    }

    // Generate function
    Function res_fcn("res_fcn", res_fcn_in, res_fcn_out);
    if (verbose_) {
      uout() << "Generated residual function ( " << res_fcn.n_nodes() << " nodes)." << endl;
    }

    // Declare difference vector d and substitute out p and v
    stringstream ss;
    casadi_int i=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      ss.str(string());
      ss << "d" << i++;
      it->d = MX::sym(ss.str(), it->v.sparsity());
      it->d_def = it->v_def - it->d;
    }

    // Declare difference vector lam_d and substitute out lam
    if (!gauss_newton_) {
      casadi_int i=0;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        ss.str(string());
        ss << "d_lam" << i++;
        it->d_lam = MX::sym(ss.str(), it->v.sparsity());
        it->d_defL = it->v_defL - it->d_lam;
      }
    }

    // Variables to be substituted and their definitions
    vector<MX> svar, sdef;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      svar.push_back(it->v);
      sdef.push_back(it->d_def);
    }
    if (!gauss_newton_) {
      for (vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it) {
        svar.push_back(it->v_lam);
        sdef.push_back(it->d_defL);
      }
    }

    vector<MX> ex(4);
    ex[0] = f;
    ex[1] = vdef_out[1];
    ex[2] = gL_defL;
    ex[3] = p_defL;

    substitute_inplace(svar, sdef, ex, false);
    i=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      it->d_def = sdef[i++];
    }
    if (!gauss_newton_) {
      for (vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it) {
        it->d_defL = sdef[i++];
      }
    }

    MX f_z = ex[0];
    MX g_z = ex[1];
    MX gL_z = ex[2];
    MX p_z = ex[3];

    // Modified function inputs
    vector<MX> mfcn_in;
    n=0;
    mfcn_in.push_back(p);                               mod_p_ = n++;
    mfcn_in.push_back(x);                               mod_x_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_in.push_back(it->d);                          it->mod_var = n++;
    }

    // Modified function outputs
    n=0;
    vector<MX> mfcn_out;
    mfcn_out.push_back(g_z);                             mod_g_ = n++;

    // Add multipliers to function inputs
    if (!gauss_newton_) {
      n = mfcn_in.size();
      mfcn_in.push_back(g_lam);                          mod_g_lam_ = n++;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        mfcn_in.push_back(it->d_lam);                    it->mod_lam = n++;
      }
    }

    // Add gradient of the Lagrangian
    n = mfcn_out.size();
    mfcn_out.push_back(f_z);                             mod_f_ = n++;
    mfcn_out.push_back(gL_z);                            mod_gl_ = n++;

    // Jacobian of the constraints
    MX jac = MX::jacobian(mfcn_out[mod_g_], mfcn_in[mod_x_]);
    if (verbose_) casadi_message("Formed Jacobian of the constraints.");

    // Hessian of the Lagrangian
    MX hes = MX::jacobian(mfcn_out[mod_gl_], mfcn_in[mod_x_], {{"symmetric", !gauss_newton_}});
    if (gauss_newton_) {
      if (verbose_) casadi_message("Formed square root of Gauss-Newton Hessian.");
    } else {
      if (verbose_) casadi_message("Formed Hessian of the Lagrangian.");
    }

    // Matrices in the reduced QP
    n=0;
    vector<MX> mat_out;
    mat_out.push_back(jac);                             mat_jac_ = n++;
    mat_out.push_back(hes);                             mat_hes_ = n++;
    Function mat_fcn("mat_fcn", mfcn_in, mat_out);

    // Definition of intermediate variables
    n = mfcn_out.size();
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_out.push_back(it->d_def);         it->mod_def = n++;
      if (!gauss_newton_) {
        mfcn_out.push_back(it->d_defL);      it->mod_defL = n++;
      }
    }

    // Modifier function
    Function mfcn("mfcn", mfcn_in, mfcn_out);

    // Directional derivative of Z
    vector<vector<MX> > mfcn_fwdSeed(1, mfcn_in), mfcn_fwdSens(1, mfcn_out);

    // Linearization in the d-direction (see Equation (2.12) in Alberspeyer2010)
    fill(mfcn_fwdSeed[0].begin(), mfcn_fwdSeed[0].end(), MX());
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_fwdSeed[0][it->mod_var] = it->d;
      if (!gauss_newton_) {
        mfcn_fwdSeed[0][it->mod_lam] = it->d_lam;
      }
    }
    mfcn->call_forward(mfcn_in, mfcn_out, mfcn_fwdSeed, mfcn_fwdSens, true, false);

    // Vector(s) b in Lifted Newton
    MX b_gf = densify(mfcn_fwdSens[0][mod_gl_]);
    MX b_g = densify(mfcn_fwdSens[0][mod_g_]);

    // Tangent function
    vector<MX> vec_fcn_out;
    n=0;
    vec_fcn_out.push_back(b_gf);                              vec_gf_ = n++;
    vec_fcn_out.push_back(b_g);                               vec_g_ = n++;
    casadi_assert_dev(n==vec_fcn_out.size());

    Function vec_fcn("vec_fcn", mfcn_in, vec_fcn_out);
    if (verbose_) {
      uout() << "Generated linearization function ( " << vec_fcn.n_nodes()
           << " nodes)." << endl;
    }

    // Expression a + A*du in Lifted Newton (Section 2.1 in Alberspeyer2010)
    MX du = MX::sym("du", nx_);   // Step in u
    MX g_dlam;               // Step lambda_g
    if (!gauss_newton_) {
      g_dlam = MX::sym("g_dlam", g_lam.sparsity());
    }

    // Interpret the Jacobian-vector multiplication as a forward directional derivative
    fill(mfcn_fwdSeed[0].begin(), mfcn_fwdSeed[0].end(), MX());
    mfcn_fwdSeed[0][mod_x_] = du;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_fwdSeed[0][it->mod_var] = -it->d;
    }
    if (!gauss_newton_) {
      mfcn_fwdSeed[0][mod_g_lam_] = g_dlam;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        mfcn_fwdSeed[0][it->mod_lam] = -it->d_lam;
      }
    }
    mfcn->call_forward(mfcn_in, mfcn_out, mfcn_fwdSeed, mfcn_fwdSens, true, false);

    // Step expansion function inputs
    n = mfcn_in.size();
    mfcn_in.push_back(du);                                 mod_du_ = n++;
    if (!gauss_newton_) {
      mfcn_in.push_back(g_dlam);                           mod_dlam_g_ = n++;
    }

    // Step expansion function outputs
    vector<MX> exp_fcn_out;
    n=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_def]); it->exp_def = n++;
    }

    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_defL]); it->exp_defL = n++;
      }
    }

    // Step expansion function
    Function exp_fcn("exp_fcn", mfcn_in, exp_fcn_out);
    if (verbose_) {
      uout() << "Generated step expansion function ( " << exp_fcn.n_nodes() << " nodes)."
           << endl;
    }

    // Generate c code and load as DLL
    if (codegen_) {
      // Name of temporary file
      string cname = temporary_file("tmp_casadi_scpgen", ".c");

      // Codegen the functions
      CodeGenerator gen(cname);
      gen.add(res_fcn);
      gen.add(mat_fcn);
      gen.add(vec_fcn);
      gen.add(exp_fcn);

      // Generate code
      if (verbose_) {
        uout() << "Generating \"" << cname << "\""  << endl;
      }
      string name = cname.substr(0, cname.find_first_of('.'));
      gen.generate();

      // Complile and run
      if (verbose_) {
        uout() << "Starting compilation"  << endl;
      }
      time_t time1 = time(nullptr);
      compiler_ = Importer(cname, compilerplugin_, jit_options_);
      time_t time2 = time(nullptr);
      double comp_time = difftime(time2, time1);
      if (verbose_) {
        uout() << "Compilation completed after " << comp_time << " s."  << endl;
      }

      // Load the generated code
      res_fcn_ = external("res_fcn", compiler_);
      mat_fcn_ = external("mat_fcn", compiler_);
      vec_fcn_ = external("vec_fcn", compiler_);
      exp_fcn_ = external("exp_fcn", compiler_);
    } else {
      mat_fcn_ = mat_fcn;
      res_fcn_ = res_fcn;
      vec_fcn_ = vec_fcn;
      exp_fcn_ = exp_fcn;
    }

    // Allocate a QP solver
    spL_ = mat_fcn_.sparsity_out(mat_hes_);
    spH_ = mtimes(spL_.T(), spL_);
    spA_ = mat_fcn_.sparsity_out(mat_jac_);
    casadi_assert(!qpsol_plugin.empty(), "'qpsol' option has not been set");
    qpsol_ = conic("qpsol", qpsol_plugin, {{"h", spH_}, {"a", spA_}},
                   qpsol_options);
    if (verbose_) {
      uout() << "Allocated QP solver." << endl;
    }

    if (verbose_) {
      uout() << "NLP preparation completed" << endl;
    }

    // Header
    if (print_header_) {
      uout() << "-------------------------------------------" << endl;
      uout() << "This is casadi::SCPgen." << endl;
      if (gauss_newton_) {
        uout() << "Using Gauss-Newton Hessian" << endl;
      } else {
        uout() << "Using exact Hessian" << endl;
      }

      // Count the total number of variables
      casadi_int n_lifted = 0;
      for (vector<Var>::const_iterator i=v_.begin(); i!=v_.end(); ++i) {
        n_lifted += i->n;
      }

      uout()
        << endl
        << "Number of reduced variables:               " << setw(9) << nx_ << endl
        << "Number of reduced constraints:             " << setw(9) << ng_ << endl
        << "Number of lifted variables/constraints:    " << setw(9) << n_lifted << endl
        << "Number of parameters:                      " << setw(9) << np_ << endl
        << "Total number of variables:                 " << setw(9) << (nx_+n_lifted) << endl
        << "Total number of constraints:               " << setw(9) << (ng_+n_lifted) << endl
        << endl;

      uout()
        << "Iteration options:" << endl
        << "{ \"max_iter\":" << max_iter_ << ", "
        << "\"max_iter_ls\":" << max_iter_ls_ << ", "
        << "\"c1\":" << c1_ << ", "
        << "\"beta\":" << beta_ << ", "
        << "\"merit_memsize\":" << merit_memsize_ << ", "
        << "\"merit_start\":" << merit_start_ << ", "
        << "\"regularize\":" << regularize_ << ", "
        << endl << "  "
        << "\"tol_pr\":" << tol_pr_ << ", "
        << "\"tol_du\":" << tol_du_ << ", "
        << "\"tol_reg\":" << tol_reg_ << ", "
        << "\"reg_threshold\":" << reg_threshold_ << "}" << endl
        << endl;
    }

    // Allocate memory, nonlfted problem
    alloc_w(ng_, true); // gk_
    alloc_w(nx_, true); // dxk_
    alloc_w(nx_, true); // lam_xk_
    alloc_w(nx_, true); // dlam_xk_
    alloc_w(ng_, true); // dlam_gk_
    alloc_w(nx_, true); // gfk_
    alloc_w(nx_, true); // gL_
    if (gauss_newton_) {
      alloc_w(ngn_, true); // b_gn_
    }

    // Allocate memory, lifted problem
    for (casadi_int i=0; i<v_.size(); ++i) {
      casadi_int n = v_[i].n;
      alloc_w(n, true); // dx
      alloc_w(n, true); // x0
      alloc_w(n, true); // x
      alloc_w(n, true); // res
      if (!gauss_newton_) {
        alloc_w(n, true); // lam
        alloc_w(n, true); // dlam
        alloc_w(n, true); // resL
      }
    }

    // Allocate QP
    alloc_w(spH_.nnz(), true); // qpH_
    alloc_w(spA_.nnz(), true); // qpA_
    alloc_w(ng_, true); // qpB_
    alloc_w(nx_, true); // qpH_times_du_
    if (gauss_newton_) {
      alloc_w(spL_.nnz(), true); // qpL_
      alloc_w(ngn_, true); // qpG_
    } else {
      alloc_w(nx_, true); // qpG_
    }

    // QP solver
    alloc_w(nx_, true); // qp_lbx_
    alloc_w(nx_, true); // qp_ubx_
    alloc_w(ng_, true); // qp_lba_
    alloc_w(ng_, true); // qp_uba_

    // Line-search memory
    alloc_w(merit_memsize_, true);

    // Temporary work vectors
    alloc(mat_fcn_);
    alloc(res_fcn_);
    alloc(vec_fcn_);
    alloc(exp_fcn_);
    if (gauss_newton_) {
      alloc_w(ngn_); // casadi_mul to get GN Hessian
    }
    alloc(qpsol_);
  }

  int Scpgen::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<ScpgenMemory*>(mem);

    // Lifted memory
    m->lifted_mem.resize(v_.size());
    for (casadi_int i=0; i<v_.size(); ++i) {
      m->lifted_mem[i].n = v_[i].n;
    }

    return 0;
  }

  void Scpgen::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<ScpgenMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Get work vectors, nonlifted problem
    m->gk = w; w += ng_;
    m->dxk = w; w += nx_;
    m->lam_xk = w; w += nx_;
    m->dlam_xk = w; w += nx_;
    m->dlam_gk = w; w += ng_;
    m->gfk = w; w += nx_;
    m->gL = w; w += nx_;
    if (gauss_newton_) {
      m->b_gn = w; w += ngn_;
    }

    // Get work vectors, lifted problem
    for (auto&& v : m->lifted_mem) {
      v.dx = w; w += v.n;
      v.x0 = w; w += v.n;
      v.x = w; w += v.n;
      v.res = w; w += v.n;
      if (!gauss_newton_) {
        v.lam = w; w += v.n;
        v.dlam = w; w += v.n;
        v.resL = w; w += v.n;
      }
    }

    // QP
    m->qpH = w; w += spH_.nnz();
    m->qpA = w; w += spA_.nnz();
    m->qpB = w; w += ng_;
    if (gauss_newton_) {
      m->qpL = w; w += spL_.nnz();
      m->qpG = w; w += ngn_;
    } else {
      m->qpG = w; w += nx_;
    }
    m->qpH_times_du = w; w += nx_;

    // QP solver
    m->qp_lbx = w; w += nx_;
    m->qp_ubx = w; w += nx_;
    m->qp_lba = w; w += ng_;
    m->qp_uba = w; w += ng_;

    // merit_mem
    m->merit_mem = w; w += merit_memsize_;

    // Residual
    for (auto&& v : m->lifted_mem) casadi_fill(v.res, v.n, 0.);
    if (!gauss_newton_) {
      for (auto&& v : m->lifted_mem) casadi_fill(v.resL, v.n, 0.);
    }
  }

  int Scpgen::solve(void* mem) const {
    auto m = static_cast<ScpgenMemory*>(mem);

    if (v_.size()>0) {
      // Initialize lifted variables using the generated function
      fill_n(m->arg, vinit_fcn_.n_in(), nullptr);
      m->arg[0] = m->x;
      m->arg[1] = m->p;
      fill_n(m->res, vinit_fcn_.n_out(), nullptr);
      for (casadi_int i=0; i<v_.size(); ++i) {
        m->res[i] = m->lifted_mem[i].x0;
      }
      vinit_fcn_(m->arg, m->res, m->iw, m->w, 0);
    }
    if (verbose_) {
      uout() << "Passed initial guess" << endl;
    }

    // Reset dual guess
    casadi_fill(m->dlam_gk, ng_, 0.);
    casadi_fill(m->lam_xk, nx_, 0.);
    casadi_fill(m->dlam_xk, nx_, 0.);
    if (!gauss_newton_) {
      for (auto&& v : m->lifted_mem) {
        casadi_fill(v.lam, v.n, 0.);
        casadi_fill(v.dlam, v.n, 0.);
      }
    }

    // Reset line-search
    m->merit_ind = 0;

    // Current guess for the primal solution
    for (auto&& v : m->lifted_mem) {
      casadi_copy(v.x0, v.n, v.x);
    }

    // Get current time and reset timers
    double time1 = clock();
    m->t_eval_mat = m->t_eval_res = m->t_eval_vec = m->t_eval_exp = m->t_solve_qp = 0;

    // Initial evaluation of the residual function
    eval_res(m);

    // Number of SQP iterations
    m->iter_count = 0;

    // Reset last step-size
    m->pr_step = 0;
    m->du_step = 0;

    // Reset line-search
    casadi_int ls_iter = 0;
    bool ls_success = true;

    // Reset regularization
    m->reg = 0;

    // Reset iteration message
    m->iteration_note = nullptr;

    // MAIN OPTIMZATION LOOP
    while (true) {

      // Evaluate the vectors in the condensed QP
      eval_vec(m);

      // Evaluate the matrices in the condensed QP
      eval_mat(m);

      // 1-norm of the primal infeasibility
      double pr_inf = primalInfeasibility(m);

      // 1-norm of the dual infeasibility
      double du_inf = dualInfeasibility(m);

      // Print header occasionally
      if (m->iter_count % 10 == 0) printIteration(m, uout());

      // Printing information about the actual iterate
      printIteration(m, uout(), m->iter_count, m->f, pr_inf, du_inf, m->reg,
                     ls_iter, ls_success);

      // Checking convergence criteria
      bool converged = pr_inf <= tol_pr_ && m->pr_step <= tol_pr_step_ && m->reg <= tol_reg_;
      converged = converged && du_inf <= tol_du_;
      if (converged) {
        uout() << endl << "casadi::SCPgen: Convergence achieved after "
                  << m->iter_count << " iterations." << endl;
        break;
      }

      if (m->iter_count >= max_iter_) {
        uout() << endl;
        uout() << "casadi::SCPgen: Maximum number of iterations reached." << endl;
        break;
      }

      // Check if not-a-number
      if (m->f!=m->f || m->pr_step != m->pr_step || pr_inf != pr_inf) {
        uout() << "casadi::SCPgen: Aborted, nan detected" << endl;
        break;
      }

      // Start a new iteration
      m->iter_count++;

      // Regularize the QP
      if (regularize_) {
        regularize(m);
      }

      // Solve the condensed QP
      solve_qp(m);

      // Expand the step
      eval_exp(m);

      // Line-search to take the step
      line_search(m, ls_iter, ls_success);
    }

    double time2 = clock();
    m->t_mainloop = (time2-time1)/CLOCKS_PER_SEC;

    // Store optimal value
    uout() << "optimal cost = " << m->f << endl;

    // Save results to outputs
    casadi_copy(m->lam_xk, nx_, m->lam_x);
    casadi_copy(m->gk, ng_, m->g);

    // Write timers
    if (print_time_) {
      uout() << endl;
      uout() << "time spent in eval_mat:    " << setw(9) << m->t_eval_mat << " s." << endl;
      uout() << "time spent in eval_res:    " << setw(9) << m->t_eval_res << " s." << endl;
      uout() << "time spent in eval_vec:    " << setw(9) << m->t_eval_vec << " s." << endl;
      uout() << "time spent in eval_exp:    " << setw(9) << m->t_eval_exp << " s." << endl;
      uout() << "time spent in solve_qp:    " << setw(9) << m->t_solve_qp << " s." << endl;
      uout() << "time spent in main loop:   " << setw(9) << m->t_mainloop << " s." << endl;
    }

    uout() << endl;
    return 0;
  }

  double Scpgen::primalInfeasibility(ScpgenMemory* m) const {
    // L1-norm of the primal infeasibility
    double pr_inf = 0;

    // Simple bounds
    pr_inf += casadi_sum_viol(nx_, m->x, m->lbx, m->ubx);

    // Lifted variables
    for (auto&& v : m->lifted_mem) pr_inf += casadi_norm_1(v.n, v.res);

    // Nonlinear bounds
    pr_inf += casadi_sum_viol(ng_, m->gk, m->lbg, m->ubg);

    return pr_inf;
  }

  double Scpgen::dualInfeasibility(ScpgenMemory* m) const {
    // L1-norm of the dual infeasibility
    return casadi_norm_1(nx_, m->gL);
  }

  void Scpgen::printIteration(ScpgenMemory* m, std::ostream &stream) const {
    stream << setw(4)  << "iter";
    stream << setw(14) << "objective";
    stream << setw(11) << "inf_pr";
    stream << setw(11) << "inf_du";
    stream << setw(11) << "pr_step";
    stream << setw(11) << "du_step";
    stream << setw(8) << "lg(rg)";
    stream << setw(3) << "ls";
    stream << ' ';

    // Print variables
    for (vector<casadi_int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i) {
      stream << setw(9) << name_x_.at(*i);
    }

    stream << endl;
    stream.unsetf(std::ios::floatfield);
  }

  void Scpgen::printIteration(ScpgenMemory* m, std::ostream &stream, casadi_int iter, double obj,
                              double pr_inf, double du_inf, double rg, casadi_int ls_trials,
                              bool ls_success) const {
    stream << setw(4) << iter;
    stream << scientific;
    stream << setw(14) << setprecision(6) << obj;
    stream << setw(11) << setprecision(2) << pr_inf;
    stream << setw(11);
    stream << setprecision(2) << du_inf;
    stream << setw(11) << setprecision(2) << m->pr_step;
    stream << setw(11);
    stream << setprecision(2) << m->du_step;
    stream << fixed;
    if (rg>0) {
      stream << setw(8) << setprecision(2) << log10(rg);
    } else {
      stream << setw(8) << "-";
    }
    stream << setw(3) << ls_trials;
    stream << (ls_success ? ' ' : 'F');

    // Print variables
    for (vector<casadi_int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i) {
      stream << setw(9) << setprecision(4) << m->x[*i];
    }

    // Print note
    if (m->iteration_note) {
      stream << "   " << m->iteration_note;
      m->iteration_note = nullptr;
    }

    stream.unsetf(std::ios::floatfield);
    stream << endl;
  }

  void Scpgen::eval_mat(ScpgenMemory* m) const {
    // Get current time
    double time1 = clock();

    // Inputs
    fill_n(m->arg, mat_fcn_.n_in(), nullptr);
    m->arg[mod_p_] = m->p; // Parameters
    m->arg[mod_x_] = m->x; // Primal step/variables
    for (size_t i=0; i<v_.size(); ++i) {
      m->arg[v_[i].mod_var] = m->lifted_mem[i].res;
    }
    if (!gauss_newton_) { // Dual steps/variables
      m->arg[mod_g_lam_] = m->lam_g;
      for (size_t i=0; i<v_.size(); ++i) {
        m->arg[v_[i].mod_lam] = m->lifted_mem[i].resL;
      }
    }

    // Outputs
    fill_n(m->res, mat_fcn_.n_out(), nullptr);
    m->res[mat_jac_] = m->qpA; // Condensed Jacobian
    m->res[mat_hes_] = gauss_newton_ ? m->qpL : m->qpH; // Condensed Hessian

    // Calculate condensed QP matrices
    mat_fcn_(m->arg, m->res, m->iw, m->w, 0);

    if (gauss_newton_) {
      // Gauss-Newton Hessian
      casadi_fill(m->qpH, spH_.nnz(), 0.);
      casadi_mtimes(m->qpL, spL_, m->qpL, spL_, m->qpH, spH_, m->w, true);

      // Gradient of the objective in Gauss-Newton
      casadi_fill(m->gfk, nx_, 0.);
      casadi_mv(m->qpL, spL_, m->b_gn, m->gfk, true);
    }

    // Calculate the gradient of the lagrangian
    casadi_copy(m->gfk, nx_, m->gL);
    casadi_axpy(nx_, 1., m->lam_xk, m->gL);
    casadi_mv(m->qpA, spA_, m->lam_g, m->gL, true);

    double time2 = clock();
    m->t_eval_mat += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::eval_res(ScpgenMemory* m) const {
    // Get current time
    double time1 = clock();

    // Inputs
    fill_n(m->arg, res_fcn_.n_in(), nullptr);
    m->arg[res_p_] = m->p; // Parameters
    m->arg[res_x_] = m->x; // Non-lifted primal variables
    for (size_t i=0; i<v_.size(); ++i) { // Lifted primal variables
      m->arg[v_[i].res_var] = m->lifted_mem[i].x;
    }
    if (!gauss_newton_) {
      m->arg[res_g_lam_] = nullptr; // Non-lifted dual variables
      for (size_t i=0; i<v_.size(); ++i) { // Lifted dual variables
        m->arg[v_[i].res_lam] = m->lifted_mem[i].lam;
      }
    }

    // Outputs
    fill_n(m->res, res_fcn_.n_out(), nullptr);
    m->res[res_f_] = &m->f; // Objective
    m->res[res_gl_] = gauss_newton_ ? m->b_gn : m->gfk; // Objective gradient
    m->res[res_g_] = m->gk; // Constraints
    for (size_t i=0; i<v_.size(); ++i) {
      m->res[v_[i].res_d] = m->lifted_mem[i].res;
      if (!gauss_newton_) {
        m->res[v_[i].res_lam_d] = m->lifted_mem[i].resL;
      }
    }
    m->res[res_p_d_] = m->lam_p; // Parameter sensitivities

    // Evaluate residual function
    res_fcn_(m->arg, m->res, m->iw, m->w, 0);

    double time2 = clock();
    m->t_eval_res += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::eval_vec(ScpgenMemory* m) const {
    // Get current time
    double time1 = clock();

    // Inputs
    fill_n(m->arg, vec_fcn_.n_in(), nullptr);
    m->arg[mod_p_] = m->p; // Parameters
    m->arg[mod_x_] = m->x; // Primal step/variables
    for (size_t i=0; i<v_.size(); ++i) {
      m->arg[v_[i].mod_var] = m->lifted_mem[i].res;
    }
    if (!gauss_newton_) {
      m->arg[mod_g_lam_] = nullptr; // Dual steps/variables
      for (size_t i=0; i<v_.size(); ++i) {
        m->arg[v_[i].mod_lam] = m->lifted_mem[i].resL;
      }
    }

    // Outputs
    fill_n(m->res, vec_fcn_.n_out(), nullptr);
    m->res[vec_gf_] = m->qpG;
    m->res[vec_g_] = m->qpB;

    // Calculate condensed QP vectors
    vec_fcn_(m->arg, m->res, m->iw, m->w, 0);

    // Linear offset in the reduced QP
    casadi_scal(ng_, -1., m->qpB);
    casadi_axpy(ng_, 1., m->gk, m->qpB);

    // Gradient of the objective in the reduced QP
    if (gauss_newton_) {
      casadi_axpy(ngn_, -1., m->qpG, m->b_gn);
    } else {
      casadi_axpy(nx_, -1., m->qpG, m->gfk);
    }

    double time2 = clock();
    m->t_eval_vec += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::regularize(ScpgenMemory* m) const {
    casadi_assert_dev(nx_==2 && spH_.is_dense());

    // Regularization
    m->reg = 0;

    // Check the smallest eigenvalue of the Hessian
    double a = m->qpH[0];
    double b = m->qpH[2];
    double c = m->qpH[1];
    double d = m->qpH[3];

    // Make sure no not a numbers
    casadi_assert_dev(a==a && b==b && c==c &&  d==d);

    // Make sure symmetric
    if (b!=c) {
      if (fabs(b-c)>=1e-10) casadi_warning("Hessian is not symmetric: "
                                           + str(b) + " != " + str(c));
      m->qpH[1] = c = b;
    }

    double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
    if (eig_smallest<reg_threshold_) {
      // Regularization
      m->reg = reg_threshold_-eig_smallest;
      m->qpH[0] += m->reg;
      m->qpH[3] += m->reg;
    }
  }

  void Scpgen::solve_qp(ScpgenMemory* m) const {
    // Get current time
    double time1 = clock();

    // Get bounds on step
    casadi_copy(m->lbx, nx_, m->qp_lbx);
    casadi_copy(m->ubx, nx_, m->qp_ubx);
    casadi_copy(m->lbg, ng_, m->qp_lba);
    casadi_copy(m->ubg, ng_, m->qp_uba);
    casadi_axpy(nx_, -1., m->x, m->qp_lbx);
    casadi_axpy(nx_, -1., m->x, m->qp_ubx);
    casadi_axpy(ng_, -1., m->qpB, m->qp_lba);
    casadi_axpy(ng_, -1., m->qpB, m->qp_uba);

    // Inputs
    fill_n(m->arg, qpsol_.n_in(), nullptr);
    m->arg[CONIC_H] = m->qpH;
    m->arg[CONIC_G] = m->gfk;
    m->arg[CONIC_A] = m->qpA;
    m->arg[CONIC_LBX] = m->qp_lbx;
    m->arg[CONIC_UBX] = m->qp_ubx;
    m->arg[CONIC_LBA] = m->qp_lba;
    m->arg[CONIC_UBA] = m->qp_uba;

    // Outputs
    fill_n(m->res, qpsol_.n_out(), nullptr);
    m->res[CONIC_X] = m->dxk; // Condensed primal step
    m->res[CONIC_LAM_X] = m->dlam_xk; // Multipliers (simple bounds)
    m->res[CONIC_LAM_A] = m->dlam_gk; // Multipliers (linear bounds)

    // Solve the QP
    qpsol_(m->arg, m->res, m->iw, m->w, 0);

    // Calculate penalty parameter of merit function
    m->sigma = merit_start_;
    m->sigma = std::max(m->sigma, 1.01*casadi_norm_inf(nx_, m->dlam_xk));
    m->sigma = std::max(m->sigma, 1.01*casadi_norm_inf(ng_, m->dlam_gk));

    // Calculate step in multipliers
    casadi_axpy(nx_, -1., m->lam_xk, m->dlam_xk);
    casadi_axpy(ng_, -1., m->lam_g, m->dlam_gk);

    double time2 = clock();
    m->t_solve_qp += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::line_search(ScpgenMemory* m, casadi_int& ls_iter, bool& ls_success) const {
    // Make sure that we have a decent direction
    if (!gauss_newton_) {
      // Get the curvature in the step direction
      double gain = casadi_bilin(m->qpH, spH_, m->dxk, m->dxk);
      if (gain < 0) {
        m->iteration_note = "Hessian indefinite in the search direction";
      }
    }

    // Calculate L1-merit function in the actual iterate
    double l1_infeas = primalInfeasibility(m);

    // Right-hand side of Armijo condition
    double F_sens = casadi_dot(nx_, m->dxk, m->gfk);
    double L1dir = F_sens - m->sigma * l1_infeas;
    double L1merit = m->f + m->sigma * l1_infeas;

    // Storing the actual merit function value in a list
    m->merit_mem[m->merit_ind] = L1merit;
    ++m->merit_ind %= merit_memsize_;

    // Calculating maximal merit function value so far
    double meritmax = m->merit_mem[0];
    for (size_t i=1; i<merit_memsize_ && i<m->iter_count; ++i) {
      if (meritmax < m->merit_mem[i]) meritmax = m->merit_mem[i];
    }

    // Stepsize
    double t = 1.0, t_prev = 0.0;

    // Merit function value in candidate
    double L1merit_cand = 0;

    // Reset line-search counter, success marker
    ls_iter = 0;
    ls_success = false;

    // Line-search
    //if (verbose_) casadi_message("Starting line-search");

    // Line-search loop
    while (true) {
      // Take the primal step
      double dt = t-t_prev;
      casadi_axpy(nx_, dt, m->dxk, m->x);
      for (auto&& v : m->lifted_mem) casadi_axpy(v.n, dt, v.dx, v.x);

      // Take the dual step
      casadi_axpy(ng_, dt, m->dlam_gk, m->lam_g);
      casadi_axpy(nx_, dt, m->dlam_xk, m->lam_xk);
      if (!gauss_newton_) {
        for (auto&& v : m->lifted_mem) casadi_axpy(v.n, dt, v.dlam, v.lam);
      }

      // Evaluate residual function to get objective and constraints
      // (and residuals for the next iteration)
      eval_res(m);
      ls_iter++;

      // Calculating merit-function in candidate
      l1_infeas = primalInfeasibility(m);
      L1merit_cand = m->f + m->sigma * l1_infeas;
      if (L1merit_cand <= meritmax + t * c1_ * L1dir) {

        // Accepting candidate
        ls_success = true;
        //if (verbose_) casadi_message("Line-search completed, candidate accepted");
        break;
      }

      // Line-search not successful, but we accept it.
      if (ls_iter == max_iter_ls_) {
        //if (verbose_) casadi_message("Line-search completed, maximum number of iterations");
        break;
      }

      // Backtracking
      t_prev = t;
      t = beta_ * t;
    }

    // Calculate primal step-size
    m->pr_step = casadi_norm_1(nx_, m->dxk);
    for (auto&& v : m->lifted_mem) m->pr_step += casadi_norm_1(v.n, v.dx);
    m->pr_step *= t;

    // Calculate the dual step-size
    m->du_step = casadi_norm_1(ng_, m->dlam_gk) + casadi_norm_1(nx_, m->dlam_xk);
    for (auto&& v : m->lifted_mem) m->du_step += casadi_norm_1(v.n, v.dlam);
    m->du_step *= t;
  }

  void Scpgen::eval_exp(ScpgenMemory* m) const {
    // Get current time
    double time1 = clock();

    // Inputs
    fill_n(m->arg, exp_fcn_.n_in(), nullptr);
    m->arg[mod_p_] = m->p; // Parameter
    m->arg[mod_du_] = m->dxk; // Primal step
    m->arg[mod_x_] = m->x; // Primal variables
    for (size_t i=0; i<v_.size(); ++i) {
      m->arg[v_[i].mod_var] = m->lifted_mem[i].res;
    }
    if (!gauss_newton_) {
      m->arg[mod_dlam_g_] = m->dlam_gk; // Dual variables
      m->arg[mod_g_lam_] = m->lam_g; // Dual step
      for (size_t i=0; i<v_.size(); ++i) {
        m->arg[v_[i].mod_lam] = m->lifted_mem[i].resL;
      }
    }

    // Outputs
    fill_n(m->res, exp_fcn_.n_out(), nullptr);
    for (casadi_int i=0; i<v_.size(); ++ i) {
      m->res[v_[i].exp_def] = m->lifted_mem[i].dx; // Expanded primal step
      if (!gauss_newton_) {
        m->res[v_[i].exp_defL] = m->lifted_mem[i].dlam; // Expanded dual step
      }
    }

    // Perform the step expansion
    exp_fcn_(m->arg, m->res, m->iw, m->w, 0);

    double time2 = clock();
    m->t_eval_exp += (time2-time1)/CLOCKS_PER_SEC;
  }

  Dict Scpgen::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<ScpgenMemory*>(mem);
    stats["t_eval_mat"] = m->t_eval_mat;
    stats["t_eval_res"] = m->t_eval_res;
    stats["t_eval_vec"] = m->t_eval_vec;
    stats["t_eval_exp"] = m->t_eval_exp;
    stats["t_solve_qp"] = m->t_solve_qp;
    stats["t_mainloop"] = m->t_mainloop;
    stats["iter_count"] = m->iter_count;
    return stats;
  }



} // namespace casadi
