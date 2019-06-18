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


#include "finite_differences.hpp"

using namespace std;

namespace casadi {

  FiniteDiff::FiniteDiff(const std::string& name, casadi_int n)
    : FunctionInternal(name), n_(n) {
  }

  FiniteDiff::~FiniteDiff() {
  }

  Options FiniteDiff::options_
  = {{&FunctionInternal::options_},
     {{"second_order_stepsize",
       {OT_DOUBLE,
        "Second order perturbation size [default: 1e-3]"}},
      {"h_max",
       {OT_DOUBLE,
        "Maximum step size [default 0]"}},
      {"h_min",
       {OT_DOUBLE,
        "Minimum step size [default inf]"}},
      {"smoothing",
       {OT_DOUBLE,
        "Smoothing regularization [default: machine precision]"}},
      {"reltol",
       {OT_DOUBLE,
        "Accuracy of function inputs [default: query object]"}},
      {"abstol",
        {OT_DOUBLE,
        "Accuracy of function outputs [default: query object]"}},
      {"u_aim",
        {OT_DOUBLE,
        "Target ratio of roundoff error to truncation error [default: 100.]"}},
      {"h_iter",
        {OT_INT,
        "Number of iterations to improve on the step-size "
        "[default: 1 if error estimate available, otherwise 0]"}},
     }
  };

  void FiniteDiff::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Default options
    h_min_ = 0;
    h_max_ = inf;
    m_.smoothing = eps;
    m_.reltol = derivative_of_->get_reltol();
    m_.abstol = derivative_of_->get_abstol();
    h_ = calc_stepsize(m_.abstol);
    u_aim_ = 100;
    h_iter_ = has_err() ? 1 : 0;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="h") {
        h_ = op.second;
      } else if (op.first=="h_min") {
        h_min_ = op.second;
      } else if (op.first=="h_max") {
        h_max_ = op.second;
      } else if (op.first=="reltol") {
        m_.reltol = op.second;
      } else if (op.first=="abstol") {
        m_.abstol = op.second;
      } else if (op.first=="smoothing") {
        m_.smoothing = op.second;
      } else if (op.first=="u_aim") {
        u_aim_ = op.second;
      } else if (op.first=="h_iter") {
        h_iter_ = op.second;
      }
    }

    // Check h_iter for consistency
    if (h_iter_!=0 && !has_err()) {
      casadi_error("Perturbation size refinement requires an error estimate, "
      "which is not available for the class '" + class_name() + "'. "
      "Choose a different differencing scheme.");
    }

    // Allocate work vector for (perturbed) inputs and outputs
    n_z_ = derivative_of_.nnz_in();
    n_y_ = derivative_of_.nnz_out();
    alloc_res(n_pert(), true); // yk
    alloc_w((n_pert() + 3) * n_y_, true); // yk[:], y0, y, J
    alloc_w(n_z_, true); // z

    // Dimensions
    if (verbose_) {
      casadi_message("Finite differences (" + class_name() + ") with "
                     + str(n_z_) + " inputs, " + str(n_y_)
                     + " outputs and " + str(n_) + " directional derivatives.");
    }

    // Allocate sufficient temporary memory for function evaluation
    alloc(derivative_of_);
  }

  Sparsity FiniteDiff::get_sparsity_in(casadi_int i) {
    casadi_int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    if (i<n_in) {
      // Non-differentiated input
      return derivative_of_.sparsity_in(i);
    } else if (i<n_in+n_out) {
      // Non-differentiated output
      return derivative_of_.sparsity_out(i-n_in);
    } else {
      // Seeds
      return repmat(derivative_of_.sparsity_in(i-n_in-n_out), 1, n_);
    }
  }

  Sparsity FiniteDiff::get_sparsity_out(casadi_int i) {
    return repmat(derivative_of_.sparsity_out(i), 1, n_);
  }

  double FiniteDiff::get_default_in(casadi_int ind) const {
    if (ind<derivative_of_.n_in()) {
      return derivative_of_.default_in(ind);
    } else {
      return 0;
    }
  }

  size_t FiniteDiff::get_n_in() {
    return derivative_of_.n_in() + derivative_of_.n_out() + derivative_of_.n_in();
  }

  size_t FiniteDiff::get_n_out() {
    return derivative_of_.n_out();
  }

  std::string FiniteDiff::get_name_in(casadi_int i) {
    casadi_int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    if (i<n_in) {
      return derivative_of_.name_in(i);
    } else if (i<n_in+n_out) {
      return "out_" + derivative_of_.name_out(i-n_in);
    } else {
      return "fwd_" + derivative_of_.name_in(i-n_in-n_out);
    }
  }

  std::string FiniteDiff::get_name_out(casadi_int i) {
    return "fwd_" + derivative_of_.name_out(i);
  }

  Function CentralDiff::get_forward(casadi_int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    // Commented out, does not work well
#if 0
    // The second order derivative is calculated as the backwards derivative
    // of the forward derivative, which is equivalent to central differences
    // of second order
    string f_name = "fd_" + name;
    Dict f_opts = {{"derivative_of", derivative_of_}};
    Function f = Function::create(new ForwardDiff(f_name, n_, h_), f_opts);
    // Calculate backwards derivative of f
    f_opts["derivative_of"] = f;
    return Function::create(new ForwardDiff(name, nfwd, -h_), f_opts);
#endif
    return Function::create(new CentralDiff(name, nfwd), opts);
  }

  int FiniteDiff::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    // Shorthands
    casadi_int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    casadi_int n_pert = this->n_pert();

    // Non-differentiated input
    const double** x0 = arg;
    arg += n_in;

    // Non-differentiated output
    double* y0 = w;
    for (casadi_int j=0; j<n_out; ++j) {
      const casadi_int nnz = derivative_of_.nnz_out(j);
      casadi_copy(*arg++, nnz, w);
      w += nnz;
    }

    // Forward seeds
    const double** seed = arg;
    arg += n_in;

    // Forward sensitivities
    double** sens = res;
    res += n_out;

    // Finite difference approximation
    double* J = w;
    w += n_y_;

    // Perturbed function values
    double** yk = res;
    res += n_pert;
    for (casadi_int j=0; j<n_pert; ++j) {
      yk[j] = w, w += n_y_;
    }

    // Setup arg and z for evaluation
    double *z = w;
    for (casadi_int j=0; j<n_in; ++j) {
      arg[j] = w;
      w += derivative_of_.nnz_in(j);
    }

    // Setup res and y for evaluation
    double *y = w;
    for (casadi_int j=0; j<n_out; ++j) {
      res[j] = w;
      w += derivative_of_.nnz_out(j);
    }

    // For all sensitivity directions
    for (casadi_int i=0; i<n_; ++i) {
      // Initial stepsize
      double h = h_;
      // Perform finite difference algorithm with different step sizes
      for (casadi_int iter=0; iter<1+h_iter_; ++iter) {
        // Calculate perturbed function values
        for (casadi_int k=0; k<n_pert; ++k) {
          // Perturb inputs
          casadi_int off = 0;
          for (casadi_int j=0; j<n_in; ++j) {
            casadi_int nnz = derivative_of_.nnz_in(j);
            casadi_copy(x0[j], nnz, z + off);
            //cout << "k = " << k << ": pert(k, h) = " << pert(k, h) << endl;
            if (seed[j]) casadi_axpy(nnz, pert(k, h), seed[j] + i*nnz, z + off);
            off += nnz;
          }
          // Evaluate
          if (derivative_of_(arg, res, iw, w)) return 1;
          // Save outputs
          casadi_copy(y, n_y_, yk[k]);
        }
        // Finite difference calculation with error estimate
        double u = calc_fd(yk, y0, J, h);
        if (iter==h_iter_) break;

        // Update step size
        if (u < 0) {
          // Perturbation failed, try a smaller step size
          h /= u_aim_;
        } else {
          // Update h to get u near the target ratio
          h *= sqrt(u_aim_ / fmax(1., u));
        }
        // Make sure h stays in the range [h_min_,h_max_]
        h = fmin(fmax(h, h_min_), h_max_);
      }

      // Gather sensitivities
      casadi_int off = 0;
      for (casadi_int j=0; j<n_out; ++j) {
        casadi_int nnz = derivative_of_.nnz_out(j);
        if (sens[j]) casadi_copy(J + off, nnz, sens[j] + i*nnz);
        off += nnz;
      }
    }
    return 0;
  }

  double ForwardDiff::calc_fd(double** yk, double* y0, double* J, double h) const {
    return casadi_forward_diff(yk, y0, J, h, n_y_, &m_);
  }

  double CentralDiff::calc_fd(double** yk, double* y0, double* J, double h) const {
    return casadi_central_diff(yk, y0, J, h, n_y_, &m_);
  }

  void FiniteDiff::codegen_declarations(CodeGenerator& g) const {
    g.add_dependency(derivative_of_);
    g.add_auxiliary(CodeGenerator::AUX_FINITE_DIFF);
  }

  void FiniteDiff::codegen_body(CodeGenerator& g) const {
    // Shorthands
    casadi_int n_in = derivative_of_.n_in(), n_out = derivative_of_.n_out();
    casadi_int n_pert = this->n_pert();

    g.comment("Non-differentiated input");
    g.local("x0", "const casadi_real", "**");
    g << "x0 = arg, arg += " << n_in << ";\n";

    g.comment("Non-differentiated output");
    g.local("y0", "casadi_real", "*");
    g << "y0 = w;\n";
    for (casadi_int j=0; j<n_out; ++j) {
      const casadi_int nnz = derivative_of_.nnz_out(j);
      g << g.copy("*arg++", nnz, "w") << " w += " << nnz << ";\n";
    }

    g.comment("Forward seeds");
    g.local("seed", "const casadi_real", "**");
    g << "seed = arg, arg += " << n_in << ";\n";

    g.comment("Forward sensitivities");
    g.local("sens", "casadi_real", "**");
    g << "sens = res, res += " << n_out << ";\n";

    g.comment("Finite difference approximation");
    g.local("J", "casadi_real", "*");
    g << "J = w, w += " << n_y_ << ";\n";

    g.comment("Perturbed function value");
    g.local("yk", "casadi_real", "**");
    g << "yk = res, res += " << n_pert << ";\n";
    g.local("j", "casadi_int");
    g << "for (j=0; j<" << n_pert << "; ++j) yk[j] = w, w += " << n_y_ << ";\n";

    g.comment("Setup arg and z for evaluation");
    g.local("z", "casadi_real", "*");
    g << "z = w;\n";
    for (casadi_int j=0; j<n_in; ++j) {
      g << "arg[" << j << "] = w, w += " << derivative_of_.nnz_in(j) << ";\n";
    }

    g.comment("Setup res and y for evaluation");
    g.local("y", "casadi_real", "*");
    g << "y = w;\n";
    for (casadi_int j=0; j<n_out; ++j) {
      g << "res[" << j << "] = w, w += " << derivative_of_.nnz_out(j) << ";\n";
    }

    g.comment("For all sensitivity directions");
    g.local("i", "casadi_int");
    g << "for (i=0; i<" << n_ << "; ++i) {\n";

    g.comment("Initial stepsize");
    g.local("h", "casadi_real");
    g << "h = " << h_ << ";\n";

    g.comment("Perform finite difference algorithm with different step sizes");
    g.local("iter", "casadi_int");
    g << "for (iter=0; iter<" << 1+h_iter_ << "; ++iter) {\n";

    g.comment("Calculate perturbed function values");
    g.local("k", "casadi_int");
    g << "for (k=0; k<" << n_pert << "; ++k) {\n";

    g.comment("Perturb inputs");
    casadi_int off=0;
    for (casadi_int j=0; j<n_in; ++j) {
      casadi_int nnz = derivative_of_.nnz_in(j);
      string s = "seed[" + str(j) + "]";
      g << g.copy("x0[" + str(j) + "]", nnz, "z+" + str(off)) << "\n"
        << "if ("+s+") " << g.axpy(nnz, pert("k"),
                                   s+"+i*"+str(nnz), "z+" + str(off)) << "\n";
      off += nnz;
    }

    g.comment("Evaluate");
    g << "if (" << g(derivative_of_, "arg", "res", "iw", "w") << ") return 1;\n";

    g.comment("Save outputs");
    g << g.copy("y", n_y_, "yk[k]") << "\n";

    g << "}\n"; // for (k=0, ...)

    g.comment("Finite difference calculation with error estimate");
    g.local("u", "casadi_real");
    g.local("m", "const struct casadi_finite_diff_mem");
    g.init_local("m", "{" + str(m_.reltol) + ", "
                          + str(m_.abstol) + ", "
                          + str(m_.smoothing) + "}");
    g << "u = " << calc_fd() << "(yk, y0, J, h, " << n_y_ << ", &m);\n";
    g << "if (iter==" << h_iter_ << ") break;\n";

    g.comment("Update step size");
    g << "if (u < 0) {\n";
    // Perturbation failed, try a smaller step size
    g << "h /= " << u_aim_ << ";\n";
    g << "} else {\n";
    // Update h to get u near the target ratio
    g << "h *= sqrt(" << u_aim_ << " / fmax(1., u));\n";
    g << "}\n";
    // Make sure h stays in the range [h_min_,h_max_]
    if (h_min_>0 || isfinite(h_max_)) {
      string h = "h";
      if (h_min_>0) h = "fmax(" + h + ", " + str(h_min_) + ")";
      if (isfinite(h_max_)) h = "fmin(" + h + ", " + str(h_max_) + ")";
      g << "h = " << h << ";\n";
    }

    g << "}\n"; // for (iter=0, ...)

    g.comment("Gather sensitivities");
    off = 0;
    for (casadi_int j=0; j<n_out; ++j) {
      casadi_int nnz = derivative_of_.nnz_out(j);
      string s = "sens[" + str(j) + "]";
      g << "if (" << s << ") " << g.copy("J+" + str(off), nnz, s + "+i*" + str(nnz)) << "\n";
      off += nnz;
    }
    g << "}\n"; // for (i=0, ...)
  }

  std::string Smoothing::pert(const std::string& k) const {
    string sign = "(2*(" + k + "/2)-1)";
    string len = "(" + k + "%%2+1)";
    return len + "*" + sign + "*" + str(h_);
  }

  double Smoothing::pert(casadi_int k, double h) const {
    casadi_int sign = 2*(k/2)-1;
    casadi_int len = k%2+1;
    return static_cast<double>(len*sign)*h;
  }

  double Smoothing::calc_fd(double** yk, double* y0, double* J, double h) const {
    return casadi_smoothing_diff(yk, y0, J, h, n_y_, &m_);
  }

  Function ForwardDiff::get_forward(casadi_int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    return Function::create(new ForwardDiff(name, nfwd), opts);
  }

  Function Smoothing::get_forward(casadi_int nfwd, const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    return Function::create(new Smoothing(name, nfwd), opts);
  }

} // namespace casadi
