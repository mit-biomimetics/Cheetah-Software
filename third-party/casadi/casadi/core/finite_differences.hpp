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


#ifndef CASADI_FINITE_DIFFERENCES_HPP
#define CASADI_FINITE_DIFFERENCES_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {
  /** Calculate derivative using finite differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT FiniteDiff : public FunctionInternal {
  public:
    // Constructor (protected, use create function)
    FiniteDiff(const std::string& name, casadi_int n);

    /** \brief Destructor */
    ~FiniteDiff() override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    /** \brief Get default input value */
    double get_default_in(casadi_int ind) const override;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    // Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief Is the scheme using the (nondifferentiated) output? */
    bool uses_output() const override {return true;}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

  protected:
    // Number of function evaluations needed
    virtual casadi_int n_pert() const = 0;

    // Get perturbation expression
    virtual std::string pert(const std::string& k) const = 0;

    // Get perturbation expression
    virtual double pert(casadi_int k, double h) const = 0;

    // Calculate finite difference approximation
    virtual double calc_fd(double** yk, double* y0, double* J, double h) const = 0;

    // Codegen finite difference approximation
    virtual std::string calc_fd() const = 0;

    // Is an error estimate available?
    virtual casadi_int has_err() const = 0;

    // Calculate step size from absolute tolerance
    virtual double calc_stepsize(double abstol) const = 0;

    // Number of directional derivatives
    casadi_int n_;

    // Iterations to improve h
    casadi_int h_iter_;

    // Perturbation
    double h_;

    // Dimensions
    casadi_int n_z_, n_y_;

    // Target ratio of truncation error to roundoff error
    double u_aim_;

    // Allowed step size range
    double h_min_, h_max_;

    // Memory object
    casadi_finite_diff_mem<double> m_;
  };

  /** Calculate derivative using forward differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT ForwardDiff : public FiniteDiff {
  public:
    // Constructor
    ForwardDiff(const std::string& name, casadi_int n) : FiniteDiff(name, n) {}

    /** \brief Destructor */
    ~ForwardDiff() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "ForwardDiff";}

    // Number of function evaluations needed
    casadi_int n_pert() const override {return 1;};

    // Get perturbation expression
    std::string pert(const std::string& k) const override {
      return str(h_);
    }

    // Get perturbation expression
    double pert(casadi_int k, double h) const override {
      return h;
    }

    // Calculate finite difference approximation
    double calc_fd(double** yk, double* y0, double* J, double h) const override;

    // Codegen finite difference approximation
    std::string calc_fd() const override {return "casadi_forward_diff";}

    // Is an error estimate available?
    casadi_int has_err() const override {return false;}

    // Calculate step size from absolute tolerance
    double calc_stepsize(double abstol) const override { return sqrt(abstol);}

    /** \brief Get absolute tolerance */
    double get_abstol() const override { return h_;}

    ///@{
    /** \brief Second order derivatives */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}
  };

  /** Calculate derivative using backward differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT BackwardDiff : public ForwardDiff {
  public:
    // Constructor
    BackwardDiff(const std::string& name, casadi_int n) : ForwardDiff(name, n) {}

    /** \brief Destructor */
    ~BackwardDiff() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "BackwardDiff";}

    // Calculate step size from absolute tolerance
    double calc_stepsize(double abstol) const override {
      return -ForwardDiff::calc_stepsize(abstol);
    }
  };

  /** Calculate derivative using central differences
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT CentralDiff : public FiniteDiff {
  public:
    // Constructor
    CentralDiff(const std::string& name, casadi_int n) : FiniteDiff(name, n) {}

    /** \brief Destructor */
    ~CentralDiff() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "CentralDiff";}

    // Number of function evaluations needed
    casadi_int n_pert() const override {return 2;};

    // Get perturbation expression
    std::string pert(const std::string& k) const override {
      return "(2*" + k + "-1)*" + str(h_);
    }

    // Get perturbation expression
    double pert(casadi_int k, double h) const override {
      return (2*static_cast<double>(k)-1)*h;
    }

    // Calculate finite difference approximation
    double calc_fd(double** yk, double* y0, double* J, double h) const override;

    // Codegen finite difference approximation
    std::string calc_fd() const override {return "casadi_central_diff";}

    // Is an error estimate available?
    casadi_int has_err() const override {return true;}

    // Calculate step size from absolute tolerance
    double calc_stepsize(double abstol) const override { return pow(abstol, 1./3);}

    /** \brief Get absolute tolerance */
    double get_abstol() const override { return h_*h_;}

    ///@{
    /** \brief Second order derivatives */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}
  };

  /** Calculate derivative using 3th order smoothing scheme
    * \author Joel Andersson
    * \date 2017
  */
  class CASADI_EXPORT Smoothing : public FiniteDiff {
  public:
    // Constructor
    Smoothing(const std::string& name, casadi_int n) : FiniteDiff(name, n) {}

    /** \brief Destructor */
    ~Smoothing() override {}

    /** \brief Get type name */
    std::string class_name() const override {return "Smoothing";}

    // Number of function evaluations needed
    casadi_int n_pert() const override {return 4;};

    // Get perturbation expression
    std::string pert(const std::string& k) const override;

    // Get perturbation expression
    double pert(casadi_int k, double h) const override;

    // Calculate finite difference approximation
    double calc_fd(double** yk, double* y0, double* J, double h) const override;

    // Codegen finite difference approximation
    std::string calc_fd() const override {return "casadi_smoothing_diff";}

    // Is an error estimate available?
    casadi_int has_err() const override {return true;}

    // Calculate step size from absolute tolerance
    double calc_stepsize(double abstol) const override { return pow(abstol, 1./3);}

    /** \brief Get absolute tolerance */
    double get_abstol() const override { return h_*h_;}

    ///@{
    /** \brief Second order derivatives */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}
  };


} // namespace casadi
/// \endcond

#endif // CASADI_FINITE_DIFFERENCES_HPP
