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


#ifndef CASADI_NLPSOL_IMPL_HPP
#define CASADI_NLPSOL_IMPL_HPP

#include "nlpsol.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT NlpsolMemory : public OracleMemory {
    // Bounds, given parameter values
    const double *lbx, *ubx, *lbg, *ubg, *p;

    // Current primal solution
    double *x;

    // Current dual solution
    double *lam_g, *lam_x, *lam_p;

    // Outputs
    double f, *g;

    // number of iterations
    casadi_int n_iter;

    // Success?
    bool success;
  };

  /** \brief NLP solver storage class

      @copydoc Nlpsol_doc
      \author Joel Andersson
      \date 2010-2013
  */
  class CASADI_EXPORT
  Nlpsol : public OracleFunction, public PluginInterface<Nlpsol> {
  public:
    /// Number of variables
    casadi_int nx_;

    /// Number of constraints
    casadi_int ng_;

    /// Number of parameters
    casadi_int np_;

    /// callback function, executed at each iteration
    Function fcallback_;

    /// Execute the callback function only after this amount of iterations
    casadi_int callback_step_;

    /// Throw an exception on failure?
    bool error_on_fail_;

    ///@{
    /** \brief Options */
    bool eval_errors_fatal_;
    bool warn_initial_bounds_;
    bool iteration_callback_ignore_errors_;
    bool calc_multipliers_;
    bool calc_lam_x_, calc_lam_p_, calc_f_, calc_g_;
    bool bound_consistency_;
    bool no_nlp_grad_;
    std::vector<bool> discrete_;
    ///@}

    // Mixed integer problem?
    bool mi_;

    /// Cache for KKT function
    mutable WeakRef kkt_;

    /// Constructor
    Nlpsol(const std::string& name, const Function& oracle);

    /// Destructor
    ~Nlpsol() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return NLPSOL_NUM_IN;}
    size_t get_n_out() override { return NLPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override { return nlpsol_in(i);}
    std::string get_name_out(casadi_int i) override { return nlpsol_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Print description */
    void disp_more(std::ostream& stream) const override;

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new NlpsolMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<NlpsolMemory*>(mem);}

    /** \brief Check if the inputs correspond to a well-posed problem */
    virtual void check_inputs(void* mem) const;

    /** \brief Get default input value */
    double get_default_in(casadi_int ind) const override { return nlpsol_default_in(ind);}

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    // Solve the NLP
    virtual int solve(void* mem) const = 0;

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    bool uses_output() const override {return true;}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    ///@{
    /** \brief Generate a function that calculates forward mode derivatives */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Generate a function that calculates reverse mode derivatives */
    bool has_reverse(casadi_int nadj) const override { return true;}
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    // Call the callback function
    int callback(void* mem, const double* x, const double* f, const double* g,
                 const double* lam_x, const double* lam_g, const double* lam_p) const;

    // Get KKT function
    Function kkt() const;

    // Make sure primal-dual solution is consistent with bounds
    static void bound_consistency(casadi_int n, double* x, double* lam,
                                  const double* lbx, const double* ubx);

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const Function& oracle);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "nlpsol";}

    // Get reduced Hessian
    virtual DM getReducedHessian();

    /// Read options from parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// WORKAROUND: Add an element to an std::vector stored in a GenericType:
    template<typename Type> static void append_to_vec(GenericType& t, Type el) {
      std::vector<Type> v = t;
      v.push_back(el);
      t = v;
    }

    /// Convert dictionary to Problem
    template<typename XType>
      static Function create_oracle(const std::map<std::string, XType>& d,
                                    const Dict& opts);
  };

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_IMPL_HPP
