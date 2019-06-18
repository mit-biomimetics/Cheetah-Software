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


#ifndef CASADI_ROOTFINDER_IMPL_HPP
#define CASADI_ROOTFINDER_IMPL_HPP

#include "rootfinder.hpp"
#include "oracle_function.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT RootfinderMemory : public OracleMemory {
    // Inputs
    const double** iarg;

    // Outputs
    double** ires;

    // Success?
    bool success;
  };

  /// Internal class
  class CASADI_EXPORT
  Rootfinder : public OracleFunction, public PluginInterface<Rootfinder> {
  public:
    /** \brief Constructor
     *
     * \param f   Function mapping from (n+1) inputs to 1 output.
     */
    Rootfinder(const std::string& name, const Function& oracle);

    /// Destructor
    ~Rootfinder() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return oracle_.n_in();}
    size_t get_n_out() override { return oracle_.n_out();}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override { return oracle_.sparsity_in(i);}
    Sparsity get_sparsity_out(casadi_int i) override { return oracle_.sparsity_out(i);}
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override { return oracle_.name_in(i);}
    std::string get_name_out(casadi_int i) override { return oracle_.name_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    // Solve the NLP
    virtual int solve(void* mem) const = 0;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res,
                    casadi_int* iw, bvec_t* w, void* mem) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

    ///@{
    /// Is the class able to propagate seeds through the algorithm?
    bool has_spfwd() const override { return true;}
    bool has_sprev() const override { return true;}
    ///@}

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    bool uses_output() const override {return true;}

    ///@{
    /** \brief Generate a function that calculates \a nfwd forward derivatives */
    bool has_forward(casadi_int nfwd) const override { return true;}
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Generate a function that calculates \a nadj adjoint derivatives */
    bool has_reverse(casadi_int nadj) const override { return true;}
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /** \brief Create call to (cached) derivative function, forward mode  */
    virtual void ad_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens,
                         bool always_inline, bool never_inline) const;

    /** \brief Create call to (cached) derivative function, reverse mode  */
    virtual void ad_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens,
                         bool always_inline, bool never_inline) const;

    /// Number of equations
    casadi_int n_;

    /// Linear solver
    Linsol linsol_;
    Sparsity sp_jac_;

    /// Constraints on decision variables
    std::vector<casadi_int> u_c_;

    /// Indices of the input and output that correspond to the actual root-finding
    casadi_int iin_, iout_;

    /// Throw an exception on failure?
    bool error_on_fail_;

    // Creator function for internal class
    typedef Rootfinder* (*Creator)(const std::string& name, const Function& oracle);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Short name
    static std::string shortname() { return "rootfinder";}

    /// Infix
    static const std::string infix_;

    /// Convert dictionary to Problem
    template<typename XType>
      static Function create_oracle(const std::map<std::string, XType>& d,
                                    const Dict& opts);
  };



} // namespace casadi
/// \endcond

#endif // CASADI_ROOTFINDER_IMPL_HPP
