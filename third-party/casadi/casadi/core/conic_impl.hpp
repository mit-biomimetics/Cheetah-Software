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


#ifndef CASADI_CONIC_IMPL_HPP
#define CASADI_CONIC_IMPL_HPP

#include "conic.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"
#include "timing.hpp"

/// \cond INTERNAL
namespace casadi {

  struct CASADI_EXPORT ConicMemory {
    // Function specific statistics
    std::map<std::string, FStats> fstats;
  };

  /// Internal class
  class CASADI_EXPORT Conic : public FunctionInternal, public PluginInterface<Conic> {
  public:

    // Constructor
    Conic(const std::string& name, const std::map<std::string, Sparsity> &st);

    // Destructor
    ~Conic() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return CONIC_NUM_IN;}
    size_t get_n_out() override { return CONIC_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override { return conic_in(i);}
    std::string get_name_out(casadi_int i) override { return conic_out(i);}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize
    void init(const Dict& opts) override;

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void check_inputs(const double* lbx, const double* ubx,
                             const double* lba, const double* uba) const;

    /** Generate native code in the interfaced language for debugging */
    virtual void generateNativeCode(std::ostream& file) const;

    // Creator function for internal class
    typedef Conic* (*Creator)(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "conic";}

    /** \brief Check if the function is of a particular type */
    bool is_a(const std::string& type, bool recursive) const override;

    /** \brief Get default input value */
    double get_default_in(casadi_int ind) const override;

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /// Can psd constraints be treated
    virtual bool psd_support() const { return false;}

    /// Print statistics
    void print_fstats(const ConicMemory* m) const;

  protected:
    /// Options
    std::vector<bool> discrete_;

    /// Problem structure
    Sparsity H_, A_, Q_, P_;

    /// Number of decision variables
    casadi_int nx_;

    /// The number of constraints (counting both equality and inequality) == A.size1()
    casadi_int na_;

    /// The shape of psd constraint matrix
    casadi_int np_;
  };


} // namespace casadi
/// \endcond
#endif // CASADI_CONIC_IMPL_HPP
