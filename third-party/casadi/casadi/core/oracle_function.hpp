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


#ifndef CASADI_ORACLE_FUNCTION_HPP
#define CASADI_ORACLE_FUNCTION_HPP

#include "function_internal.hpp"
#include "timing.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Function memory with temporary work vectors */
  struct CASADI_EXPORT OracleMemory {
    // Work vectors
    const double** arg;
    double** res;
    casadi_int* iw;
    double* w;

    // Function specific statistics
    std::map<std::string, FStats> fstats;

    // Add a statistic
    void add_stat(const std::string& s) {
      bool added = fstats.insert(std::make_pair(s, FStats())).second;
      casadi_assert(added, "Duplicate stat: '" + s + "'");
    }
  };

  /** \brief Base class for functions that perform calculation with an oracle
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT OracleFunction : public FunctionInternal {
  protected:
    /// Oracle: Used to generate other functions
    Function oracle_;

    /// Options for creating functions
    Dict common_options_;
    Dict specific_options_;

    // Information about one function
    struct RegFun {
      Function f;
      bool jit;
      bool monitored = false;
    };

    // All NLP functions
    std::map<std::string, RegFun> all_functions_;
  public:
    /** \brief  Constructor */
    OracleFunction(const std::string& name, const Function& oracle);

    /** \brief  Destructor */
    ~OracleFunction() override = 0;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** Initialize  */
    void init(const Dict& opts) override;

    /// Finalize initialization
    void finalize(const Dict& opts) override;

    /** \brief Get oracle */
    const Function& oracle() const override { return oracle_;}

    // Replace MX oracle with SX oracle?
    void expand();

    /** Create an oracle function */
    Function
    create_function(const std::string& fname,
                    const std::vector<std::string>& s_in,
                    const std::vector<std::string>& s_out,
                    const Function::AuxOut& aux=Function::AuxOut());

    /** Register the function for evaluation and statistics gathering */
    void set_function(const Function& fcn, const std::string& fname, bool jit=false);

    /** Register the function for evaluation and statistics gathering */
    void set_function(const Function& fcn) { set_function(fcn, fcn.name()); }

    // Calculate an oracle function
    casadi_int calc_function(OracleMemory* m, const std::string& fcn,
                      const double* const* arg=nullptr) const;

    // Get list of dependency functions
    std::vector<std::string> get_function() const override;

    // Get a dependency function
    const Function& get_function(const std::string &name) const override;

    // Is a function monitored?
    virtual bool monitored(const std::string &name) const;

    // Check if a particular dependency exists
    bool has_function(const std::string& fname) const override;

    /** \brief Export / Generate C code for the generated functions */
    std::string generate_dependencies(const std::string& fname, const Dict& opts) const override;

    /** \brief JIT for dependencies */
    void jit_dependencies(const std::string& fname) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new OracleMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<OracleMemory*>(mem);}

    /** \brief Set the work vectors */
    void set_temp(void* mem, const double** arg, double** res,
                          casadi_int* iw, double* w) const override;

    /// Print statistics
    void print_fstats(const OracleMemory* m) const;

    /// Get all statistics
    Dict get_stats(void* mem) const override;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_ORACLE_FUNCTION_HPP
