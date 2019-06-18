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

#ifndef CASADI_CPLEX_INTERFACE_HPP
#define CASADI_CPLEX_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/cplex/casadi_conic_cplex_export.h>
#include "ilcplex/cplexx.h"

#include <string>

/** \defgroup plugin_Conic_cplex

      Interface to Cplex solver for sparse Quadratic Programs
*/

/** \pluginsection{Conic,cplex} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_CPLEX_EXPORT CplexMemory : public ConicMemory {
    /// Indicates if we have to warm-start
    bool is_warm;

    /// Nature of problem (always minimization)
    casadi_int objsen;

    /// Determines relation >,<, = in the linear constraints
    std::vector<char> sense;

    /// Coefficients of matrix A (constraint Jacobian)
    std::vector<int> matcnt;

    /// Right-hand side of constraints
    std::vector<double> rhs;

    /// Range of constraints
    std::vector<double> rngval;

    /// Coefficients of matrix H (objective Hessian)
    std::vector<int> qmatcnt;

    /// Storage for basis info of primal variables
    std::vector<int> cstat;

    /// Storage for basis info of slack variables
    std::vector<int> rstat;

    /// CPLEX environment
    CPXENVptr env;
    CPXLPptr lp;

    std::vector<CPXDIM> a_row, h_row;
    std::vector<CPXNNZ> a_colind, h_colind;


    int return_status;
    bool success;

    /// Constructor
    CplexMemory();

    /// Destructor
    ~CplexMemory();
  };

  /** \brief \pluginbrief{Conic,cplex}

      @copydoc Conic_doc
      @copydoc plugin_Conic_cplex

      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_CONIC_CPLEX_EXPORT CplexInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new CplexInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit CplexInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~CplexInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "cplex";}

    // Get name of the class
    std::string class_name() const override { return "CplexInterface";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new CplexMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<CplexMemory*>(mem);}

    // Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// Can discrete variables be treated
    bool integer_support() const override { return true;}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /// All CPLEX options
    Dict opts_;

    ///@{
    /// Options
    casadi_int qp_method_;
    bool dump_to_file_;
    std::string dump_filename_;
    double tol_;
    casadi_int dep_check_;
    bool warm_start_;
    ///@}

    // Are we solving a mixed-integer problem?
    bool mip_;
    std::vector<char> ctype_;

    /// A documentation string
    static const std::string meta_doc;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CPLEX_INTERFACE_HPP
