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


#ifndef CASADI_GUROBI_INTERFACE_HPP
#define CASADI_GUROBI_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/gurobi/casadi_conic_gurobi_export.h>

// GUROBI header
extern "C" {
#include "gurobi_c.h" // NOLINT(build/include)
}

/** \defgroup plugin_Conic_gurobi
    Interface to the GUROBI Solver for quadratic programming
*/

/** \pluginsection{Conic,gurobi} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_CONIC_GUROBI_EXPORT GurobiMemory : public ConicMemory {
    // Gurobi environment
    GRBenv *env;

    int return_status;
    bool success;

    /// Constructor
    GurobiMemory();

    /// Destructor
    ~GurobiMemory();
  };

  /** \brief \pluginbrief{Conic,gurobi}

      @copydoc Conic_doc
      @copydoc plugin_Conic_gurobi

  */
  class CASADI_CONIC_GUROBI_EXPORT GurobiInterface : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit GurobiInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new GurobiInterface(name, st);
    }

    /** \brief  Destructor */
    ~GurobiInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "gurobi";}

    // Get name of the class
    std::string class_name() const override { return "GurobiInterface";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new GurobiMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<GurobiMemory*>(mem);}

    /// Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// Can discrete variables be treated
    bool integer_support() const override { return true;}

    /// Can psd constraints be treated
    bool psd_support() const override { return true;}

    /// A documentation string
    static const std::string meta_doc;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Variable types
    std::vector<char> vtype_;

    /// Gurobi options
    Dict opts_;

    // Block partition vector for SOCP (block i runs from r_[i] to r_[i+1])
    std::vector<casadi_int> r_;

    // Tranpose of A, and corresponding mapping
    Sparsity AT_;
    std::vector<casadi_int> A_mapping_;

    // Aggregate SOCP helper constraints (lhs)
    IM map_Q_;

    // Aggregate SOCP helper constraints (rhs)
    std::vector<casadi_int> map_P_;

    // Maximum size of ind/val vectors
    casadi_int indval_size_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_GUROBI_INTERFACE_HPP
