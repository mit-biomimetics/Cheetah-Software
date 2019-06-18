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

#ifndef CASADI_CLP_INTERFACE_HPP
#define CASADI_CLP_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/clp/casadi_conic_clp_export.h>

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpEventHandler.hpp"

#include <string>

/** \defgroup plugin_Conic_clp

      Interface to Clp solver for sparse Quadratic Programs
*/

/** \pluginsection{Conic,clp} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_CLP_EXPORT ClpMemory : public ConicMemory {
    /// Constructor
    ClpMemory();

    /// Destructor
    ~ClpMemory();

    std::vector<int> colind, row;

    int return_status;
    int secondary_return_status;
    bool success;

  };

  /** \brief \pluginbrief{Conic,clp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_clp

      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_CONIC_CLP_EXPORT ClpInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new ClpInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit ClpInterface(const std::string& name,
                            const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~ClpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "clp";}

    // Get name of the class
    std::string class_name() const override { return "ClpInterface";}

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new ClpMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<ClpMemory*>(mem);}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CLP_INTERFACE_HPP
