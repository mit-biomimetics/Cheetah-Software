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


#ifndef CASADI_CSPARSE_INTERFACE_HPP
#define CASADI_CSPARSE_INTERFACE_HPP

/** \defgroup plugin_Linsol_csparse
 * Linsol with CSparse Interface
*/

/** \pluginsection{Linsol,csparse} */

/// \cond INTERNAL
#include <cs.h>
#include "casadi/core/linsol_internal.hpp"
#include <casadi/interfaces/csparse/casadi_linsol_csparse_export.h>

namespace casadi {
  struct CASADI_LINSOL_CSPARSE_EXPORT CsparseMemory : public LinsolMemory {
    // Destructor
    ~CsparseMemory();

    // The linear system CSparse form (CCS)
    cs A;

    // The symbolic factorization
    css *S;

    // The numeric factorization
    csn *N;

    // Temporary
    std::vector<double> temp_;

    std::vector<int> colind, row;
  };

  /** \brief \pluginbrief{LinsolInternal,csparse}
   * @copydoc LinsolInternal_doc
   * @copydoc plugin_LinsolInternal_csparse
   */
  class CASADI_LINSOL_CSPARSE_EXPORT CsparseInterface : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    CsparseInterface(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new LinsolInternal */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new CsparseInterface(name, sp);
    }

    // Destructor
    ~CsparseInterface() override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new CsparseMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<CsparseMemory*>(mem);}

    // Symbolic factorization
    int sfact(void* mem, const double* A) const override;

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    const char* plugin_name() const override { return "csparse";}

    // Get name of the class
    std::string class_name() const override { return "CsparseInterface";}
  };

} // namespace casadi

/// \endcond

#endif // CASADI_CSPARSE_INTERFACE_HPP
