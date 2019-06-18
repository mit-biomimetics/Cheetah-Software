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


#ifndef CASADI_SYMBOLIC_QR_HPP
#define CASADI_SYMBOLIC_QR_HPP

#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_symbolicqr_export.h>

/** \defgroup plugin_Linsol_symbolicqr

       Linsol based on QR factorization with sparsity pattern based reordering
      _without_ partial pivoting
*/

/** \pluginsection{Linsol,symbolicqr} */

/// \cond INTERNAL

namespace casadi {
  typedef SX* SXPtr;
  typedef std::vector<SXPtr> SXPtrV;

  /** \brief Memory for SymbolicQR  */
  struct CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQrMemory : public LinsolMemory {
    // Work vectors
    std::vector<const double*> arg;
    std::vector<double*> res;
    std::vector<casadi_int> iw;
    std::vector<double> w;

    // Allocate memory for a function
    void alloc(const Function& f);

    // Storage for QR factorization
    std::vector<double> q, r;
  };

  /** \brief \pluginbrief{Linsol,symbolicqr}

      @copydoc Linsol_doc
      @copydoc plugin_Linsol_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQr : public LinsolInternal {
  public:
    // Constructor
    SymbolicQr(const std::string& name, const Sparsity& sp);

    // Destructor
    ~SymbolicQr() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "symbolicqr";}

    // Name of the class
    std::string class_name() const override { return "SymbolicQr";}

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new SymbolicQr(name, sp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SymbolicQrMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SymbolicQrMemory*>(mem);}

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /** \brief Evaluate symbolically (SX) */
    void linsol_eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem,
                        bool tr, casadi_int nrhs) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Functions for factorization and (optionally transposed) solve
    Function factorize_, solve_, solveT_;

    // Generated function options
    Dict fopts_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SYMBOLIC_QR_HPP
