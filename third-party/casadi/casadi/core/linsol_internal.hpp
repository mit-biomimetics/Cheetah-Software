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


#ifndef CASADI_LINSOL_INTERNAL_HPP
#define CASADI_LINSOL_INTERNAL_HPP

#include "linsol.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  struct CASADI_EXPORT LinsolMemory {
    // Current state of factorization
    bool is_sfact, is_nfact;

    // Constructor
    LinsolMemory() : is_sfact(false), is_nfact(false) {}
  };

  /** Internal class
      @copydoc Linsol_doc
  */
  class CASADI_EXPORT LinsolInternal
    : public ProtoFunction, public PluginInterface<LinsolInternal> {
  public:
    /// Constructor
    LinsolInternal(const std::string& name, const Sparsity& sp);

    /// Destructor
    ~LinsolInternal() override;

    /** \brief Display object */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Print more */
    virtual void disp_more(std::ostream& stream) const {}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LinsolMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LinsolMemory*>(mem);}

    /// Evaluate SX, possibly transposed
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res,
                                casadi_int* iw, SXElem* w, void* mem,
                                bool tr, casadi_int nrhs) const;

#if 0
    // (Re)factorize the system
    casadi_int factorize(void* mem, const double* A) const;

    // Needs symbolic factorization
    virtual bool needs_sfact(void* mem, const double* A) const;

    // Needs numeric factorization
    virtual bool needs_nfact(void* mem, const double* A) const;
#endif

    // Symbolic factorization
    virtual int sfact(void* mem, const double* A) const { return 0;}

    /// Numeric factorization
    virtual int nfact(void* mem, const double* A) const;

    // Solve numerically
    virtual int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const;

    /// Number of negative eigenvalues
    virtual casadi_int neig(void* mem, const double* A) const;

    /// Matrix rank
    virtual casadi_int rank(void* mem, const double* A) const;

    /// Generate C code
    virtual void generate(CodeGenerator& g, const std::string& A, const std::string& x,
                          casadi_int nrhs, bool tr) const;

    // Creator function for internal class
    typedef LinsolInternal* (*Creator)(const std::string& name, const Sparsity& sp);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    // Get name of the plugin
    const char* plugin_name() const override = 0;

    /// Get sparsity pattern
    casadi_int nrow() const { return sp_.size1();}
    casadi_int ncol() const { return sp_.size2();}
    const casadi_int* colind() const { return sp_.colind();}
    const casadi_int* row() const { return sp_.row();}
    casadi_int nnz() const { return sp_.nnz();}

    // Sparsity pattern of the linear system
    Sparsity sp_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_LINSOL_INTERNAL_HPP
