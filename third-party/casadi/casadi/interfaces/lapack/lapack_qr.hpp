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


#ifndef CASADI_LAPACK_QR_HPP
#define CASADI_LAPACK_QR_HPP

#include "casadi/core/linsol_internal.hpp"
#include <casadi/interfaces/lapack/casadi_linsol_lapackqr_export.h>

extern "C" {
  /// QR-factorize dense matrix (lapack)
  void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
               double *work, int *lwork, int *info);

  /// Multiply right hand side with Q-transpose (lapack)
  void dormqr_(char *side, char *trans, int *n, int *m, int *k, double *a,
               int *lda, double *tau, double *c, int *ldc,
               double *work, int *lwork, int *info);

  /// Solve upper triangular system (lapack)
  void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
                         double *alpha, double *a, int *lda, double *b, int *ldb);
}

/** \defgroup plugin_Linsol_lapackqr
*
* This class solves the linear system <tt>A.x=b</tt> by making an QR factorization of A: \n
* <tt>A = Q.R</tt>, with Q orthogonal and R upper triangular
*/

/** \pluginsection{Linsol,lapackqr} */

/// \cond INTERNAL
namespace casadi {
  struct CASADI_LINSOL_LAPACKQR_EXPORT LapackQrMemory : public LinsolMemory {
    // Matrix
    std::vector<double> mat;

    // The scalar factors of the elementary reflectors
    std::vector<double> tau;

    // qr work array
    std::vector<double> work;
  };

  /** \brief  \pluginbrief{Linsol,lapackqr}
   *
   @copydoc Linsol_doc
   @copydoc plugin_Linsol_lapackqr
   *
   */
  class CASADI_LINSOL_LAPACKQR_EXPORT LapackQr : public LinsolInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackQr(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new LapackQr(name, sp);
    }

    // Destructor
    ~LapackQr() override;

    // Initialize the solver
    void init(const Dict& opts) override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LapackQrMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LapackQrMemory*>(mem);}

    // Factorize the linear system
    int nfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve_batch(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    const char* plugin_name() const override { return "lapackqr";}

    // Get name of the class
    std::string class_name() const override { return "LapackQr";}

    // Maximum number of right-hand-sides
    casadi_int max_nrhs_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_LAPACK_QR_HPP
