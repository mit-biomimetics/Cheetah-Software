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


#ifndef CASADI_SLICOT_DPLE_HPP
#define CASADI_SLICOT_DPLE_HPP

#include "../../core/dple_impl.hpp"
#include "../../core/linsol.hpp"
#include <casadi/interfaces/slicot/casadi_dple_slicot_export.h>

/** \defgroup plugin_Dple_slicot
 *
 * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT

 * Uses Periodic Schur Decomposition ('psd') and does not assume positive definiteness.
 * Based on Periodic Lyapunov equations: some applications and new algorithms.
 * Int. J. Control, vol. 67, pp. 69-87, 1997.
 *
 * Overview of the method:
 *   J. Gillis
 *   Practical Methods for Approximate Robust Periodic Optimal Control ofNonlinear Mechanical Systems,
 *   PhD Thesis, KULeuven, 2015
*/

/** \pluginsection{Dple,slicot} */

/// \cond INTERNAL
namespace casadi {


  // Forward declaration
  class SlicotDple;

  struct CASADI_DPLE_SLICOT_EXPORT SlicotDpleMemory {

    /// T Hessenberg-triangular data
    /// Z Schur form multiplier data
    /// X Schur form multiplier data
    // Xbar Schur form multiplier data
    /// VZ Transformed V data
    /// nnKa Temp data  (n x n) x K
    /// nnKb Temp data  (n x n) x K
    /// eig_real Real parts of eigenvalues
    /// eig_imag Imaginary parts of eigenvalues

    /// Temp data  F
    /// Temp data  FF
    /// dwork Work vector for periodic Schur form

    double *VZ, *T, *Z, *X, *Xbar, *nnKa, *nnKb, *eig_real, *eig_imag, *F, *FF, *A, *B;
    double *dwork, *wruntime;
    casadi_int* partition, *iwruntime;

    /// Solvers for low-order Discrete Periodic Sylvester Equations
    std::vector< std::vector< Linsol> > dpse_solvers;

    /// Constructor
    SlicotDpleMemory() {}

    /// Destructor
    ~SlicotDpleMemory() {}
  };

  /** \brief \pluginbrief{Dple,slicot}
   *
   * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT
   *
   * @copydoc Dple_doc
   * @copydoc plugin_Dple_slicot

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLE_SLICOT_EXPORT SlicotDple : public Dple {
  public:
    /** \brief  Constructor */
    explicit SlicotDple();

    /** \brief  Constructor
     * \param st \structargument{Dple}
     */
    SlicotDple(const std::string& name, const SpDict & st);

    /** \brief  Create a new QP Solver */
    static Dple* creator(const std::string& name,
                          const SpDict& st) {
      return new SlicotDple(name, st);
    }

    /** \brief  Destructor */
    ~SlicotDple() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "slicot";}

    // Get name of the class
    std::string class_name() const override { return "SlicotDple";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new SlicotDpleMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SlicotDpleMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /** \brief  Evaluate numerically */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    SlicotDple(const SpDict & st);

  private:
    /// Dimension of state-space
    casadi_int n_;

    inline casadi_int partindex(const SlicotDpleMemory* m, casadi_int i, casadi_int j, casadi_int k,
                                casadi_int r, casadi_int c) const;

    /// Numerical zero, used in periodic Schur form
    double psd_num_zero_;

    /// Linear solver name
    std::string linear_solver_;

    /// Options to be passed to linear solver constructor
    Dict linear_solver_options_;

    /// Has the plugin been loaded already?
    static bool has_loaded_;

  };


  void slicot_periodic_schur(casadi_int n, casadi_int K, const double* a,
                             double* t,  double * z,
                             double* dwork, double* eig_real,
                             double *eig_imag, double num_zero=0);
} // namespace casadi

/// \endcond
#endif // CASADI_SLICOT_DPLE_HPP
