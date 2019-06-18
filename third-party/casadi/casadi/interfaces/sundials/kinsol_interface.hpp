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


#ifndef CASADI_KINSOL_INTERFACE_HPP
#define CASADI_KINSOL_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_rootfinder_kinsol_export.h>
#include "casadi/core/rootfinder_impl.hpp"
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h>  /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>  /* definition of type double */
#include <kinsol/kinsol.h>            /* prototypes for CVode fcts. and consts. */
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_impl.h> /* Needed for the provided linear solver */
#include <kinsol/kinsol_spils_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Rootfinder_kinsol
 KINSOL interface from the Sundials suite
*/
/** \pluginsection{Rootfinder,kinsol} */

/// \cond INTERNAL
namespace casadi {
  /// Forward declaration
  class KinsolInterface;

  // Memory
  struct CASADI_ROOTFINDER_KINSOL_EXPORT KinsolMemory : public RootfinderMemory {
    /// Function object
    const KinsolInterface& self;

    /// Constructor
    KinsolMemory(const KinsolInterface& s);

    /// Destructor
    ~KinsolMemory();

    /// KINSOL memory block
    void* mem;

    /// Variable
    N_Vector u;

    // Current Jacobian
    double* jac;
  };

  /** \brief \pluginbrief{Rootfinder,kinsol}
  *
  * @copydoc Rootfinder_doc
  * @copydoc plugin_Rootfinder_kinsol
  */
  class CASADI_ROOTFINDER_KINSOL_EXPORT KinsolInterface : public Rootfinder {
  public:
    /** \brief  Constructor */
    explicit KinsolInterface(const std::string& name, const Function& f);

    /** \brief  Destructor */
    ~KinsolInterface() override;

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new KinsolInterface(name, f);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize stage */
    void init(const Dict& opts) override;

    /// Solve the system of equations and calculate derivatives
    int solve(void* mem) const override;

    // Get name of the plugin
    const char* plugin_name() const override { return "kinsol";}

    // Get name of the class
    std::string class_name() const override { return "KinsolInterface";}

    // Scaling
    N_Vector u_scale_, f_scale_;

    /// Globalization strategy
    casadi_int strategy_;

    // Should KINSOL internal warning messages be ignored
    bool disable_internal_warnings_;

    // Maximum number of iterations
    casadi_int max_iter_;

    // Use exact Jacobian?
    bool exact_jac_;

    // Type of linear solver
    enum LinsolType { DENSE, BANDED, ITERATIVE, USER_DEFINED};
    LinsolType linear_solver_type_;

    // Bandwidth (for banded solvers)
    casadi_int upper_bandwidth_, lower_bandwidth_;

    // Krylov subspace size (for iterative solvers)
    casadi_int maxl_;

    // Iterative solver
    enum IterativeSolver { GMRES, BCGSTAB, TFQMR};
    IterativeSolver iterative_solver_;

    // Should a preconditioner be used
    bool use_preconditioner_;

    // Absolute tolerance
    double abstol_;

    // Jacobian times vector function
    Function jtimes_;

    // Get jtimes_
    void get_jtimes();

    // Raise an error specific to KinSol
    void kinsol_error(const std::string& module, int flag, bool fatal=true) const;

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new KinsolMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<KinsolMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /** \brief Cast to memory object */
    static KinsolMemory* to_mem(void *mem) {
      KinsolMemory* m = static_cast<KinsolMemory*>(mem);
      casadi_assert_dev(m);
      return m;
    }

    /** \brief Callback functions (to be updated) */
    void func(KinsolMemory& m, N_Vector u, N_Vector fval) const;
    void djac(KinsolMemory& m, long N, N_Vector u, N_Vector fu,
              DlsMat J, N_Vector tmp1, N_Vector tmp2) const;
    void bjac(KinsolMemory& m, long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J,
              N_Vector tmp1, N_Vector tmp2) const;
    void jtimes(KinsolMemory& m, N_Vector v, N_Vector Jv, N_Vector u, int* new_u) const;
    void psetup(KinsolMemory& m, N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                N_Vector tmp1, N_Vector tmp2) const;
    void psolve(KinsolMemory& m, N_Vector u, N_Vector uscale,
                N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp) const;

    /** \brief Wrappers to callback functions*/
    static int func_wrapper(N_Vector u, N_Vector fval, void *user_data);
    static int djac_wrapper(long N, N_Vector u, N_Vector fu, DlsMat J, void *user_data,
                            N_Vector tmp1, N_Vector tmp2);
    static int bjac_wrapper(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J,
                            void *user_data, N_Vector tmp1, N_Vector tmp2);
    static int jtimes_wrapper(N_Vector v, N_Vector Jv, N_Vector u, int* new_u, void *user_data);
    static int psetup_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                              void* user_data, N_Vector tmp1, N_Vector tmp2);
    static int psolve_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                              N_Vector v, void* user_data, N_Vector tmp);

    /** \brief Callback functions (updated) */
    static int lsetup(KINMem kin_mem);
    static int lsolve(KINMem kin_mem, N_Vector x, N_Vector b, double *sJpnorm, double *sFdotJp);
    static void ehfun(int error_code, const char *module, const char *function,
                      char *msg, void *eh_data);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_KINSOL_INTERFACE_HPP
