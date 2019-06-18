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


#ifndef CASADI_IDAS_INTERFACE_HPP
#define CASADI_IDAS_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_integrator_idas_export.h>
#include "sundials_interface.hpp"
#include <idas/idas.h>            /* prototypes for CVODE fcts. and consts. */
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Integrator_idas
      Interface to IDAS from the Sundials suite.
*/

/** \pluginsection{Integrator,idas} */

/// \cond INTERNAL

namespace casadi {
  // Forward declaration
  class IdasInterface;

  // IdasMemory
  struct CASADI_INTEGRATOR_IDAS_EXPORT IdasMemory : public SundialsMemory {
    /// Function object
    const IdasInterface& self;

    // Idas memory block
    void* mem;

    // Ids of backward problem
    int whichB;

    /// Constructor
    IdasMemory(const IdasInterface& s);

    /// Destructor
    ~IdasMemory();
  };

  /** \brief \pluginbrief{Integrator,idas}

      @copydoc IdasIntegrator_doc
      @copydoc plugin_Integrator_idas

      \author Joel Andersson
      \date 2010
  */
  class CASADI_INTEGRATOR_IDAS_EXPORT IdasInterface : public SundialsInterface {

  public:

    /** \brief  Constructor */
    explicit IdasInterface(const std::string& name, const Function& dae);

    /** \brief  Create a new integrator */
    static Integrator* creator(const std::string& name, const Function& dae) {
      return new IdasInterface(name, dae);
    }

    /** \brief  Destructor */
    ~IdasInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "idas";}

    // Get name of the class
    std::string class_name() const override { return "IdasInterface";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Initialize the taping */
    void initTaping(IdasMemory* m) const;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new IdasMemory(*this);}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<IdasMemory*>(mem);}

    /** \brief  Reset the forward problem and bring the time back to t0 */
    void reset(IntegratorMemory* mem, double t, const double* x,
                       const double* z, const double* p) const override;

    /** \brief  Advance solution in time */
    void advance(IntegratorMemory* mem, double t, double* x,
                         double* z, double* q) const override;

    /** \brief  Reset the backward problem and take time to tf */
    void resetB(IntegratorMemory* mem, double t, const double* rx,
                        const double* rz, const double* rp) const override;

    /** \brief  Retreat solution in time */
    void retreat(IntegratorMemory* mem, double t, double* rx,
                         double* rz, double* rq) const override;

    /** \brief  Set the stop time of the forward integration */
    void setStopTime(IntegratorMemory* mem, double tf) const override;

    /** \brief Cast to memory object */
    static IdasMemory* to_mem(void *mem) {
      IdasMemory* m = static_cast<IdasMemory*>(mem);
      casadi_assert_dev(m);
      return m;
    }

    ///@{
    // Get system Jacobian
    Function getJ(bool backward) const override;
    template<typename MatType> Function getJ(bool backward) const;
    ///@}

    /// A documentation string
    static const std::string meta_doc;

  protected:

    // Sundials callback functions
    static int res(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, void *user_data);
    static int resB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                    N_Vector rrB, void *user_data);
    static void ehfun(int error_code, const char *module, const char *function, char *msg,
                      void *eh_data);
    static int jtimes(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v,
                      N_Vector Jv, double cj, void *user_data,
                      N_Vector tmp1, N_Vector tmp2);
    static int jtimesB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                       N_Vector resvalB, N_Vector vB, N_Vector JvB, double cjB,
                       void *user_data, N_Vector tmp1B, N_Vector tmp2B);
    static int rhsQ(double t, N_Vector xz, N_Vector xzdot, N_Vector qdot, void *user_data);
    static int rhsQB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                     N_Vector qdotA, void *user_data);
    static int psolve(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec,
                      N_Vector zvec, double cj, double delta, void *user_data,
                      N_Vector tmp);
    static int psetup(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj,
                      void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int psolveB(double t,
                       N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                       N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, double cjB,
                       double deltaB, void *user_dataB, N_Vector tmpB);
    static int psetupB(double t,
                       N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                       N_Vector resvalB, double cjB, void *user_dataB,
                       N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int lsetup(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
                      N_Vector xzdotcur, N_Vector rescur);
    static int lsetupB(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                       N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolveB(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
                       N_Vector xzdotcur, N_Vector rescur);
  public:

    // Throw error
    static void idas_error(const char* module, int flag);

    // Options
    bool cj_scaling_;
    bool calc_ic_;
    bool calc_icB_;
    bool suppress_algebraic_;
    double max_step_size_;
    std::vector<double> abstolv_;
    double first_time_;

    //  Initial values for \p xdot and \p z
    std::vector<double> init_xdot_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_IDAS_INTERFACE_HPP
