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


#ifndef CASADI_COLLOCATION_HPP
#define CASADI_COLLOCATION_HPP

#include "casadi/core/integrator_impl.hpp"
#include "casadi/core/integration_tools.hpp"
#include <casadi/solvers/casadi_integrator_collocation_export.h>

/** \defgroup plugin_Integrator_collocation

     Fixed-step implicit Runge-Kutta integrator
     ODE/DAE integrator based on collocation schemes

     The method is still under development

*/

/** \pluginsection{Integrator,collocation} */

/// \cond INTERNAL
namespace casadi {

  /**
     \brief \pluginbrief{Integrator,collocation}

     @copydoc DAE_doc
     @copydoc plugin_Integrator_collocation

     \author Joel Andersson
     \date 2014
  */
  class CASADI_INTEGRATOR_COLLOCATION_EXPORT Collocation :
        public ImplicitFixedStepIntegrator {
  public:

    /// Constructor
    explicit Collocation(const std::string& name, const Function& dae);

    /** \brief  Create a new integrator */
    static Integrator* creator(const std::string& name, const Function& dae) {
      return new Collocation(name, dae);
    }

    /// Destructor
    ~Collocation() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "collocation";}

    // Get name of the class
    std::string class_name() const override { return "Collocation";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Initialize stage
    void init(const Dict& opts) override;

    /// Setup F and G
    void setupFG() override;

    // Return zero if smaller than machine epsilon
    static double zeroIfSmall(double x);

    /** \brief Reset the forward problem */
    void reset(IntegratorMemory* mem, double t, const double* x,
                       const double* z, const double* p) const override;

    /// Reset the backward problem and take time to tf
    void resetB(IntegratorMemory* mem, double t, const double* rx,
                        const double* rz, const double* rp) const override;

    // Interpolation order
    casadi_int deg_;

    // Collocation scheme
    std::string collocation_scheme_;

    /// A documentation string
    static const std::string meta_doc;

    /// Continuous time dynamics
    Function f_, g_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_COLLOCATION_HPP
