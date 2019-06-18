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


#ifndef CASADI_FAST_NEWTON_HPP
#define CASADI_FAST_NEWTON_HPP

#include "casadi/core/rootfinder_impl.hpp"
#include <casadi/solvers/casadi_rootfinder_fast_newton_export.h>

/** \defgroup plugin_Rootfinder_fast_newton
     Implements simple newton iterations to solve an implicit function.
*/

/** \pluginsection{Rootfinder,fast_newton} */

/// \cond INTERNAL
namespace casadi {

  // Memory
  struct CASADI_ROOTFINDER_FAST_NEWTON_EXPORT FastNewtonMemory
    : public RootfinderMemory {
    // Return status
    int return_status;
    // Number of iterations
    casadi_int iter;

    casadi_newton_mem<double> M;
  };

  /** \brief \pluginbrief{Rootfinder,fast_newton}

      @copydoc Rootfinder_doc
      @copydoc plugin_Rootfinder_fast_newton

      \author Joris Gillis
      \date 2018
  */
  class CASADI_ROOTFINDER_FAST_NEWTON_EXPORT FastNewton : public Rootfinder {
  public:
    /** \brief  Constructor */
    explicit FastNewton(const std::string& name, const Function& f);

    /** \brief  Destructor */
    ~FastNewton() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "fast_newton";}

    // Name of the class
    std::string class_name() const override { return "FastNewton";}

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new FastNewton(name, f);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new FastNewtonMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<FastNewtonMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// Solve the system of equations and calculate derivatives
    int solve(void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

  protected:
    /// Maximum number of Newton iterations
    casadi_int max_iter_;

    /// Absolute tolerance that should be met on residual
    double abstol_;

    /// Absolute tolerance that should be met on step
    double abstolStep_;

    /// Reference to jacobian function
    Function jac_f_z_;

    /// Data for qr
    Sparsity sp_v_;
    Sparsity sp_r_;
    std::vector<casadi_int> prinv_;
    std::vector<casadi_int> pc_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_FAST_NEWTON_HPP
