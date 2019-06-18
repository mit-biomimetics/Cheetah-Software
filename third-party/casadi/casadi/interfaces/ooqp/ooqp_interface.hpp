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


#ifndef CASADI_OOQP_INTERFACE_HPP
#define CASADI_OOQP_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/ooqp/casadi_conic_ooqp_export.h>

/** \defgroup plugin_Conic_ooqp
 Interface to the OOQP Solver for quadratic programming
  The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.

  NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
*/

/** \pluginsection{Conic,ooqp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Conic,ooqp}

      @copydoc Conic_doc
      @copydoc plugin_Conic_ooqp

  */
  class CASADI_CONIC_OOQP_EXPORT OoqpInterface : public Conic {
  public:
    /** \brief  Create a new Solver */
    explicit OoqpInterface(const std::string& name,
                           const std::map<std::string, Sparsity>& st);

    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                                     const std::map<std::string, Sparsity>& st) {
      return new OoqpInterface(name, st);
    }

    /** \brief  Destructor */
    ~OoqpInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "ooqp";}

    // Get name of the class
    std::string class_name() const override { return "OoqpInterface";}

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /// Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /// Throw error
    static const char* errFlag(int flag);

    /// Print an OOQP bounds vector
    static std::string printBounds(const std::vector<double>& b,
                                   const std::vector<char>& ib, casadi_int n, const char *sign);

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Transpose of linear constraints
    Sparsity spAT_;

    // Number of nonzeros in upper triangular half of Hessian
    casadi_int nQ_;

    // Number of nonzeros in Hessian
    casadi_int nH_;

    // Number of nonzeros in constraint matrix
    casadi_int nA_;

    // Print level
    casadi_int print_level_;

    // Tolerances
    double mutol_, artol_;

    /// A documentation string
    static const std::string meta_doc;

    mutable int return_status_;
    mutable bool success_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_OOQP_INTERFACE_HPP
