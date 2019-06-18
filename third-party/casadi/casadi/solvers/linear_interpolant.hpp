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


#ifndef CASADI_LINEAR_INTERPOLANT_HPP
#define CASADI_LINEAR_INTERPOLANT_HPP

#include "casadi/core/interpolant_impl.hpp"
#include <casadi/solvers/casadi_interpolant_linear_export.h>

/** \defgroup plugin_Interpolant_linear
*/

/** \pluginsection{Interpolant,linear} */

/// \cond INTERNAL

namespace casadi {
  /** \brief \pluginbrief{Interpolant,linear}
    Implements a multilinear interpolant: For 1D, the interpolating polynomial
    will be linear. For 2D, the interpolating polynomial will be bilinear, etc.

    @copydoc Interpolant_doc
    @copydoc plugin_Interpolant_linear
    \author Joel Andersson
    \date 2016
  */
  class CASADI_INTERPOLANT_LINEAR_EXPORT LinearInterpolant : public Interpolant {
  public:
    // Constructor
    LinearInterpolant(const std::string& name,
                      const std::vector<double>& grid,
                      const std::vector<casadi_int>& offset,
                      const std::vector<double>& values,
                      casadi_int m);

    // Destructor
    ~LinearInterpolant() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "linear";}

    // Name of the class
    std::string class_name() const override { return "LinearInterpolant";}

    /** \brief  Create a new Interpolant */
    static Interpolant* creator(const std::string& name,
                                const std::vector<double>& grid,
                                const std::vector<casadi_int>& offset,
                                const std::vector<double>& values,
                                casadi_int m) {
      return new LinearInterpolant(name, grid, offset, values, m);
    }

    // Initialize
    void init(const Dict& opts) override;

    /// Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    ///@{
    /** \brief Full Jacobian */
    bool has_jacobian() const override { return true;}
    Function get_jacobian(const std::string& name,
                                      const std::vector<std::string>& inames,
                                      const std::vector<std::string>& onames,
                                      const Dict& opts) const override;
    ///@}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    /// A documentation string
    static const std::string meta_doc;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    std::vector<casadi_int> lookup_mode_;
  };

  /** First order derivatives */
  class CASADI_INTERPOLANT_LINEAR_EXPORT LinearInterpolantJac : public FunctionInternal {
  public:
    /// Constructor
    LinearInterpolantJac(const std::string& name) : FunctionInternal(name) {}

    /// Destructor
    ~LinearInterpolantJac() override {}

    /** \brief Get type name */
    std::string class_name() const override { return "LinearInterpolantJac";}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    // Initialize
    void init(const Dict& opts) override;

    /// Evaluate numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    ///@{
    /** \brief Full Jacobian */
    bool has_jacobian() const override { return true;}
    Function get_jacobian(const std::string& name,
                                      const std::vector<std::string>& inames,
                                      const std::vector<std::string>& onames,
                                      const Dict& opts) const override;
    ///@}

  };

} // namespace casadi

/// \endcond
#endif // CASADI_LINEAR_INTERPOLANT_HPP
