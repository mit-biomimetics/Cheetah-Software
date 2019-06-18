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


#ifndef CASADI_JIT_FUNCTION_HPP
#define CASADI_JIT_FUNCTION_HPP

#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi {
  class CASADI_EXPORT JitFunction : public FunctionInternal {
  public:
    /** \brief Constructor */
    JitFunction(const std::string& name, const std::string& body,
              const std::vector<std::string>& name_in,
              const std::vector<std::string>& name_out,
              const std::vector<Sparsity>& sparsity_in,
              const std::vector<Sparsity>& sparsity_out);

    /** \brief Get type name */
    std::string class_name() const override { return "JitFunction";}

    /** \brief Destructor */
    ~JitFunction() override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override { return name_in_.size();}
    size_t get_n_out() override { return name_out_.size();}
    ///@}

    /** \brief Is codegen supported? */
    bool has_codegen() const override { return true;}

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    ///@{
    /** \brief Jacobian of all outputs with respect to all inputs */
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    // Buffer the function calls
    bool buffered_;

    // Function body
    std::string body_;

    // Jacobian function body
    std::string jac_body_;

    // Hessian function body
    std::string hess_body_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_JIT_FUNCTION_HPP
