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

#ifndef CASADI_CALLBACK_INTERNAL_HPP
#define CASADI_CALLBACK_INTERNAL_HPP

#include "callback.hpp"
#include "function_internal.hpp"

namespace casadi {

  class CASADI_EXPORT CallbackInternal : public FunctionInternal {
    friend class CallbackFunction;
  public:

    /** \brief Constructor */
    explicit CallbackInternal(const std::string& name, Callback* self);

    /** \brief Destructor */
    ~CallbackInternal() override;

    /** \brief Get type name */
    std::string class_name() const override {return "CallbackInternal";}

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    /// @}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Finalize the object creation */
    void finalize(const Dict& opts) override;

    /** \brief  Evaluate numerically, work vectors given */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief  Evaluate symbolically, work vectors given */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw,
                SXElem* w, void* mem) const override;

    /** \brief Evaluate with DM matrices */
    std::vector<DM> eval_dm(const std::vector<DM>& arg) const override;

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    bool uses_output() const override;

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                          const std::vector<std::string>& inames,
                          const std::vector<std::string>& onames,
                          const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives */
    bool has_forward(casadi_int nfwd) const override;
    Function get_forward(casadi_int nfwd, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives */
    bool has_reverse(casadi_int nadj) const override;
    Function get_reverse(casadi_int nadj, const std::string& name,
                         const std::vector<std::string>& inames,
                         const std::vector<std::string>& onames,
                         const Dict& opts) const override;
    ///@}

    /** \brief Pointer to the public class */
    Callback* self_;
  };

} // namespace casadi

#endif // CASADI_CALLBACK_INTERNAL_HPP
