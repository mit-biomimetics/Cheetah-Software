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


#ifndef CASADI_EXTERNAL_IMPL_HPP
#define CASADI_EXTERNAL_IMPL_HPP

#include "external.hpp"
#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {
  class CASADI_EXPORT External : public FunctionInternal {
  protected:
    /** \brief Information about the library */
    Importer li_;

    /** \brief Increase/decrease reference counter */
    signal_t incref_, decref_;

    /** \brief Number of inputs and outputs */
    getint_t get_n_in_, get_n_out_;

    /** \brief Names of inputs and outputs */
    name_t get_name_in_, get_name_out_;

    /** \brief Work vector sizes */
    work_t work_;

    ///@{
    /** \brief Data vectors */
    std::vector<casadi_int> int_data_;
    std::vector<double> real_data_;
    std::string string_data_;
    ///@}
  public:

    /** \brief Constructor */
    External(const std::string& name, const Importer& li);

    /** \brief Destructor */
    ~External() override = 0;

    // Factory
    Function factory(const std::string& name,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const Function::AuxOut& aux,
                             const Dict& opts) const override;

    /** \brief Get type name */
    std::string class_name() const override { return "External";}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    ///@{
    /** \brief Number of function inputs and outputs */
    size_t get_n_in() override;
    size_t get_n_out() override;
    ///@}

    ///@{
    /** \brief Names of function input and outputs */
    std::string get_name_in(casadi_int i) override;
    std::string get_name_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Forward mode derivatives */
    Function get_forward(casadi_int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const override;
    bool has_forward(casadi_int nfwd) const override;
    ///@}

    ///@{
    /** \brief Reverse mode derivatives */
    Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const override;
    bool has_reverse(casadi_int nadj) const override;
    ///@}

    ///@{
    /** \brief Full Jacobian */
    bool has_jacobian() const override;
    Function get_jacobian(const std::string& name,
                                     const std::vector<std::string>& inames,
                                     const std::vector<std::string>& onames,
                                     const Dict& opts) const override;
    ///@}
  };

  class CASADI_EXPORT GenericExternal : public External {
    // Sparsities
    sparsity_t get_sparsity_in_, get_sparsity_out_;

    // Memory allocation
    alloc_mem_t alloc_mem_;
    init_mem_t init_mem_;
    free_mem_t free_mem_;

  public:
    /** \brief Constructor */
    GenericExternal(const std::string& name, const Importer& li);

    /** \brief  Destructor */
    ~GenericExternal() override { this->clear_mem();}

    /// Initialize
    void init(const Dict& opts) override;

    /// @{
    /** \brief Retreive sparsities */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    /** \brief Create memory block */
    void* alloc_mem() const override;

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_EXTERNAL_IMPL_HPP
