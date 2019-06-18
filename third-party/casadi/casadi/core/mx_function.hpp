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


#ifndef CASADI_MX_FUNCTION_HPP
#define CASADI_MX_FUNCTION_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "x_function.hpp"
#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {

#ifndef SWIG
  /** \brief  An element of the algorithm, namely an MX node */
  struct MXAlgEl {
    /// Operator index
    casadi_int op;

    /// Data associated with the operation
    MX data;

    /// Work vector indices of the arguments
    std::vector<casadi_int> arg;

    /// Work vector indices of the results
    std::vector<casadi_int> res;
  };
#endif // SWIG

  /** \brief  Internal node class for MXFunction
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT MXFunction :
        public XFunction<MXFunction, MX, MXNode>{
  public:
    /** \brief  An element of the algorithm, namely an MX node */
    typedef MXAlgEl AlgEl;

    /** \brief  All the runtime elements in the order of evaluation */
    std::vector<AlgEl> algorithm_;

    /** \brief Offsets for elements in the w_ vector */
    std::vector<casadi_int> workloc_;

    /// Free variables
    std::vector<MX> free_vars_;

    /// Default input values
    std::vector<double> default_in_;

    /** \brief Constructor */
    MXFunction(const std::string& name,
      const std::vector<MX>& input, const std::vector<MX>& output,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

    /** \brief  Destructor */
    ~MXFunction() override;

    /** \brief  Evaluate numerically, work vectors given */
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

    /** \brief  Print description */
    void disp_more(std::ostream& stream) const override;

    /** \brief Get type name */
    std::string class_name() const override {return "MXFunction";}

    /** \brief Check if the function is of a particular type */
    bool is_a(const std::string& type, bool recursive) const override;

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /** \brief Codegen incref for dependencies */
    void codegen_incref(CodeGenerator& g) const override;

    /** \brief Codegen decref for dependencies */
    void codegen_decref(CodeGenerator& g) const override;

    /** \brief Generate code for the body of the C function */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Extract the residual function G and the modified function Z out of an expression
     * (see Albersmeyer2010 paper) */
    void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const override;

    /** Inline calls? */
    bool should_inline(bool always_inline, bool never_inline) const override;

    /** \brief Evaluate symbolically, SX type*/
    int eval_sx(const SXElem** arg, SXElem** res,
                casadi_int* iw, SXElem* w, void* mem) const override;

    /** \brief Evaluate symbolically, MX type */
    void eval_mx(const MXVector& arg, MXVector& res,
                 bool always_inline, bool never_inline) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fwdSeed,
                        std::vector<std::vector<MX> >& fwdSens) const;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& adjSeed,
                        std::vector<std::vector<MX> >& adjSens) const;

    /// Get a vector of symbolic variables corresponding to the outputs
    std::vector<MX> symbolic_output(const std::vector<MX>& arg) const override;

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res,
                  casadi_int* iw, bvec_t* w, void* mem) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

    // print an element of an algorithm
    std::string print(const AlgEl& el) const;

    ///@{
    /** \brief Get function input(s) and output(s)  */
    const MX mx_in(casadi_int ind) const override;
    const std::vector<MX> mx_in() const override;
    ///@}

    /// Get free variables (MX)
    std::vector<MX> free_mx() const override {return free_vars_;}

    /** \brief Does the function have free variables */
    bool has_free() const override { return !free_vars_.empty();}

    /** \brief Print free variables */
    std::vector<std::string> get_free() const override {
      std::vector<std::string> ret;
      for (auto&& e : free_vars_) ret.push_back(e.name());
      return ret;
    }

    /** \brief Number of nodes in the algorithm */
    casadi_int n_nodes() const override { return algorithm_.size();}

    casadi_int n_instructions() const override { return algorithm_.size();}

    /** *\brief get MX expression associated with instruction */
    MX instruction_MX(casadi_int k) const override;

    /** \brief Get an atomic operation operator index */
    casadi_int instruction_id(casadi_int k) const override { return algorithm_.at(k).op;}

    /** \brief Get default input value */
    double get_default_in(casadi_int ind) const override { return default_in_.at(ind);}

    /** \brief Get the (integer) input arguments of an atomic operation */
    std::vector<casadi_int> instruction_input(casadi_int k) const override;

    /** \brief Get the (integer) output argument of an atomic operation */
    std::vector<casadi_int> instruction_output(casadi_int k) const override;

    /** \brief Export function in a specific language */
    void export_code_body(const std::string& lang,
      std::ostream &stream, const Dict& options) const override;

    /// Substitute inplace, internal implementation
    void substitute_inplace(std::vector<MX>& vdef, std::vector<MX>& ex) const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MX_FUNCTION_HPP
