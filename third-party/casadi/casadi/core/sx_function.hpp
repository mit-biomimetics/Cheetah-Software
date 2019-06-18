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


#ifndef CASADI_SX_FUNCTION_HPP
#define CASADI_SX_FUNCTION_HPP

#include "x_function.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief  An atomic operation for the SXElem virtual machine */
  struct ScalarAtomic {
    int op;     /// Operator index
    int i0;
    union {
      double d;
      struct { int i1, i2; };
    };
  };

/** \brief  Internal node class for SXFunction
    Do not use any internal class directly - always use the public Function
    \author Joel Andersson
    \date 2010-2015
*/
class CASADI_EXPORT SXFunction :
        public XFunction<SXFunction, Matrix<SXElem>, SXNode>{
  public:
    /** \brief Constructor */
    SXFunction(const std::string& name,
               const std::vector<Matrix<SXElem> >& inputv,
               const std::vector<Matrix<SXElem> >& outputv,
               const std::vector<std::string>& name_in,
               const std::vector<std::string>& name_out);

  /** \brief  Destructor */
  ~SXFunction() override;

  /** \brief  Evaluate numerically, work vectors given */
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  /** \brief  evaluate symbolically while also propagating directional derivatives */
  int eval_sx(const SXElem** arg, SXElem** res,
              casadi_int* iw, SXElem* w, void* mem) const override;

  /** Inline calls? */
  bool should_inline(bool always_inline, bool never_inline) const override {
    return true;
  }

  /** \brief Calculate forward mode directional derivatives */
  void ad_forward(const std::vector<std::vector<SX> >& fseed,
                            std::vector<std::vector<SX> >& fsens) const;

  /** \brief Calculate reverse mode directional derivatives */
  void ad_reverse(const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens) const;

  /** \brief  Check if smooth */
  bool is_smooth() const;

  /** \brief  Print the algorithm */
  void disp_more(std::ostream& stream) const override;

  /** \brief Get type name */
  std::string class_name() const override {return "SXFunction";}

  /** \brief Check if the function is of a particular type */
  bool is_a(const std::string& type, bool recursive) const override;

  ///@{
  /** \brief Get function input(s) and output(s)  */
  const SX sx_in(casadi_int ind) const override;
  const std::vector<SX> sx_in() const override;
  ///@}

  /// Get free variables (SX)
  std::vector<SX> free_sx() const override {
    std::vector<SX> ret(free_vars_.size());
    std::copy(free_vars_.begin(), free_vars_.end(), ret.begin());
    return ret;
  }

  /** \brief Does the function have free variables */
  bool has_free() const override { return !free_vars_.empty();}

  /** \brief Print free variables */
  std::vector<std::string> get_free() const override {
    std::vector<std::string> ret;
    for (auto&& e : free_vars_) ret.push_back(e.name());
    return ret;
  }

  /** \brief Hessian (forward over adjoint) via source code transformation */
  SX hess(casadi_int iind=0, casadi_int oind=0);

  /** \brief Get the number of atomic operations */
  casadi_int n_instructions() const override { return algorithm_.size();}

  /** \brief Get an atomic operation operator index */
  casadi_int instruction_id(casadi_int k) const override { return algorithm_.at(k).op;}

  /** \brief Get the (integer) input arguments of an atomic operation */
  std::vector<casadi_int> instruction_input(casadi_int k) const override {
    auto e = algorithm_.at(k);
    if (casadi_math<double>::ndeps(e.op)==2 || e.op==OP_INPUT) {
      return {e.i1, e.i2};
    } else if (casadi_math<double>::ndeps(e.op)==1) {
      return {e.i1};
    } else {
      return {};
    }
  }

  /** \brief Get the floating point output argument of an atomic operation */
  double instruction_constant(casadi_int k) const override {
    return algorithm_.at(k).d;
  }

  /** \brief Get the (integer) output argument of an atomic operation */
  std::vector<casadi_int> instruction_output(casadi_int k) const override {
    auto e = algorithm_.at(k);
    if (e.op==OP_OUTPUT) {
      return {e.i0, e.i2};
    } else {
      return {e.i0};
    }
  }

  /** \brief Number of nodes in the algorithm */
  casadi_int n_nodes() const override { return algorithm_.size() - nnz_out();}

  /** \brief  DATA MEMBERS */

  /** \brief  An element of the algorithm, namely a binary operation */
  typedef ScalarAtomic AlgEl;

  /** \brief  An element of the tape */
  template<typename T>
  struct TapeEl {
    T d[2];
  };

  /** \brief  all binary nodes of the tree in the order of execution */
  std::vector<AlgEl> algorithm_;

  // Work vector size
  size_t worksize_;

  /// Free variables
  std::vector<SXElem> free_vars_;

  /// The expressions corresponding to each binary operation
  std::vector<SXElem> operations_;

  /// The expressions corresponding to each constant
  std::vector<SXElem> constants_;

  /// Default input values
  std::vector<double> default_in_;

  ///@{
  /** \brief Options */
  static Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /** \brief  Initialize */
  void init(const Dict& opts) override;

  /** \brief Generate code for the declarations of the C function */
  void codegen_declarations(CodeGenerator& g) const override;

  /** \brief Generate code for the body of the C function */
  void codegen_body(CodeGenerator& g) const override;

  /** \brief  Propagate sparsity forward */
  int sp_forward(const bvec_t** arg, bvec_t** res,
                  casadi_int* iw, bvec_t* w, void* mem) const override;

  /** \brief  Propagate sparsity backwards */
  int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const override;

  /** \brief Return Jacobian of all input elements with respect to all output elements */
  Function get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const override;

  /** \brief Get default input value */
  double get_default_in(casadi_int ind) const override { return default_in_.at(ind);}

  /** \brief Export function in a specific language */
  void export_code_body(const std::string& lang,
    std::ostream &stream, const Dict& options) const override;

  /** \brief Serialize */
  void serialize(std::ostream &stream) const override;

  /** \brief Build function from serialization */
  static Function deserialize(std::istream &stream);

  /// With just-in-time compilation using OpenCL
  bool just_in_time_opencl_;

  /// With just-in-time compilation for the sparsity propagation
  bool just_in_time_sparsity_;
};


} // namespace casadi

/// \endcond
#endif // CASADI_SX_FUNCTION_HPP
