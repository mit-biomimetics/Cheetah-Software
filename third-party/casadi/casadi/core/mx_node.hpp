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


#ifndef CASADI_MX_NODE_HPP
#define CASADI_MX_NODE_HPP

#include "mx.hpp"
#include "shared_object_internal.hpp"
#include "sx_elem.hpp"
#include "calculus.hpp"
#include "code_generator.hpp"
#include "linsol.hpp"
#include <vector>
#include <stack>

namespace casadi {
  /** \brief Node class for MX objects
      \author Joel Andersson
      \date 2010
      Internal class.
  */
  class CASADI_EXPORT MXNode : public SharedObjectInternal {
    friend class MX;

  public:
    /// Constructor
    MXNode();

    /** \brief  Destructor */
    ~MXNode() override=0;

    /** \brief Check the truth value of this node
     */
    virtual bool __nonzero__() const;

    /** \brief Check if identically zero */
    virtual bool is_zero() const { return false;}

    /** \brief Check if identically one */
    virtual bool is_one() const { return false;}

    /** \brief Check if identically  minus one */
    virtual bool is_minus_one() const { return false;}

    /** \brief Check if a certain value */
    virtual bool is_value(double val) const { return false;}

    /** \brief Check if identity matrix */
    virtual bool is_eye() const { return false;}

    /** \brief Check if unary operation */
    virtual bool is_unary() const { return false;}

    /** \brief Check if binary operation */
    virtual bool is_binary() const { return false;}

    /** \brief Find out which nodes can be inlined */
    void can_inline(std::map<const MXNode*, casadi_int>& nodeind) const;

    /** \brief Print compact */
    std::string print_compact(std::map<const MXNode*, casadi_int>& nodeind,
                             std::vector<std::string>& intermed) const;

    /** \brief  Print expression */
    virtual std::string disp(const std::vector<std::string>& arg) const = 0;

    /** \brief Add a dependent function */
    virtual void add_dependency(CodeGenerator& g) const {}

    /** \brief Is reference counting needed in codegen? */
    virtual bool has_refcount() const { return false;}

    /** \brief Codegen incref */
    virtual void codegen_incref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Codegen decref */
    virtual void codegen_decref(CodeGenerator& g, std::set<void*>& added) const {}

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res) const;

    /** \brief  Evaluate numerically */
    virtual int eval(const double** arg, double** res, casadi_int* iw, double* w) const;

    /** \brief  Evaluate symbolically (SX) */
    virtual int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const;

    /** \brief Calculate forward mode directional derivatives */
    virtual void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const;

    /** \brief Calculate reverse mode directional derivatives */
    virtual void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const;

    /** \brief  Propagate sparsity forward */
    virtual int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

    /** \brief  Propagate sparsity backwards */
    virtual int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

    /** \brief  Get the name */
    virtual const std::string& name() const;

    /** \brief Get name of public class */
    std::string class_name() const override;

    /** \brief  Print a description */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Check if valid function input */
    virtual bool is_valid_input() const { return false;}

    /** \brief Get the number of symbolic primitives */
    virtual casadi_int n_primitives() const;

    /** \brief Get symbolic primitives */
    virtual void primitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool has_duplicates() const;

    /** \brief Reset the marker for an input expression */
    virtual void reset_input() const;

    /** \brief  Check if evaluation output */
    virtual bool is_output() const {return false;}

    /** \brief  Check if a multiple output node */
    virtual bool has_output() const {return false;}

    /** \brief  Get function output */
    virtual casadi_int which_output() const;

    /** \brief  Get called function */
    virtual const Function& which_function() const;

    /** \brief Get the operation */
    virtual casadi_int op() const = 0;

    /** Obtain information about node */
    virtual Dict info() const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    static bool is_equal(const MXNode* x, const MXNode* y, casadi_int depth);
    virtual bool is_equal(const MXNode* node, casadi_int depth) const { return false;}

    /** \brief Get equality checking depth */
    inline static bool maxDepth() { return MX::get_max_depth();}

    /** \brief Checks if two nodes have the same operation and have
     * equivalent dependencies up to a given depth */
    bool sameOpAndDeps(const MXNode* node, casadi_int depth) const;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    const MX& dep(casadi_int ind=0) const { return dep_.at(ind);}

    /** \brief  Number of dependencies */
    casadi_int n_dep() const;

    /** \brief  Number of outputs */
    virtual casadi_int nout() const { return 1;}

    /** \brief  Get an output */
    virtual MX get_output(casadi_int oind) const;

    /// Get the sparsity
    const Sparsity& sparsity() const { return sparsity_;}

    /// Get the sparsity of output oind
    virtual const Sparsity& sparsity(casadi_int oind) const;

    /// Get shape
    casadi_int numel() const { return sparsity().numel(); }
    casadi_int nnz(casadi_int i=0) const { return sparsity(i).nnz(); }
    casadi_int size1() const { return sparsity().size1(); }
    casadi_int size2() const { return sparsity().size2(); }
    std::pair<casadi_int, casadi_int> size() const { return sparsity().size();}

    // Get IO index
    virtual casadi_int ind() const;

    // Get IO segment
    virtual casadi_int segment() const;

    // Get IO offset
    virtual casadi_int offset() const;

    /// Set the sparsity
    void set_sparsity(const Sparsity& sparsity);

    /** \brief Get required length of arg field */
    virtual size_t sz_arg() const { return n_dep();}

    /** \brief Get required length of res field */
    virtual size_t sz_res() const { return nout();}

    /** \brief Get required length of iw field */
    virtual size_t sz_iw() const { return 0;}

    /** \brief Get required length of w field */
    virtual size_t sz_w() const { return 0;}

    /// Set unary dependency
    void set_dep(const MX& dep);

    /// Set binary dependencies
    void set_dep(const MX& dep1, const MX& dep2);

    /// Set ternary dependencies
    void set_dep(const MX& dep1, const MX& dep2, const MX& dep3);

    /// Set multiple dependencies
    void set_dep(const std::vector<MX>& dep);

    /// Convert scalar to matrix
    inline static MX to_matrix(const MX& x, const Sparsity& sp) {
      if (x.size()==sp.size()) {
        return x;
      } else {
        return MX(sp, x);
      }
    }

    /// Get the value (only for scalar constant nodes)
    virtual double to_double() const;

    /// Get the value (only for constant nodes)
    virtual DM get_DM() const;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual casadi_int n_inplace() const { return 0;}

    /// Get an IM representation of a GetNonzeros or SetNonzeros node
    virtual Matrix<casadi_int> mapping() const;

    /// Create a horizontal concatenation node
    virtual MX get_horzcat(const std::vector<MX>& x) const;

    /// Create a horizontal split node
    virtual std::vector<MX> get_horzsplit(const std::vector<casadi_int>& output_offset) const;

    /// Create a repeated matrix node
    virtual MX get_repmat(casadi_int m, casadi_int n) const;

    /// Create a repeated sum node
    virtual MX get_repsum(casadi_int m, casadi_int n) const;

    /// Create a vertical concatenation node (vectors only)
    virtual MX get_vertcat(const std::vector<MX>& x) const;

    /// Create a vertical split node (vectors only)
    virtual std::vector<MX> get_vertsplit(const std::vector<casadi_int>& output_offset) const;

    /// Create a diagonal concatenation node
    virtual MX get_diagcat(const std::vector<MX>& x) const;

    /// Create a diagonal split node
    virtual std::vector<MX> get_diagsplit(const std::vector<casadi_int>& offset1,
                                         const std::vector<casadi_int>& offset2) const;

    /// Transpose
    virtual MX get_transpose() const;

    /// Reshape
    virtual MX get_reshape(const Sparsity& sp) const;

    /** \brief Matrix multiplication and addition */
    virtual MX get_mac(const MX& y, const MX& z) const;

    /** \brief Einstein product and addition */
    virtual MX get_einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_c, const std::vector<casadi_int>& dim_a,
      const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& c, const std::vector<casadi_int>& a,
      const std::vector<casadi_int>& b) const;

    /** \brief Bilinear form */
    virtual MX get_bilin(const MX& x, const MX& y) const;

    /** \brief Bilinear form */
    virtual MX get_rank1(const MX& alpha, const MX& x, const MX& y) const;

    /** \brief Solve a system of linear equations
    *
    *      For system Ax = b:
    *
    *      A->get_solve(b)
    *
    */
    virtual MX get_solve(const MX& r, bool tr, const Linsol& linear_solver) const;

    /** \brief Get the nonzeros of matrix
    *
    *   a->get_nzref(sp,nz)
    *
    *   returns Matrix(sp,a[nz])
    */
    virtual MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const;

    /** \brief Assign the nonzeros of a matrix to another matrix
    *
    *   a->get_nzassign(b,nz)
    *   returns b with b[nz]=a
    */
    virtual MX get_nzassign(const MX& y, const std::vector<casadi_int>& nz) const;

    /** \brief Add the nonzeros of a matrix to another matrix
    *
    *   a->get_nzadd(b,nz)
    *   returns b with b[nz]+=a
    */
    virtual MX get_nzadd(const MX& y, const std::vector<casadi_int>& nz) const;

    /// Get submatrix reference
    virtual MX get_subref(const Slice& i, const Slice& j) const;

    /// Get submatrix assignment
    virtual MX get_subassign(const MX& y, const Slice& i, const Slice& j) const;

    /// Create set sparse
    virtual MX get_project(const Sparsity& sp) const;

    /// Get a unary operation
    virtual MX get_unary(casadi_int op) const;

    /// Get a binary operation operation
    MX get_binary(casadi_int op, const MX& y) const;

    /// Get a binary operation operation (matrix-matrix)
    virtual MX _get_binary(casadi_int op, const MX& y, bool scX, bool scY) const;

    /// Determinant
    virtual MX get_det() const;

    /// Inverse
    virtual MX get_inv() const;

    /// Inner product
    virtual MX get_dot(const MX& y) const;

    /// Frobenius norm
    virtual MX get_norm_fro() const;

    /// Spectral norm
    virtual MX get_norm_2() const;

    /// Infinity norm
    virtual MX get_norm_inf() const;

    /// 1-norm
    virtual MX get_norm_1() const;

    /// Min
    virtual MX get_mmin() const;

    /// Max
    virtual MX get_mmax() const;

    /// Assertion
    MX get_assert(const MX& y, const std::string& fail_message) const;

    /// Monitor
    MX get_monitor(const std::string& comment) const;

    /// Find
    MX get_find() const;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    mutable casadi_int temp;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;

    /** \brief  The sparsity pattern */
    Sparsity sparsity_;

    /** \brief Propagate sparsities forward through a copy operation */
    static void copy_fwd(const bvec_t* arg, bvec_t* res, casadi_int len);

    /** \brief Propagate sparsities backwards through a copy operation */
    static void copy_rev(bvec_t* arg, bvec_t* res, casadi_int len);
  };

  /// \endcond
} // namespace casadi

#endif // CASADI_MX_NODE_HPP
