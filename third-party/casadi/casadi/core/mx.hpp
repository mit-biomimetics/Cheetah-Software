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


#ifndef CASADI_MX_HPP
#define CASADI_MX_HPP
#include "shared_object.hpp"
#include "matrix.hpp"
#include "generic_expression.hpp"
#include "generic_type.hpp"
#include "printable.hpp"
#include <vector>
namespace casadi {

  /** \brief  Forward declaration */
  class MXNode;
  class Function;

  /** \brief MX - Matrix expression

      The MX class is used to build up trees made up from MXNodes. It is a more general
      graph representation than the scalar expression, SX, and much less efficient for small
      objects. On the other hand, the class allows much more general operations than does SX,
      in particular matrix valued operations and calls to arbitrary differentiable functions.

      The MX class is designed to have identical syntax with the Matrix<> template class,
      and uses DM (i.e. Matrix<double>) as its internal representation of the values at a node. By keeping
      the syntaxes identical, it is possible to switch from one class to the other,
      as well as inlining MX functions to SXElem functions.

      Note that an operation is always "lazy", making a matrix multiplication will create a
      matrix multiplication node, not perform the actual multiplication.

      \author Joel Andersson
      \date 2010-2011
  */
  class CASADI_EXPORT MX :
    public SWIG_IF_ELSE(GenericExpressionCommon, GenericExpression<MX>),
    public SWIG_IF_ELSE(PrintableCommon, Printable<MX>),
    public GenericMatrix<MX>,
    public SharedObject {
  public:
    /** \brief Get type name */
    static std::string type_name() {return "MX";}

    /** \brief  Default constructor */
    MX();

    /** \brief Create a sparse matrix with all structural zeros */
    MX(casadi_int nrow, casadi_int ncol);

#ifndef SWIG
    /** \brief Create a sparse matrix with all structural zeros */
    explicit MX(const std::pair<casadi_int, casadi_int>& rc);
#endif // SWIG

    /** \brief Create a sparse matrix from a sparsity pattern.
        Same as MX::ones(sparsity)
     */
    explicit MX(const Sparsity& sp);

    /** \brief Construct matrix with a given sparsity and nonzeros */
    MX(const Sparsity& sp, const MX& val);

    /** \brief  Create scalar constant (also implicit type conversion) */
    MX(double x);

    /** \brief  Copy constructor */
    MX(const MX& x);

    /** \brief  Create vector constant (also implicit type conversion) */
    MX(const std::vector<double> &x);

    /** \brief  Create sparse matrix constant (also implicit type conversion) */
    MX(const Matrix<double> &x);

/// \cond INTERNAL
    /** \brief  Destructor */
    ~MX();
/// \endcond

#ifndef SWIG
/// \cond INTERNAL
    /** \brief  Create from node */
    static MX create(MXNode* node);

    /** \brief  Create from node (multiple-outputs) */
    static std::vector<MX> createMultipleOutput(MXNode* node);
/// \endcond

    /** \brief Get the sparsity pattern */
    const Sparsity& sparsity() const;

    /// \cond INTERNAL
    /// Scalar type
    typedef MX ScalarType;
    /// \endcond

    /// Base class
    typedef GenericMatrix<MX> B;

    /// Expose base class functions
    using B::horzsplit;
    using B::diagsplit;
    using B::vertsplit;
    using B::mtimes;
    using B::repmat;
#endif // SWIG

    /// Returns the truth value of an MX expression
    bool __nonzero__() const;

    /** \brief Get an owning reference to the sparsity pattern */
    Sparsity get_sparsity() const { return sparsity();}

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase rows and/or columns of a matrix */
    void erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                bool ind1=false);

    /** \brief Erase a submatrix (leaving structural zeros in its place)
        Erase elements of a matrix */
    void erase(const std::vector<casadi_int>& rr, bool ind1=false);

    /** \brief Enlarge matrix
        Make the matrix larger by inserting empty rows and columns, keeping the existing non-zeros */
    void enlarge(casadi_int nrow, casadi_int ncol,
                  const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc,
                  bool ind1=false);

    MX operator-() const;

#ifndef SWIG
    /// \cond INTERNAL
    ///@{
    /** \brief  Access a member of the node */
    MXNode* operator->();

    /** \brief  Const access a member of the node */
    const MXNode* operator->() const;
    ///@}
    /// \endcond
#endif // SWIG

    /** \brief Get the nth dependency as MX */
    MX dep(casadi_int ch=0) const;

    /** \brief  Number of outputs */
    casadi_int n_out() const;

    /** \brief  Get an output */
    MX get_output(casadi_int oind) const;

    /** \brief Get the number of dependencies of a binary SXElem */
    casadi_int n_dep() const;

    /// Get the name.
    std::string name() const;

    /// Get the value (only for scalar constant nodes)
    explicit operator double() const;

    /// Get the value (only for constant nodes)
    explicit operator Matrix<double>() const;

    /// Check if symbolic
    bool is_symbolic() const;

    /// Check if constant
    bool is_constant() const;

    /// Check if evaluation
    bool is_call() const;

    /// Get function - only valid when is_call() is true
    Function which_function() const;

    /// Check if evaluation output
    bool is_output() const;

    /// Get the index of evaluation output - only valid when is_output() is true
    casadi_int which_output() const;

    /// Is it a certain operation
    bool is_op(casadi_int op) const;

    /// Check if multiplication
    bool is_multiplication() const;

    /// Check if commutative operation
    bool is_commutative() const;

    /// Check if norm
    bool is_norm() const;

    /** \brief Check if matrix can be used to define function inputs.
        Valid inputs for MXFunctions are combinations of Reshape, concatenations and SymbolicMX
    */
    bool is_valid_input() const;

    /** \brief Get the number of primitives for MXFunction inputs/outputs
    */
    casadi_int n_primitives() const;

    /** \brief Get primitives */
    std::vector<MX> primitives() const;

    /** \brief Split up an expression along symbolic primitives */
    std::vector<MX> split_primitives(const MX& x) const;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(const std::vector<MX>& v) const;

    /// \cond INTERNAL
    /** \brief Detect duplicate symbolic expressions
        If there are symbolic primitives appearing more than once, the function will return
        true and the names of the duplicate expressions will be passed to casadi_warning.
        Note: Will mark the node using MX::set_temp.
        Make sure to call reset_input() after usage.
    */
    bool has_duplicates() const;

    /** \brief Reset the marker for an input expression */
    void reset_input() const;
  /// \endcond

    /** \brief  check if identity */
    bool is_eye() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool is_zero() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool is_one() const;

    /** \brief  check if zero (note that false negative answers are possible) */
    bool is_minus_one() const;

    /** \brief  Is the expression a transpose? */
    bool is_transpose() const;

    /// Checks if expression does not contain NaN or Inf
    bool is_regular() const;

    /// Is binary operation
    bool is_binary() const;

    /// Is unary operation
    bool is_unary() const;

    /// Get operation type
    casadi_int op() const;

    /** Obtain information about node */
    Dict info() const;

    /// \cond INTERNAL
    /// Get the temporary variable
    casadi_int get_temp() const;

    /// Set the temporary variable
    void set_temp(casadi_int t) const;
    /// \endcond

    ///@{
    /** \brief  Create nodes by their ID */
    static MX binary(casadi_int op, const MX &x, const MX &y);
    static MX unary(casadi_int op, const MX &x);
    ///@}

    ///@{
    /** \brief  create a matrix with all inf */
    static MX inf(const Sparsity& sp);
    static MX inf(casadi_int nrow=1, casadi_int ncol=1);
    static MX inf(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    ///@{
    /** \brief  create a matrix with all nan */
    static MX nan(const Sparsity& sp);
    static MX nan(casadi_int nrow=1, casadi_int ncol=1);
    static MX nan(const std::pair<casadi_int, casadi_int>& rc);
    ///@}

    /** \brief  Identity matrix */
    static MX eye(casadi_int ncol);

#ifndef SWIG
    /// Get a const pointer to the node
    MXNode* get() const;
#endif // SWIG

    ///@{
    /// Get a submatrix, single argument
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Sparsity& sp) const;
    ///@}

    /// Get a submatrix, two arguments
    ///@{
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Slice& rr, const Matrix<casadi_int>& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr, const Slice& cc) const;
    void get(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& rr,
                                            const Matrix<casadi_int>& cc) const;
    ///@}

    ///@{
    /// Set a submatrix, single argument
    void set(const MX& m, bool ind1, const Slice& rr);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr);
    void set(const MX& m, bool ind1, const Sparsity& sp);
    ///@}

    ///@{
    /// Set a submatrix, two arguments
    ///@}
    void set(const MX& m, bool ind1, const Slice& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc);
    void set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc);
    ///@}

    ///@{
    /// Get a set of nonzeros
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Slice& kk) const;
    void get_nz(MX& SWIG_OUTPUT(m), bool ind1, const Matrix<casadi_int>& kk) const;
    ///@}

    ///@{
    /// Set a set of nonzeros
    void set_nz(const MX& m, bool ind1, const Slice& kk);
    void set_nz(const MX& m, bool ind1, const Matrix<casadi_int>& kk);
    ///@}

    ///@{
    /** \brief Computes an einstein dense tensor contraction

        Computes the product:
        C_c = A_a + B_b
          where a b c are index/einstein notation in an encoded form

        For example, an matrix-matrix product may be written as:
        C_ij = A_ik B_kj

        The encoded form uses strictly negative numbers to indicate labels.
        For the above example, we would have:
        a {-1, -3} b {-3, -2} c {-1 -2}

    */
    static MX einstein(const MX& A, const MX& B, const MX& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);

    static MX einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c);
    ///@}

#ifndef SWIG
    /// \cond CLUTTER
    ///@{
    /// Functions called by friend functions defined for GenericExpression
    static bool is_equal(const MX& x, const MX& y, casadi_int depth=0);
    static MX mmin(const MX &x);
    static MX mmax(const MX &x);
    ///@}

    ///@{
    /// Functions called by friend functions defined for SparsityInterface
    static MX horzcat(const std::vector<MX>& x);
    static MX diagcat(const std::vector<MX>& x);
    static MX vertcat(const std::vector<MX>& x);
    static std::vector<MX> horzsplit(const MX& x, const std::vector<casadi_int>& offset);
    static std::vector<MX> diagsplit(const MX& x, const std::vector<casadi_int>& offset1,
                                     const std::vector<casadi_int>& offset2);
    static std::vector<MX> vertsplit(const MX& x, const std::vector<casadi_int>& offset);
    static MX blockcat(const std::vector< std::vector<MX > > &v);
    static MX mtimes(const MX& x, const MX& y);
    static MX mac(const MX& x, const MX& y, const MX& z);
    static MX reshape(const MX& x, casadi_int nrow, casadi_int ncol);
    static MX reshape(const MX& x, const Sparsity& sp);
    static MX kron(const MX& x, const MX& b);
    static MX repmat(const MX& x, casadi_int n, casadi_int m=1);
    ///@}

    ///@{
    /// Functions called by friend functions defined for GenericMatrix
    static MX jacobian(const MX& f, const MX& x, const Dict& opts = Dict());
    static MX hessian(const MX& f, const MX& x);
    static MX hessian(const MX& f, const MX& x, MX& g);
    static std::vector<std::vector<MX> >
    forward(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<std::vector<MX> >
    reverse(const std::vector<MX> &ex,
            const std::vector<MX> &arg,
            const std::vector<std::vector<MX> > &v,
            const Dict& opts = Dict());
    static std::vector<bool> which_depends(const MX &expr, const MX &var,
        casadi_int order=1, bool tr=false);
    static MX substitute(const MX& ex, const MX& v, const MX& vdef);
    static std::vector<MX> substitute(const std::vector<MX> &ex,
                                         const std::vector<MX> &v,
                                         const std::vector<MX> &vdef);
    static void substitute_inplace(const std::vector<MX>& v,
                                  std::vector<MX>& vdef,
                                  std::vector<MX>& ex, bool reverse);
    static MX solve(const MX& A, const MX& b, const std::string& lsolver="qr",
                    const Dict& dict = Dict());
    static MX inv_minor(const MX& A);
    static MX inv_node(const MX& A);
    static MX inv(const MX& A, const std::string& lsolver="qr", const Dict& dict = Dict());
    static MX pinv(const MX& A, const std::string& lsolver="qr",
               const Dict& dict = Dict());
    static MX expm_const(const MX& A, const MX& t);
    static MX expm(const MX& A);
    static casadi_int n_nodes(const MX& x);
    static std::string print_operator(const MX& x, const std::vector<std::string>& args);
    static void shared(std::vector<MX>& ex, std::vector<MX>& v,
                              std::vector<MX>& vdef, const std::string& v_prefix,
                              const std::string& v_suffix);
    static MX if_else(const MX& cond, const MX& if_true,
                      const MX& if_false, bool short_circuit=false);
    static MX conditional(const MX& ind, const std::vector<MX> &x, const MX& x_default,
                          bool short_circuit=false);
    static bool depends_on(const MX& x, const MX& arg);
    static MX simplify(const MX& x);
    static MX dot(const MX& x, const MX& y);
    static MX mrdivide(const MX& a, const MX& b);
    static MX mldivide(const MX& a, const MX& b);
    static MX norm_2(const MX& x);
    static MX norm_fro(const MX& x);
    static MX norm_1(const MX& x);
    static MX norm_inf(const MX& x);
    static MX unite(const MX& A, const MX& B);
    static MX trace(const MX& x);
    static MX diag(const MX& x);
    static MX sum2(const MX& x);
    static MX sum1(const MX& x);
    static MX polyval(const MX& p, const MX& x);
    static MX det(const MX& x);
    static std::vector<MX> symvar(const MX& x);
    static MX nullspace(const MX& A);
    static MX repsum(const MX& x, casadi_int n, casadi_int m=1);
    static MX densify(const MX& x, const MX& val=0);
    static MX _bilin(const MX& A, const MX& x, const MX& y);
    static MX _rank1(const MX& A, const MX& alpha, const MX& x, const MX& y);
    static MX project(const MX& x, const Sparsity& sp, bool intersect=false);
    static MX cumsum(const MX &x, casadi_int axis=-1);
    ///@}

    ///@{
    /// Functions called by friend functions defined for this class
    static MX find(const MX& x);
    static MX graph_substitute(const MX& x, const std::vector<MX> &v,
                               const std::vector<MX> &vdef);
    static std::vector<MX> graph_substitute(const std::vector<MX> &ex,
                                            const std::vector<MX> &expr,
                                            const std::vector<MX> &exprs);
    static MX matrix_expand(const MX& e, const std::vector<MX> &boundary,
                            const Dict& options);
    static std::vector<MX> matrix_expand(const std::vector<MX>& e,
                                         const std::vector<MX>& boundary,
                                         const Dict& options);
    static MX lift(const MX& x, const MX& x_guess);
    static DM evalf(const MX& x);
    ///@}
    /// \endcond

#endif // SWIG

    MX printme(const MX& y) const;

#if !defined(SWIG) || defined(DOXYGEN)
/**
\ingroup expression_tools
@{
*/
    /** \brief Find first nonzero
     * If failed, returns the number of rows
     */
    inline friend MX find(const MX& x) {
      return MX::find(x);
    }

    /** \brief Substitute single expression in graph
     * Substitute variable v with expression vdef in an expression ex, preserving nodes
     */
    inline friend MX graph_substitute(const MX& ex, const std::vector<MX> &v,
                                      const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    /** \brief Substitute multiple expressions in graph
     * Substitute variable var with expression expr in
     * multiple expressions, preserving nodes
     */
    inline friend std::vector<MX>
      graph_substitute(const std::vector<MX> &ex,
                       const std::vector<MX> &v,
                       const std::vector<MX> &vdef) {
      return MX::graph_substitute(ex, v, vdef);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
     */
    inline friend MX
      matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>(),
        const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }

    /** \brief Expand MX graph to SXFunction call
     *
     *  Expand the given expression e, optionally
     *  supplying expressions contained in it at which expansion should stop.
     *
     */
    inline friend std::vector<MX>
      matrix_expand(const std::vector<MX>& e,
                    const std::vector<MX> &boundary = std::vector<MX>(),
                    const Dict& options = Dict()) {
      return MX::matrix_expand(e, boundary, options);
    }

    /** \brief Lift the expression
     * Experimental feature
     *
     */
    inline friend MX lift(const MX& x, const MX& x_guess) {
      return MX::lift(x, x_guess);
    }

    /** \brief Inverse node
     *
     */
    inline friend MX inv_node(const MX& x) {
      return MX::inv_node(x);
    }

    /** \brief Evaluates the expression numerically
    *
    * An error is raised when the expression contains symbols
    */
    inline friend DM evalf(const MX& expr) {
      return MX::evalf(expr);
    }

/** @} */
#endif // SWIG

    /** \brief returns itself, but with an assertion attached
    *
    *  If y does not evaluate to 1, a runtime error is raised
    */
    MX attachAssert(const MX& y, const std::string& fail_message="") const;

    /** \brief Monitor an expression
    * Returns itself, but with the side effect of printing the nonzeros along with a comment
    */
    MX monitor(const std::string& comment) const;

    /// Transpose the matrix
    MX T() const;

    /** \brief Get an IM representation of a GetNonzeros or SetNonzeros node */
    Matrix<casadi_int> mapping() const;

    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void set_max_depth(casadi_int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static casadi_int get_max_depth();

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /** \brief Get function inputs */
    static std::vector<MX> get_input(const Function& f);

    /** \brief Get free variables */
    static std::vector<MX> get_free(const Function& f);

    /// Readability typedef
    typedef std::map<std::string, MX> MXDict;

#ifndef SWIG
    ///@{
    /** \brief Called from MXFunction */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const;
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                    std::vector<std::vector<MX> >& fsens) const;
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                    std::vector<std::vector<MX> >& asens) const;
    ///@}

    /// Construct constant matrix with a given sparsity and values
    MX(const Sparsity& sp, double val, bool dummy);

    // Create matrix symbolic primitive
    static MX _sym(const std::string& name, const Sparsity& sp);

  private:

    /// Create an expression from a node: extra dummy arguments to avoid ambiguity for 0/NULL
    MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4);

    // Depth when checking equalities
    static casadi_int eq_depth_;

#endif // SWIG
  };


  ///@{
  /// Readability typedefs
  typedef std::vector<MX> MXVector;
  typedef std::initializer_list<MX> MXIList;
  typedef std::vector<MXVector> MXVectorVector;
  typedef std::map<std::string, MX> MXDict;
  ///@}

} // namespace casadi

#endif // CASADI_MX_HPP
