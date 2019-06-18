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


#ifndef CASADI_FUNCTION_INTERNAL_HPP
#define CASADI_FUNCTION_INTERNAL_HPP

#include "function.hpp"
#include <set>
#include <stack>
#include "code_generator.hpp"
#include "importer.hpp"
#include "sparse_storage.hpp"
#include "options.hpp"
#include "shared_object_internal.hpp"
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

/// \cond INTERNAL

namespace casadi {
  template<typename T>
  std::vector<std::pair<std::string, T>> zip(const std::vector<std::string>& id,
                                             const std::vector<T>& mat) {
    casadi_assert_dev(id.size()==mat.size());
    std::vector<std::pair<std::string, T>> r(id.size());
    for (casadi_uint i=0; i<r.size(); ++i) r[i] = make_pair(id[i], mat[i]);
    return r;
  }

  /// Combine two dictionaries, giving priority to first one
  Dict CASADI_EXPORT combine(const Dict& first, const Dict& second);

  /** \brief Base class for FunctionInternal and LinsolInternal
    \author Joel Andersson
    \date 2017
  */
  class CASADI_EXPORT ProtoFunction : public SharedObjectInternal {
  public:
    /** \brief Constructor */
    ProtoFunction(const std::string& name);

    /** \brief  Destructor */
    ~ProtoFunction() override = 0;

    /** \brief Construct
        Prepares the function for evaluation
     */
    void construct(const Dict& opts);

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation.
        This method is typically called after setting options but before evaluating.
        If passed to another class (in the constructor), this class should invoke
        this function when initialized. */
    virtual void init(const Dict& opts);

    /** \brief Finalize the object creation
        This function, which visits the class hierarchy in reverse order is run after
        init() has been completed.
    */
    virtual void finalize(const Dict& opts);

    /// Checkout a memory object
    casadi_int checkout() const;

    /// Release a memory object
    void release(casadi_int mem) const;

    /// Memory objects
    void* memory(casadi_int ind) const;

    /** \brief Create memory block */
    virtual void* alloc_mem() const {return nullptr;}

    /** \brief Initalize memory block */
    virtual int init_mem(void* mem) const { return 0;}

    /** \brief Free memory block */
    virtual void free_mem(void *mem) const;

    /** \brief Clear all memory (called from destructor) */
    void clear_mem();

  protected:
    /// Name
    std::string name_;

    /// Verbose printout
    bool verbose_;
  private:
    /// Memory objects
    mutable std::vector<void*> mem_;

    /// Unused memory objects
    mutable std::stack<casadi_int> unused_;

#ifdef CASADI_WITH_THREAD
    /// Mutex for thread safety
    mutable std::mutex mtx_;
#endif // CASADI_WITH_THREAD
  };

  /** \brief Internal class for Function
      \author Joel Andersson
      \date 2010-2015
  */
  class CASADI_EXPORT FunctionInternal : public ProtoFunction {
    friend class Function;
  public:
    /** \brief Constructor */
    FunctionInternal(const std::string& name);

    /** \brief  Destructor */
    ~FunctionInternal() override = 0;

    /** \brief  Obtain solver name from Adaptor */
    virtual std::string getAdaptorSolverName() const { return ""; }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize */
    void init(const Dict& opts) override;

    /** \brief Finalize the object creation */
    void finalize(const Dict& opts) override;

    /** \brief Get a public class instance */
    Function self() const { return shared_from_this<Function>();}

    // Factory
    virtual Function factory(const std::string& name,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const Function::AuxOut& aux,
                             const Dict& opts) const;

    // Get list of dependency functions
    virtual std::vector<std::string> get_function() const;

    // Get a dependency function
    virtual const Function& get_function(const std::string &name) const;

    // Check if a particular dependency exists
    virtual bool has_function(const std::string& fname) const {return false;}

    /** \brief Which variables enter with some order
    * \param[in] s_in Input name
    * \param[in] s_out Output name(s)
    * \param[in] order Only 1 (linear) and 2 (nonlinear) allowed
    * \param[in] tr   Flip the relationship. Return which expressions contain the variables
    */
    virtual std::vector<bool> which_depends(const std::string& s_in,
                                           const std::vector<std::string>& s_out,
                                           casadi_int order, bool tr=false) const;

    ///@{
    /** \brief  Is the class able to propagate seeds through the algorithm? */
    virtual bool has_spfwd() const { return false;}
    virtual bool has_sprev() const { return false;}
    ///@}

    ///@{
    /** \brief  Evaluate numerically */
    int eval_gen(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;
    virtual int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;
    ///@}

    /** \brief  Evaluate with symbolic scalars */
    virtual int eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const;

    /** \brief  Evaluate with symbolic matrices */
    virtual void eval_mx(const MXVector& arg, MXVector& res,
                         bool always_inline, bool never_inline) const;

    /** \brief Evaluate with DM matrices */
    virtual std::vector<DM> eval_dm(const std::vector<DM>& arg) const;

    ///@{
    /** \brief Evaluate a function, overloaded */
    int eval_gen(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
      return eval_sx(arg, res, iw, w, mem);
    }
    int eval_gen(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
      return sp_forward(arg, res, iw, w, mem);
    }
    ///@}

    ///@{
    /** \brief Call a function, overloaded */
    void call_gen(const MXVector& arg, MXVector& res, casadi_int npar,
                  bool always_inline, bool never_inline) const;

    template<typename D>
    void call_gen(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
                  casadi_int npar, bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Call a function, templated */
    template<typename M>
      void call(const std::vector<M>& arg, std::vector<M>& res,
               bool always_inline, bool never_inline) const;

    /** Helper function */
    static bool check_mat(const Sparsity& arg, const Sparsity& inp, casadi_int& npar);

    /** \brief Check if input arguments have correct length and dimensions
     */
    template<typename M>
    void check_arg(const std::vector<M>& arg, casadi_int& npar) const;

    /** \brief Check if output arguments have correct length and dimensions */
    template<typename M>
    void check_res(const std::vector<M>& res, casadi_int& npar) const;

    /** \brief Check if input arguments that needs to be replaced
     */
    template<typename M> bool
    matching_arg(const std::vector<M>& arg, casadi_int& npar) const;

    /** \brief Check if output arguments that needs to be replaced */
    template<typename M> bool
    matching_res(const std::vector<M>& arg, casadi_int& npar) const;

    /** \brief Replace 0-by-0 inputs
     */
    template<typename M> std::vector<M>
    replace_arg(const std::vector<M>& arg, casadi_int npar) const;

    /** \brief Replace 0-by-0 outputs */
    template<typename M> std::vector<M>
    replace_res(const std::vector<M>& res, casadi_int npar) const;

    /** \brief Replace 0-by-0 forward seeds */
    template<typename M> std::vector<std::vector<M>>
    replace_fseed(const std::vector<std::vector<M>>& fseed, casadi_int npar) const;

    /** \brief Replace 0-by-0 reverse seeds */
    template<typename M> std::vector<std::vector<M>>
    replace_aseed(const std::vector<std::vector<M>>& aseed, casadi_int npar) const;

    ///@{
    /** \brief Forward mode AD, virtual functions overloaded in derived classes */
    virtual void call_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& fseed,
                            std::vector<std::vector<MX> >& fsens,
                            bool always_inline, bool never_inline) const;
    virtual void call_forward(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& fseed,
                            std::vector<std::vector<SX> >& fsens,
                            bool always_inline, bool never_inline) const;
    ///@}

    ///@{
    /** \brief Reverse mode, virtual functions overloaded in derived classes */
    virtual void call_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                            const std::vector<std::vector<MX> >& aseed,
                            std::vector<std::vector<MX> >& asens,
                            bool always_inline, bool never_inline) const;
    virtual void call_reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                            const std::vector<std::vector<SX> >& aseed,
                            std::vector<std::vector<SX> >& asens,
                            bool always_inline, bool never_inline) const;
    ///@}

    /** \brief Parallel evaluation */
    std::vector<MX> mapsum_mx(const std::vector<MX > &arg, const std::string& parallelization);

    /** \brief Do the derivative functions need nondifferentiated outputs? */
    virtual bool uses_output() const {return false;}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    Function jacobian() const;
    virtual bool has_jacobian() const { return false;}
    virtual Function get_jacobian(const std::string& name,
                                  const std::vector<std::string>& inames,
                                  const std::vector<std::string>& onames,
                                  const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    Function jac() const;
    virtual bool has_jac() const { return false;}
    virtual Function get_jac(const std::string& name,
                             const std::vector<std::string>& inames,
                             const std::vector<std::string>& onames,
                             const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates forward derivatives
     *    forward(nfwd) returns a cached instance if available,
     *    and calls <tt>Function get_forward(casadi_int nfwd)</tt>
     *    if no cached version is available.
     */
    Function forward(casadi_int nfwd) const;
    virtual bool has_forward(casadi_int nfwd) const { return false;}
    virtual Function get_forward(casadi_int nfwd, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    ///@{
    /** \brief Return function that calculates adjoint derivatives
     *    reverse(nadj) returns a cached instance if available,
     *    and calls <tt>Function get_reverse(casadi_int nadj)</tt>
     *    if no cached version is available.
     */
    Function reverse(casadi_int nadj) const;
    virtual bool has_reverse(casadi_int nadj) const { return false;}
    virtual Function get_reverse(casadi_int nadj, const std::string& name,
                                 const std::vector<std::string>& inames,
                                 const std::vector<std::string>& onames,
                                 const Dict& opts) const;
    ///@}

    /** \brief returns a new function with a selection of inputs/outputs of the original */
    virtual Function slice(const std::string& name, const std::vector<casadi_int>& order_in,
                           const std::vector<casadi_int>& order_out, const Dict& opts) const;

    /** \brief Get oracle */
    virtual const Function& oracle() const;

    /** \brief Can derivatives be calculated in any way? */
    bool has_derivative() const;

    /** \brief  Weighting factor for chosing forward/reverse mode */
    virtual double ad_weight() const;

    /** \brief  Weighting factor for chosing forward/reverse mode,
        sparsity propagation */
    virtual double sp_weight() const;

    /** \brief Get Jacobian sparsity */
    virtual Sparsity get_jacobian_sparsity() const;

    ///@{
    /** \brief Get function input(s) and output(s)  */
    virtual const SX sx_in(casadi_int ind) const;
    virtual const SX sx_out(casadi_int ind) const;
    virtual const std::vector<SX> sx_in() const;
    virtual const std::vector<SX> sx_out() const;
    virtual const MX mx_in(casadi_int ind) const;
    virtual const MX mx_out(casadi_int ind) const;
    virtual const std::vector<MX> mx_in() const;
    virtual const std::vector<MX> mx_out() const;
    ///@}

    /// Get free variables (MX)
    virtual std::vector<MX> free_mx() const;

    /// Get free variables (SX)
    virtual std::vector<SX> free_sx() const;

    /** \brief Does the function have free variables */
    virtual bool has_free() const { return false;}

    /** \brief Extract the functions needed for the Lifted Newton method */
    virtual void generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const;

    /** \brief Get the number of atomic operations */
    virtual casadi_int n_instructions() const;

    /** \brief Get an atomic operation operator index */
    virtual casadi_int instruction_id(casadi_int k) const;

    /** \brief Get the (integer) input arguments of an atomic operation */
    virtual std::vector<casadi_int> instruction_input(casadi_int k) const;

    /** \brief Get the floating point output argument of an atomic operation */
    virtual double instruction_constant(casadi_int k) const;

    /** \brief Get the (integer) output argument of an atomic operation */
    virtual std::vector<casadi_int> instruction_output(casadi_int k) const;

    /** \brief Number of nodes in the algorithm */
    virtual casadi_int n_nodes() const;

    /** *\brief get MX expression associated with instruction */
    virtual MX instruction_MX(casadi_int k) const;

    /** \brief Wrap in an Function instance consisting of only one MX call */
    Function wrap() const;

    /** \brief Get function in cache */
    bool incache(const std::string& fname, Function& f) const;

    /** \brief Save function to cache */
    void tocache(const Function& f) const;

    /** \brief Generate code the function */
    void codegen(CodeGenerator& g, const std::string& fname) const;

    /** \brief Generate meta-information allowing a user to evaluate a generated function */
    void codegen_meta(CodeGenerator& g) const;

    /** \brief Codegen sparsities */
    void codegen_sparsities(CodeGenerator& g) const;

    /** \brief Get name in codegen */
    virtual std::string codegen_name(const CodeGenerator& g) const;

    /** \brief Codegen incref for dependencies */
    virtual void codegen_incref(CodeGenerator& g) const {}

    /** \brief Codegen decref for dependencies */
    virtual void codegen_decref(CodeGenerator& g) const {}

    /** \brief Code generate the function  */
    std::string signature(const std::string& fname) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void codegen_declarations(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    virtual void codegen_body(CodeGenerator& g) const;

    /** \brief Export / Generate C code for the dependency function */
    virtual std::string generate_dependencies(const std::string& fname, const Dict& opts) const;

    /** \brief Is codegen supported? */
    virtual bool has_codegen() const { return false;}

    /** \brief Jit dependencies */
    virtual void jit_dependencies(const std::string& fname) {}

    /** \brief Export function in a specific language */
    virtual void export_code(const std::string& lang,
      std::ostream &stream, const Dict& options) const;

    /** \brief Serialize */
    virtual void serialize(std::ostream &stream) const;

    /** \brief Serialize function header */
    void serialize_header(std::ostream &stream) const;

    /** \brief Build function from serialization */
    static void deserialize_header(std::istream& stream,
        std::string& name,
        std::vector<Sparsity>& sp_in, std::vector<Sparsity>& sp_out,
        std::vector<std::string>& names_in, std::vector<std::string>& names_out,
        casadi_int& sz_w, casadi_int& sz_iw);

    /** \brief Display object */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Print more */
    virtual void disp_more(std::ostream& stream) const {}

    /** \brief C-style formatted printing during evaluation */
    void print(const char* fmt, ...) const;

    /** \brief C-style formatted printing to string */
    void sprint(char* buf, size_t buf_sz, const char* fmt, ...) const;

    /** \brief Get function signature: name:(inputs)->(outputs) */
    std::string definition() const;

    /** \brief Print dimensions of inputs and outputs */
    void print_dimensions(std::ostream &stream) const;

    /** \brief Print list of options */
    void print_options(std::ostream &stream) const;

    /** \brief Print all information there is to know about a certain option */
    void print_option(const std::string &name, std::ostream &stream) const;

    /** \brief Print free variables */
    virtual std::vector<std::string> get_free() const;

    /** \brief Get the unidirectional or bidirectional partition */
    void get_partition(casadi_int iind, casadi_int oind, Sparsity& D1, Sparsity& D2,
                      bool compact, bool symmetric,
                      bool allow_forward, bool allow_reverse) const;

    ///@{
    /** \brief Number of input/output nonzeros */
    casadi_int nnz_in() const;
    casadi_int nnz_in(casadi_int ind) const { return sparsity_in(ind).nnz(); }
    casadi_int nnz_out() const;
    casadi_int nnz_out(casadi_int ind) const { return sparsity_out(ind).nnz(); }
    ///@}

    ///@{
    /** \brief Number of input/output elements */
    casadi_int numel_in() const;
    casadi_int numel_in(casadi_int ind) const { return sparsity_in(ind).numel(); }
    casadi_int numel_out(casadi_int ind) const { return sparsity_out(ind).numel(); }
    casadi_int numel_out() const;
    ///@}

    ///@{
    /** \brief Input/output dimensions */
    casadi_int size1_in(casadi_int ind) const { return sparsity_in(ind).size1(); }
    casadi_int size2_in(casadi_int ind) const { return sparsity_in(ind).size2(); }
    casadi_int size1_out(casadi_int ind) const { return sparsity_out(ind).size1(); }
    casadi_int size2_out(casadi_int ind) const { return sparsity_out(ind).size2(); }
    std::pair<casadi_int, casadi_int> size_in(casadi_int ind) const {
      return sparsity_in(ind).size();
    }
    std::pair<casadi_int, casadi_int> size_out(casadi_int ind) const {
      return sparsity_out(ind).size();
    }
    ///@}

    ///@{
    /** \brief Input/output sparsity */
    const Sparsity& sparsity_in(casadi_int ind) const { return sparsity_in_.at(ind); }
    const Sparsity& sparsity_out(casadi_int ind) const { return sparsity_out_.at(ind); }
    ///@}

    ///@{
    /** \brief Are all inputs and outputs scalar */
    bool all_scalar() const;

    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(casadi_int iind, casadi_int oind, bool symmetric) const;

    /// Get the sparsity pattern, forward mode
    template<bool fwd>
    Sparsity getJacSparsityGen(casadi_int iind, casadi_int oind, bool symmetric,
                                casadi_int gr_i=1, casadi_int gr_o=1) const;

    /// A flavor of getJacSparsity that does hierarchical block structure recognition
    Sparsity getJacSparsityHierarchical(casadi_int iind, casadi_int oind) const;

    /** A flavor of getJacSparsity that does hierarchical block
    * structure recognition for symmetric Jacobians
    */
    Sparsity getJacSparsityHierarchicalSymm(casadi_int iind, casadi_int oind) const;

    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& sparsity_jac(casadi_int iind, casadi_int oind, bool compact, bool symmetric) const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolic_output(const std::vector<MX>& arg) const;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in();
    virtual size_t get_n_out();
    ///@}


    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(casadi_int i);
    virtual std::string get_name_out(casadi_int i);
    ///@}

    /** \brief Get default input value */
    virtual double get_default_in(casadi_int ind) const {
      return 0;
    }

    /** \brief Get largest input value */
    virtual double get_max_in(casadi_int ind) const {
      return inf;
    }

    /** \brief Get smallest input value */
    virtual double get_min_in(casadi_int ind) const {
      return -inf;
    }

    /** \brief Get relative tolerance */
    virtual double get_reltol() const {
      return eps;
    }

    /** \brief Get absolute tolerance */
    virtual double get_abstol() const {
      return eps;
    }

    /** \brief Get sparsity of a given input */
    virtual Sparsity get_sparsity_in(casadi_int i);

    /** \brief Get sparsity of a given output */
    virtual Sparsity get_sparsity_out(casadi_int i);

    /** \brief Get input scheme index by name */
    casadi_int index_in(const std::string &name) const {
      for (casadi_int i=0; i<name_in_.size(); ++i) {
        if (name_in_[i]==name) return i;
      }
      casadi_error("FunctionInternal::index_in: could not find entry \""
                   + name + "\". Available names are: " + str(name_in_) + ".");
      return -1;
    }

    /** \brief Get output scheme index by name */
    casadi_int index_out(const std::string &name) const {
      for (casadi_int i=0; i<name_out_.size(); ++i) {
        if (name_out_[i]==name) return i;
      }
      casadi_error("FunctionInternal::index_out: could not find entry \""
                   + name + "\". Available names are: " + str(name_out_) + ".");
      return -1;
    }

    /** \brief  Propagate sparsity forward */
    virtual int sp_forward(const bvec_t** arg, bvec_t** res,
                            casadi_int* iw, bvec_t* w, void* mem) const;

    /** \brief  Propagate sparsity backwards */
    virtual int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const;

    /** \brief Get number of temporary variables needed */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

    /** \brief Get required length of arg field */
    size_t sz_arg() const { return sz_arg_per_ + sz_arg_tmp_;}

    /** \brief Get required length of res field */
    size_t sz_res() const { return sz_res_per_ + sz_res_tmp_;}

    /** \brief Get required length of iw field */
    size_t sz_iw() const { return sz_iw_per_ + sz_iw_tmp_;}

    /** \brief Get required length of w field */
    size_t sz_w() const { return sz_w_per_ + sz_w_tmp_;}

    /** \brief Ensure required length of arg field */
    void alloc_arg(size_t sz_arg, bool persistent=false);

    /** \brief Ensure required length of res field */
    void alloc_res(size_t sz_res, bool persistent=false);

    /** \brief Ensure required length of iw field */
    void alloc_iw(size_t sz_iw, bool persistent=false);

    /** \brief Ensure required length of w field */
    void alloc_w(size_t sz_w, bool persistent=false);

    /** \brief Ensure work vectors long enough to evaluate function */
    void alloc(const Function& f, bool persistent=false);

    /// Get all statistics
    virtual Dict get_stats(void* mem) const { return Dict();}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const {}

    /** \brief Set the (temporary) work vectors */
    virtual void set_temp(void* mem, const double** arg, double** res,
                          casadi_int* iw, double* w) const {}

    /** \brief Set the (persistent and temporary) work vectors */
    void setup(void* mem, const double** arg, double** res, casadi_int* iw, double* w) const;

    ///@{
    /** \brief Calculate derivatives by multiplying the full Jacobian and multiplying */
    virtual bool fwdViaJac(casadi_int nfwd) const;
    virtual bool adjViaJac(casadi_int nadj) const;
    ///@}

    /** Obtain information about function */
    virtual Dict info() const;

    /** \brief Generate/retrieve cached serial map */
    Function map(casadi_int n, const std::string& parallelization) const;

    /// Number of inputs and outputs
    size_t n_in_, n_out_;

    /// Input and output sparsity
    std::vector<Sparsity> sparsity_in_, sparsity_out_;

    /// Input and output scheme
    std::vector<std::string> name_in_, name_out_;

    /** \brief  Use just-in-time compiler */
    bool jit_;

    /** \brief Numerical evaluation redirected to a C function */
    eval_t eval_;

    /** \brief Dict of statistics (resulting from evaluate) */
    Dict stats_;

    /** \brief Reference counting in codegen? */
    bool has_refcount_;

    /// Function cache
    mutable std::map<std::string, WeakRef> cache_;

    /// Cache for full Jacobian
    mutable WeakRef jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    mutable SparseStorage<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// If the function is the derivative of another function
    Function derivative_of_;

    /// User-set field
    void* user_data_;

    /// Just-in-time compiler
    std::string compilerplugin_;
    Importer compiler_;
    Dict jit_options_;

    /// Penalty factor for using a complete Jacobian to calculate directional derivatives
    double jac_penalty_;

    // Types of derivative calculation permitted
    bool enable_forward_, enable_reverse_, enable_jacobian_, enable_fd_;

    /// Weighting factor for derivative calculation and sparsity pattern calculation
    double ad_weight_, ad_weight_sp_;

    /// Maximum number of sensitivity directions
    casadi_int max_num_dir_;

    /// Errors are thrown when NaN is produced
    bool regularity_check_;

    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;

    // Print timing statistics
    bool print_time_;

    // Finite difference step
    Dict fd_options_;

    // Finite difference step size
    double fd_step_;

    // Finite difference method
    std::string fd_method_;

    /** \brief Check if the function is of a particular type */
    virtual bool is_a(const std::string& type, bool recursive) const;

    /** \brief Can a derivative direction be skipped */
    template<typename MatType>
    static bool purgable(const std::vector<MatType>& seed);

    /** \brief Symbolic expressions for the forward seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    fwd_seed(casadi_int nfwd) const;

    /** \brief Symbolic expressions for the adjoint seeds */
    template<typename MatType>
    std::vector<std::vector<MatType> >
    symbolicAdjSeed(casadi_int nadj, const std::vector<MatType>& v) const;

  protected:
    /** \brief Populate jac_sparsity_ and jac_sparsity_compact_ during initialization */
    void set_jac_sparsity(const Sparsity& sp);

  private:
    /** \brief Memory that is persistent during a call (but not between calls) */
    size_t sz_arg_per_, sz_res_per_, sz_iw_per_, sz_w_per_;

    /** \brief Temporary memory inside a function */
    size_t sz_arg_tmp_, sz_res_tmp_, sz_iw_tmp_, sz_w_tmp_;

    /** \brief Fall back to eval_DM */
    int eval_fallback(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;
  };

  // Template implementations
  template<typename MatType>
  bool FunctionInternal::purgable(const std::vector<MatType>& v) {
    for (auto i=v.begin(); i!=v.end(); ++i) {
      if (!i->is_zero()) return false;
    }
    return true;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::
  fwd_seed(casadi_int nfwd) const {
    std::vector<std::vector<MatType>> fseed(nfwd);
    for (casadi_int dir=0; dir<nfwd; ++dir) {
      fseed[dir].resize(n_in_);
      for (casadi_int iind=0; iind<n_in_; ++iind) {
        std::string n = "f" + str(dir) + "_" +  name_in_[iind];
        fseed[dir][iind] = MatType::sym(n, sparsity_in(iind));
      }
    }
    return fseed;
  }

  template<typename MatType>
  std::vector<std::vector<MatType> >
  FunctionInternal::
  symbolicAdjSeed(casadi_int nadj, const std::vector<MatType>& v) const {
    std::vector<std::vector<MatType> > aseed(nadj, v);
    for (casadi_int dir=0; dir<nadj; ++dir) {
      // Replace symbolic inputs
      casadi_int oind=0;
      for (typename std::vector<MatType>::iterator i=aseed[dir].begin();
          i!=aseed[dir].end();
          ++i, ++oind) {
        // Name of the adjoint seed
        std::stringstream ss;
        ss << "a";
        if (nadj>1) ss << dir << "_";
        ss << oind;

        // Save to matrix
        *i = MatType::sym(ss.str(), i->sparsity());

      }
    }
    return aseed;
  }

  template<typename M>
  void FunctionInternal::call(const std::vector<M>& arg, std::vector<M>& res,
                              bool always_inline, bool never_inline) const {
    // If all inputs are scalar ...
    if (all_scalar()) {
      // ... and some arguments are matrix-valued with matching dimensions ...
      bool matrix_call = false;
      std::pair<casadi_int, casadi_int> sz;
      for (auto&& a : arg) {
        if (!a.is_scalar() && !a.is_empty()) {
          if (!matrix_call) {
            // Matrix call
            matrix_call = true;
            sz = a.size();
          } else if (a.size()!=sz) {
            // Not same dimensions
            matrix_call = false;
            break;
          }
        }
      }

      // ... then, call multiple times
      if (matrix_call) {
        // Start with zeros
        res.resize(n_out_);
        M z = M::zeros(sz);
        for (auto&& a : res) a = z;
        // Call multiple times
        std::vector<M> arg1 = arg, res1;
        for (casadi_int c=0; c<sz.second; ++c) {
          for (casadi_int r=0; r<sz.first; ++r) {
            // Get scalar arguments
            for (casadi_int i=0; i<arg.size(); ++i) {
              if (arg[i].size()==sz) arg1[i] = arg[i](r, c);
            }
            // Call recursively with scalar arguments
            call(arg1, res1, always_inline, never_inline);
            // Get results
            casadi_assert_dev(res.size() == res1.size());
            for (casadi_int i=0; i<res.size(); ++i) res[i](r, c) = res1[i];
          }
        }
        // All elements assigned
        return;
      }
    }

    // Check if inputs need to be replaced
    casadi_int npar = 1;
    if (!matching_arg(arg, npar)) {
      return call(replace_arg(arg, npar), res, always_inline, never_inline);
    }

    // Call the type-specific method
    call_gen(arg, res, npar, always_inline, never_inline);
  }

  template<typename D>
  void FunctionInternal::
  call_gen(const std::vector<Matrix<D> >& arg, std::vector<Matrix<D> >& res,
           casadi_int npar, bool always_inline, bool never_inline) const {
    casadi_assert(!never_inline, "Call-nodes only possible in MX expressions");
    casadi_assert_dev(arg.size()==n_in_);

    // Which arguments require mapped evaluation
    std::vector<bool> mapped(n_in_);
    for (casadi_int i=0; i<n_in_; ++i) {
      mapped[i] = arg[i].size2()!=size2_in(i);
    }

    // Check if matching input sparsity
    std::vector<bool> matching(n_in_);
    bool any_mismatch = false;
    for (casadi_int i=0; i<n_in_; ++i) {
      if (mapped[i]) {
        matching[i] = arg[i].sparsity().is_stacked(sparsity_in(i), npar);
      } else {
        matching[i] = arg[i].sparsity()==sparsity_in(i);
      }
      any_mismatch = any_mismatch || !matching[i];
    }

    // Correct input sparsity via recursive call if needed
    if (any_mismatch) {
      std::vector<Matrix<D> > arg2(arg);
      for (casadi_int i=0; i<n_in_; ++i) {
        if (!matching[i]) {
          if (mapped[i]) {
            arg2[i] = project(arg2[i], repmat(sparsity_in(i), 1, npar));
          } else {
            arg2[i] = project(arg2[i], sparsity_in(i));
          }
        }
      }
      return call_gen(arg2, res, npar, always_inline, never_inline);
    }

    // Allocate results
    res.resize(n_out_);
    for (casadi_int i=0; i<n_out_; ++i) {
      if (!res[i].sparsity().is_stacked(sparsity_out(i), npar)) {
        res[i] = Matrix<D>::zeros(repmat(sparsity_out(i), 1, npar));
      }
    }

    // Allocate temporary memory if needed
    std::vector<casadi_int> iw_tmp(sz_iw());
    std::vector<D> w_tmp(sz_w());

    // Get pointers to input arguments
    std::vector<const D*> argp(sz_arg());
    for (casadi_int i=0; i<n_in_; ++i) argp[i]=get_ptr(arg[i]);

    // Get pointers to output arguments
    std::vector<D*> resp(sz_res());
    for (casadi_int i=0; i<n_out_; ++i) resp[i]=get_ptr(res[i]);

    // For all parallel calls
    for (casadi_int p=0; p<npar; ++p) {
      // Call memory-less
      if (eval_gen(get_ptr(argp), get_ptr(resp),
                   get_ptr(iw_tmp), get_ptr(w_tmp), memory(0))) {
        casadi_error("Evaluation failed");
      }
      // Update offsets
      if (p==npar-1) break;
      for (casadi_int i=0; i<n_in_; ++i) if (mapped[i]) argp[i] += nnz_in(i);
      for (casadi_int i=0; i<n_out_; ++i) resp[i] += nnz_out(i);
    }
  }

  template<typename M>
  void FunctionInternal::check_arg(const std::vector<M>& arg, casadi_int& npar) const {
    casadi_assert(arg.size()==n_in_, "Incorrect number of inputs: Expected "
                          + str(n_in_) + ", got " + str(arg.size()));
    for (casadi_int i=0; i<n_in_; ++i) {
      if (!check_mat(arg[i].sparsity(), sparsity_in(i), npar)) {
        // Dimensions
        std::string d_arg = str(arg[i].size1()) + "-by-" + str(arg[i].size2());
        std::string d_in = str(size1_in(i)) + "-by-" + str(size2_in(i));
        casadi_error("Input " + str(i) + " (" + name_in_[i] + ") has mismatching shape. "
                     "Got " + d_arg + ". Allowed dimensions, in general, are:\n"
                     " - The input dimension N-by-M (here " + d_in + ")\n"
                     " - A scalar, i.e. 1-by-1\n"
                     " - M-by-N if N=1 or M=1 (i.e. a transposed vector)\n"
                     " - N-by-M1 if K*M1=M for some K (argument repeated horizontally)\n"
                     " - N-by-P*M, indicating evaluation with multiple arguments (P must be a "
                     "multiple of " + str(npar) + " for consistency with previous inputs)");
      }
    }
  }

  template<typename M>
  void FunctionInternal::check_res(const std::vector<M>& res, casadi_int& npar) const {
    casadi_assert(res.size()==n_out_, "Incorrect number of outputs: Expected "
                          + str(n_out_) + ", got " + str(res.size()));
    for (casadi_int i=0; i<n_out_; ++i) {
      casadi_assert(check_mat(res[i].sparsity(), sparsity_out(i), npar),
                    "Output " + str(i) + " (" + name_out_[i] + ") has mismatching shape. "
                    "Expected " + str(size_out(i)) + ", got " + str(res[i].size()));
    }
  }

  template<typename M>
  bool FunctionInternal::matching_arg(const std::vector<M>& arg, casadi_int& npar) const {
    check_arg(arg, npar);
    for (casadi_int i=0; i<n_in_; ++i) {
      if (arg.at(i).size1()!=size1_in(i)) return false;
      if (arg.at(i).size2()!=size2_in(i) && arg.at(i).size2()!=npar*size2_in(i)) return false;
    }
    return true;
  }

  template<typename M>
  bool FunctionInternal::matching_res(const std::vector<M>& res, casadi_int& npar) const {
    check_res(res, npar);
    for (casadi_int i=0; i<n_out_; ++i) {
      if (res.at(i).size1()!=size1_out(i)) return false;
      if (res.at(i).size2()!=size2_out(i) && res.at(i).size2()!=npar*size2_out(i)) return false;
    }
    return true;
  }

  template<typename M>
  M replace_mat(const M& arg, const Sparsity& inp, casadi_int npar) {
    if (arg.size()==inp.size()) {
      // Matching dimensions already
      return arg;
    } else if (arg.is_empty()) {
      // Empty matrix means set zero
      return M(inp.size());
    } else if (arg.is_scalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else if (arg.is_vector() && inp.size()==std::make_pair(arg.size2(), arg.size1())) {
      // Transpose vector
      return arg.T();
    } else if (arg.size1()==inp.size1() && arg.size2()>0 && inp.size2()>0
               && inp.size2()%arg.size2()==0) {
      // Horizontal repmat
      return repmat(arg, 1, inp.size2()/arg.size2());
    } else {
      // Multiple evaluation
      return repmat(arg, 1, (npar*inp.size2())/arg.size2());
    }
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  replace_arg(const std::vector<M>& arg, casadi_int npar) const {
    std::vector<M> r(arg.size());
    for (casadi_int i=0; i<r.size(); ++i) r[i] = replace_mat(arg[i], sparsity_in(i), npar);
    return r;
  }

  template<typename M>
  std::vector<M> FunctionInternal::
  replace_res(const std::vector<M>& res, casadi_int npar) const {
    std::vector<M> r(res.size());
    for (casadi_int i=0; i<r.size(); ++i) r[i] = replace_mat(res[i], sparsity_out(i), npar);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> > FunctionInternal::
  replace_fseed(const std::vector<std::vector<M> >& fseed, casadi_int npar) const {
    std::vector<std::vector<M> > r(fseed.size());
    for (casadi_int d=0; d<r.size(); ++d) r[d] = replace_arg(fseed[d], npar);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> > FunctionInternal::
  replace_aseed(const std::vector<std::vector<M> >& aseed, casadi_int npar) const {
    std::vector<std::vector<M> > r(aseed.size());
    for (casadi_int d=0; d<r.size(); ++d) r[d] = replace_res(aseed[d], npar);
    return r;
  }

} // namespace casadi

/// \endcond

#endif // CASADI_FUNCTION_INTERNAL_HPP
