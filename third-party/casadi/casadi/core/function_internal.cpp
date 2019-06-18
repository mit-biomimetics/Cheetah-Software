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


#include "function_internal.hpp"
#include "casadi_call.hpp"
#include "casadi_misc.hpp"
#include "global_options.hpp"
#include "external.hpp"
#include "finite_differences.hpp"
#include "map.hpp"

#include <typeinfo>
#include <cctype>
#ifdef WITH_DL
#include <cstdlib>
#include <ctime>
#endif // WITH_DL
#include <iomanip>

using namespace std;

namespace casadi {
  Dict combine(const Dict& first, const Dict& second) {
    if (first.empty()) return second;
    if (second.empty()) return first;

    Dict ret = second;
    for (auto&& op : first) {
      ret[op.first] = op.second;
    }

    return ret;
  }

  ProtoFunction::ProtoFunction(const std::string& name) : name_(name) {
    // Default options (can be overridden in derived classes)
    verbose_ = false;
  }

  FunctionInternal::FunctionInternal(const std::string& name) : ProtoFunction(name) {
    // Make sure valid function name
    if (!Function::check_name(name_)) {
      casadi_error("Function name is not valid. A valid function name is a string "
                   "starting with a letter followed by letters, numbers or "
                   "non-consecutive underscores. It may also not match the keywords "
                   "'null', 'jac' or 'hess'. Got '" + name_ + "'");
    }

    // By default, reverse mode is about twice as expensive as forward mode
    ad_weight_ = 0.33; // i.e. nf <= 2*na <=> 1/3*nf <= (1-1/3)*na, forward when tie
    // Both modes equally expensive by default (no "taping" needed)
    ad_weight_sp_ = 0.49; // Forward when tie
    jac_penalty_ = 2;
    max_num_dir_ = GlobalOptions::getMaxNumDir();
    user_data_ = nullptr;
    regularity_check_ = false;
    inputs_check_ = true;
    jit_ = false;
    compilerplugin_ = "clang";
    print_time_ = true;
    eval_ = nullptr;
    has_refcount_ = false;
    enable_forward_ = true;
    enable_reverse_ = true;
    enable_jacobian_ = true;
    enable_fd_ = false;
    sz_arg_tmp_ = 0;
    sz_res_tmp_ = 0;
    sz_iw_tmp_ = 0;
    sz_w_tmp_ = 0;
    sz_arg_per_ = 0;
    sz_res_per_ = 0;
    sz_iw_per_ = 0;
    sz_w_per_ = 0;
  }

  ProtoFunction::~ProtoFunction() {
    for (void* m : mem_) {
      if (m!=nullptr) casadi_warning("Memory object has not been properly freed");
    }
    mem_.clear();
  }

  FunctionInternal::~FunctionInternal() {
  }

  void ProtoFunction::construct(const Dict& opts) {
    // Sanitize dictionary is needed
    if (!Options::is_sane(opts)) {
      // Call recursively
      return construct(Options::sanitize(opts));
    }

    // Make sure all options exist
    get_options().check(opts);

    // Initialize the class hierarchy
    try {
      init(opts);
    } catch (exception& e) {
      casadi_error("Error calling " + class_name() + "::init for '" + name_ + "':\n"
        + string(e.what()));
    }

    // Revisit class hierarchy in reverse order
    try {
      finalize(opts);
    } catch (exception& e) {
      casadi_error("Error calling " + class_name() + "::finalize for '" + name_ + "':\n"
        + string(e.what()));
    }
  }

  Options ProtoFunction::options_
  = {{},
     {{"verbose",
       {OT_BOOL,
        "Verbose evaluation -- for debugging"}}
      }
  };

  Options FunctionInternal::options_
  = {{&ProtoFunction::options_},
      {{"ad_weight",
       {OT_DOUBLE,
        "Weighting factor for derivative calculation."
        "When there is an option of either using forward or reverse mode "
        "directional derivatives, the condition ad_weight*nf<=(1-ad_weight)*na "
        "is used where nf and na are estimates of the number of forward/reverse "
        "mode directional derivatives needed. By default, ad_weight is calculated "
        "automatically, but this can be overridden by setting this option. "
        "In particular, 0 means forcing forward mode and 1 forcing reverse mode. "
        "Leave unset for (class specific) heuristics."}},
      {"ad_weight_sp",
       {OT_DOUBLE,
        "Weighting factor for sparsity pattern calculation calculation."
        "Overrides default behavior. Set to 0 and 1 to force forward and "
        "reverse mode respectively. Cf. option \"ad_weight\"."}},
      {"jac_penalty",
       {OT_DOUBLE,
        "When requested for a number of forward/reverse directions,   "
        "it may be cheaper to compute first the full jacobian and then "
        "multiply with seeds, rather than obtain the requested directions "
        "in a straightforward manner. "
        "Casadi uses a heuristic to decide which is cheaper. "
        "A high value of 'jac_penalty' makes it less likely for the heurstic "
        "to chose the full Jacobian strategy. "
        "The special value -1 indicates never to use the full Jacobian strategy"}},
      {"user_data",
       {OT_VOIDPTR,
        "A user-defined field that can be used to identify "
        "the function or pass additional information"}},
      {"regularity_check",
       {OT_BOOL,
        "Throw exceptions when NaN or Inf appears during evaluation"}},
      {"inputs_check",
       {OT_BOOL,
        "Throw exceptions when the numerical values of the inputs don't make sense"}},
      {"gather_stats",
       {OT_BOOL,
        "Deprecated option (ignored): Statistics are now always collected."}},
      {"input_scheme",
       {OT_STRINGVECTOR,
        "Deprecated option (ignored)"}},
      {"output_scheme",
       {OT_STRINGVECTOR,
        "Deprecated option (ignored)"}},
      {"jit",
       {OT_BOOL,
        "Use just-in-time compiler to speed up the evaluation"}},
      {"compiler",
       {OT_STRING,
        "Just-in-time compiler plugin to be used."}},
      {"jit_options",
       {OT_DICT,
        "Options to be passed to the jit compiler."}},
      {"derivative_of",
       {OT_FUNCTION,
        "The function is a derivative of another function. "
        "The type of derivative (directional derivative, Jacobian) "
        "is inferred from the function name."}},
      {"max_num_dir",
       {OT_INT,
        "Specify the maximum number of directions for derivative functions."
        " Overrules the builtin optimized_num_dir."}},
      {"print_time",
       {OT_BOOL,
        "print information about execution time"}},
      {"enable_forward",
       {OT_BOOL,
        "Enable derivative calculation using generated functions for"
        " Jacobian-times-vector products - typically using forward mode AD"
        " - if available. [default: true]"}},
      {"enable_reverse",
        {OT_BOOL,
        "Enable derivative calculation using generated functions for"
        " transposed Jacobian-times-vector products - typically using reverse mode AD"
        " - if available. [default: true]"}},
      {"enable_jacobian",
        {OT_BOOL,
        "Enable derivative calculation using generated functions for"
        " Jacobians of all differentiable outputs with respect to all differentiable inputs"
        " - if available. [default: true]"}},
      {"enable_fd",
       {OT_BOOL,
        "Enable derivative calculation by finite differencing. [default: false]]"}},
      {"fd_options",
       {OT_DICT,
        "Options to be passed to the finite difference instance"}},
      {"fd_method",
       {OT_STRING,
        "Method for finite differencing [default 'central']"}}
     }
  };

  void ProtoFunction::init(const Dict& opts) {
    // Read options
    for (auto&& op : opts) {
      if (op.first=="verbose") {
        verbose_ = op.second;
      }
    }
  }

  void FunctionInternal::init(const Dict& opts) {
    // Call the initialization method of the base class
    ProtoFunction::init(opts);

    // Default options
    fd_step_ = 1e-8;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="jac_penalty") {
        jac_penalty_ = op.second;
      } else if (op.first=="user_data") {
        user_data_ = op.second.to_void_pointer();
      } else if (op.first=="regularity_check") {
        regularity_check_ = op.second;
      } else if (op.first=="inputs_check") {
        inputs_check_ = op.second;
      } else if (op.first=="gather_stats") {
        casadi_warning("Deprecated option \"gather_stats\": Always enabled");
      } else if (op.first=="input_scheme") {
        casadi_warning("Deprecated option: \"input_scheme\" set via constructor");
      } else if (op.first=="output_scheme") {
        casadi_warning("Deprecated option: \"output_scheme\" set via constructor");
      } else if (op.first=="jit") {
        jit_ = op.second;
      } else if (op.first=="compiler") {
        compilerplugin_ = op.second.to_string();
      } else if (op.first=="jit_options") {
        jit_options_ = op.second;
      } else if (op.first=="derivative_of") {
        derivative_of_ = op.second;
      } else if (op.first=="ad_weight") {
        ad_weight_ = op.second;
      } else if (op.first=="ad_weight_sp") {
        ad_weight_sp_ = op.second;
      } else if (op.first=="max_num_dir") {
        max_num_dir_ = op.second;
      } else if (op.first=="print_time") {
        print_time_ = op.second;
      } else if (op.first=="enable_forward") {
        enable_forward_ = op.second;
      } else if (op.first=="enable_reverse") {
        enable_reverse_ = op.second;
      } else if (op.first=="enable_jacobian") {
        enable_jacobian_ = op.second;
      } else if (op.first=="enable_fd") {
        enable_fd_ = op.second;
      } else if (op.first=="fd_options") {
        fd_options_ = op.second;
      } else if (op.first=="fd_method") {
        fd_method_ = op.second.to_string();
      }
    }

    // Verbose?
    if (verbose_) casadi_message(name_ + "::init");

    // Get the number of inputs
    n_in_ = get_n_in();
    if (n_in_>=10000) {
      casadi_warning("Function " + name_ + " has many inputs (" + str(n_in_) + "). "
                     "Changing the problem formulation is strongly encouraged.");
    }

    // Get the number of outputs
    n_out_ = get_n_out();
    if (n_out_>=10000) {
      casadi_warning("Function " + name_ + " has many outputs (" + str(n_out_) + "). "
                     "Changing the problem formulation is strongly encouraged.");
    }

    // Query input sparsities if not already provided
    if (sparsity_in_.empty()) {
      sparsity_in_.resize(n_in_);
      for (casadi_int i=0; i<n_in_; ++i) sparsity_in_[i] = get_sparsity_in(i);
    } else {
      casadi_assert_dev(sparsity_in_.size()==n_in_);
    }

    // Query output sparsities if not already provided
    if (sparsity_out_.empty()) {
      sparsity_out_.resize(n_out_);
      for (casadi_int i=0; i<n_out_; ++i) sparsity_out_[i] = get_sparsity_out(i);
    } else {
      casadi_assert_dev(sparsity_out_.size()==n_out_);
    }

    // Query input names if not already provided
    if (name_in_.empty()) {
      name_in_.resize(n_in_);
      for (casadi_int i=0; i<n_in_; ++i) name_in_[i] = get_name_in(i);
    } else {
      casadi_assert_dev(name_in_.size()==n_in_);
    }

    // Query output names if not already provided
    if (name_out_.empty()) {
      name_out_.resize(n_out_);
      for (casadi_int i=0; i<n_out_; ++i) name_out_[i] = get_name_out(i);
    } else {
      casadi_assert_dev(name_out_.size()==n_out_);
    }

    // Allocate memory for function inputs and outputs
    sz_arg_per_ += n_in_;
    sz_res_per_ += n_out_;

    // Resize the matrix that holds the sparsity of the Jacobian blocks
    jac_sparsity_ = jac_sparsity_compact_ = SparseStorage<Sparsity>(Sparsity(n_out_, n_in_));

    // Type of derivative calculations enabled
    enable_forward_ = enable_forward_ && has_forward(1);
    enable_reverse_ = enable_reverse_ && has_reverse(1);
    enable_jacobian_ = enable_jacobian_ && has_jacobian();
    enable_fd_ = enable_fd_ && !enable_forward_;

    alloc_arg(0);
    alloc_res(0);
  }

  std::string FunctionInternal::get_name_in(casadi_int i) {
    return "i" + str(i);
  }

  std::string FunctionInternal::get_name_out(casadi_int i) {
    return "o" + str(i);
  }

  void FunctionInternal::finalize(const Dict& opts) {
    if (jit_) {
      string jit_name = "jit_tmp";
      if (has_codegen()) {
        if (verbose_) casadi_message("Codegenerating function '" + name_ + "'.");
        // JIT everything
        CodeGenerator gen(jit_name);
        gen.add(self());
        if (verbose_) casadi_message("Compiling function '" + name_ + "'..");
        compiler_ = Importer(gen.generate(), compilerplugin_, jit_options_);
        if (verbose_) casadi_message("Compiling function '" + name_ + "' done.");
        // Try to load
        eval_ = (eval_t)compiler_.get_function(name_);
        casadi_assert(eval_!=nullptr, "Cannot load JIT'ed function.");
      } else {
        // Just jit dependencies
        jit_dependencies(jit_name);
      }
    }
    // Finalize base classes
    ProtoFunction::finalize(opts);
  }

  void ProtoFunction::finalize(const Dict& opts) {
    // Create memory object
    casadi_int mem = checkout();
    casadi_assert_dev(mem==0);
  }

  int FunctionInternal::
  eval_gen(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    if (eval_) {
      return eval_(arg, res, iw, w, mem);
    } else {
      return eval(arg, res, iw, w, mem);
    }
  }

  void FunctionInternal::print_dimensions(ostream &stream) const {
    stream << " Number of inputs: " << n_in_ << endl;
    for (casadi_int i=0; i<n_in_; ++i) {
      stream << "  Input " << i  << " (\"" << name_in_[i] << "\"): "
             << sparsity_in_[i].dim() << endl;
    }
    stream << " Number of outputs: " << n_out_ << endl;
    for (casadi_int i=0; i<n_out_; ++i) {
      stream << "  Output " << i  << " (\"" << name_out_[i] << "\"): "
             << sparsity_out_[i].dim() << endl;
    }
  }

  void FunctionInternal::print_options(std::ostream &stream) const {
    get_options().print_all(stream);
  }

  void FunctionInternal::print_option(const std::string &name, std::ostream &stream) const {
    get_options().print_one(name, stream);
  }

  std::vector<std::string> FunctionInternal::get_free() const {
    casadi_assert_dev(!has_free());
    return std::vector<std::string>();
  }

  std::string FunctionInternal::definition() const {
    stringstream s;

    // Print name
    s << name_ << ":(";
    // Print input arguments
    for (casadi_int i=0; i<n_in_; ++i) {
      s << name_in_[i] << sparsity_in_[i].postfix_dim() << (i==n_in_-1 ? "" : ",");
    }
    s << ")->(";
    // Print output arguments
    for (casadi_int i=0; i<n_out_; ++i) {
      s << name_out_[i] << sparsity_out_[i].postfix_dim() << (i==n_out_-1 ? "" : ",");
    }
    s << ")";

    return s.str();
  }

  void FunctionInternal::disp(ostream &stream, bool more) const {
    stream << definition() << " " << class_name();
    if (more) {
      stream << endl;
      disp_more(stream);
    }
  }

  bool FunctionInternal::incache(const std::string& fname, Function& f) const {
    auto it = cache_.find(fname);
    if (it!=cache_.end() && it->second.alive()) {
      f = shared_cast<Function>(it->second.shared());
      return true;
    } else {
      return false;
    }
  }

  void FunctionInternal::tocache(const Function& f) const {
    // Add to cache
    cache_.insert(make_pair(f.name(), f));
    // Remove a lost reference, if any, to prevent uncontrolled growth
    for (auto it = cache_.begin(); it!=cache_.end(); ++it) {
      if (!it->second.alive()) {
        cache_.erase(it);
        break; // just one dead reference is enough
      }
    }
  }

  Function FunctionInternal::map(casadi_int n, const std::string& parallelization) const {
    Function f;
    if (parallelization=="serial") {
      // Serial maps are cached
      string fname = "map" + str(n) + "_" + name_;
      if (!incache(fname, f)) {
        // Create new serial map
        f = Map::create(parallelization, self(), n);
        casadi_assert_dev(f.name()==fname);
        // Save in cache
        tocache(f);
      }
    } else {
      // Non-serial maps are not cached
      f = Map::create(parallelization, self(), n);
    }
    return f;
  }

  Function FunctionInternal::wrap() const {
    Function f;
    string fname = "wrap_" + name_;
    if (!incache(fname, f)) {
      // Options
      Dict opts;
      opts["derivative_of"] = derivative_of_;
      opts["ad_weight"] = ad_weight();
      opts["ad_weight_sp"] = sp_weight();
      opts["max_num_dir"] = max_num_dir_;
      // Wrap the function
      vector<MX> arg = mx_in();
      vector<MX> res = self()(arg);
      f = Function(fname, arg, res, name_in_, name_out_, opts);
      // Save in cache
      tocache(f);
    }
    return f;
  }

  std::vector<MX> FunctionInternal::symbolic_output(const std::vector<MX>& arg) const {
    return self()(arg);
  }

  /// \cond INTERNAL

  void bvec_toggle(bvec_t* s, casadi_int begin, casadi_int end, casadi_int j) {
    for (casadi_int i=begin; i<end; ++i) {
      s[i] ^= (bvec_t(1) << j);
    }
  }

  void bvec_clear(bvec_t* s, casadi_int begin, casadi_int end) {
    for (casadi_int i=begin; i<end; ++i) {
      s[i] = 0;
    }
  }


  void bvec_or(bvec_t* s, bvec_t & r, casadi_int begin, casadi_int end) {
    r = 0;
    for (casadi_int i=begin; i<end; ++i) r |= s[i];
  }
  /// \endcond

  // Traits
  template<bool fwd> struct JacSparsityTraits {};
  template<> struct JacSparsityTraits<true> {
    typedef const bvec_t* arg_t;
    static inline void sp(const FunctionInternal *f,
                          const bvec_t** arg, bvec_t** res,
                          casadi_int* iw, bvec_t* w, void* mem) {
      f->sp_forward(arg, res, iw, w, mem);
    }
  };
  template<> struct JacSparsityTraits<false> {
    typedef bvec_t* arg_t;
    static inline void sp(const FunctionInternal *f,
                          bvec_t** arg, bvec_t** res,
                          casadi_int* iw, bvec_t* w, void* mem) {
      f->sp_reverse(arg, res, iw, w, mem);
    }
  };

  template<bool fwd>
  Sparsity FunctionInternal::
  getJacSparsityGen(casadi_int iind, casadi_int oind, bool symmetric,
      casadi_int gr_i, casadi_int gr_o) const {
    // Number of nonzero inputs and outputs
    casadi_int nz_in = nnz_in(iind);
    casadi_int nz_out = nnz_out(oind);

    // Evaluation buffers
    vector<typename JacSparsityTraits<fwd>::arg_t> arg(sz_arg(), nullptr);
    vector<bvec_t*> res(sz_res(), nullptr);
    vector<casadi_int> iw(sz_iw());
    vector<bvec_t> w(sz_w(), 0);

    // Seeds and sensitivities
    vector<bvec_t> seed(nz_in, 0);
    arg[iind] = get_ptr(seed);
    vector<bvec_t> sens(nz_out, 0);
    res[oind] = get_ptr(sens);
    if (!fwd) std::swap(seed, sens);

    // Number of forward sweeps we must make
    casadi_int nsweep = seed.size() / bvec_size;
    if (seed.size() % bvec_size) nsweep++;

    // Print
    if (verbose_) {
      casadi_message(str(nsweep) + string(fwd ? " forward" : " reverse") + " sweeps "
                     "needed for " + str(seed.size()) + " directions");
    }

    // Progress
    casadi_int progress = -10;

    // Temporary vectors
    std::vector<casadi_int> jcol, jrow;

    // Loop over the variables, bvec_size variables at a time
    for (casadi_int s=0; s<nsweep; ++s) {

      // Print progress
      if (verbose_) {
        casadi_int progress_new = (s*100)/nsweep;
        // Print when entering a new decade
        if (progress_new / 10 > progress / 10) {
          progress = progress_new;
          casadi_message(str(progress) + " %");
        }
      }

      // Nonzero offset
      casadi_int offset = s*bvec_size;

      // Number of local seed directions
      casadi_int ndir_local = seed.size()-offset;
      ndir_local = std::min(static_cast<casadi_int>(bvec_size), ndir_local);

      for (casadi_int i=0; i<ndir_local; ++i) {
        seed[offset+i] |= bvec_t(1)<<i;
      }

      // Propagate the dependencies
      JacSparsityTraits<fwd>::sp(this, get_ptr(arg), get_ptr(res),
                                  get_ptr(iw), get_ptr(w), memory(0));

      // Loop over the nonzeros of the output
      for (casadi_int el=0; el<sens.size(); ++el) {

        // Get the sparsity sensitivity
        bvec_t spsens = sens[el];

        if (!fwd) {
          // Clear the sensitivities for the next sweep
          sens[el] = 0;
        }

        // If there is a dependency in any of the directions
        if (spsens!=0) {

          // Loop over seed directions
          for (casadi_int i=0; i<ndir_local; ++i) {

            // If dependents on the variable
            if ((bvec_t(1) << i) & spsens) {
              // Add to pattern
              jcol.push_back(el);
              jrow.push_back(i+offset);
            }
          }
        }
      }

      // Remove the seeds
      for (casadi_int i=0; i<ndir_local; ++i) {
        seed[offset+i] = 0;
      }
    }

    // Construct sparsity pattern and return
    if (!fwd) swap(jrow, jcol);
    Sparsity ret = Sparsity::triplet(nz_out, nz_in, jcol, jrow);
    if (verbose_) {
      casadi_message("Formed Jacobian sparsity pattern (dimension " + str(ret.size()) + ", "
          + str(ret.nnz()) + " (" + str(ret.density()) + " %) nonzeros.");
    }
    return ret;
  }

  Sparsity FunctionInternal::
  getJacSparsityHierarchicalSymm(casadi_int iind, casadi_int oind) const {
    casadi_assert_dev(has_spfwd());

    // Number of nonzero inputs
    casadi_int nz = nnz_in(iind);
    casadi_assert_dev(nz==nnz_out(oind));

    // Evaluation buffers
    vector<const bvec_t*> arg(sz_arg(), nullptr);
    vector<bvec_t*> res(sz_res(), nullptr);
    vector<casadi_int> iw(sz_iw());
    vector<bvec_t> w(sz_w());

    // Seeds
    vector<bvec_t> seed(nz, 0);
    arg[iind] = get_ptr(seed);

    // Sensitivities
    vector<bvec_t> sens(nz, 0);
    res[oind] = get_ptr(sens);

    // Sparsity triplet accumulator
    std::vector<casadi_int> jcol, jrow;

    // Cols/rows of the coarse blocks
    std::vector<casadi_int> coarse(2, 0); coarse[1] = nz;

    // Cols/rows of the fine blocks
    std::vector<casadi_int> fine;

    // In each iteration, subdivide each coarse block in this many fine blocks
    casadi_int subdivision = bvec_size;

    Sparsity r = Sparsity::dense(1, 1);

    // The size of a block
    casadi_int granularity = nz;

    casadi_int nsweeps = 0;

    bool hasrun = false;

    while (!hasrun || coarse.size()!=nz+1) {
      if (verbose_) casadi_message("Block size: " + str(granularity));

      // Clear the sparsity triplet acccumulator
      jcol.clear();
      jrow.clear();

      // Clear the fine block structure
      fine.clear();

      Sparsity D = r.star_coloring();

      if (verbose_) {
        casadi_message("Star coloring on " + str(r.dim()) + ": "
          + str(D.size2()) + " <-> " + str(D.size1()));
      }

      // Clear the seeds
      fill(seed.begin(), seed.end(), 0);

      // Subdivide the coarse block
      for (casadi_int k=0; k<coarse.size()-1; ++k) {
        casadi_int diff = coarse[k+1]-coarse[k];
        casadi_int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<casadi_int> temp = range(coarse[k], coarse[k+1], new_diff);
        fine.insert(fine.end(), temp.begin(), temp.end());
      }
      if (fine.back()!=coarse.back()) fine.push_back(coarse.back());

      granularity = fine[1] - fine[0];

      // The index into the bvec bit vector
      casadi_int bvec_i = 0;

      // Create lookup tables for the fine blocks
      std::vector<casadi_int> fine_lookup = lookupvector(fine, nz+1);

      // Triplet data used as a lookup table
      std::vector<casadi_int> lookup_col;
      std::vector<casadi_int> lookup_row;
      std::vector<casadi_int> lookup_value;

      // The maximum number of fine blocks contained in one coarse block
      casadi_int n_fine_blocks_max = 0;
      for (casadi_int i=0;i<coarse.size()-1;++i) {
        casadi_int del = fine_lookup[coarse[i+1]]-fine_lookup[coarse[i]];
        n_fine_blocks_max = std::max(n_fine_blocks_max, del);
      }

      // Loop over all coarse seed directions from the coloring
      for (casadi_int csd=0; csd<D.size2(); ++csd) {


        casadi_int fci_offset = 0;
        casadi_int fci_cap = bvec_size-bvec_i;

        // Flag to indicate if all fine blocks have been handled
        bool f_finished = false;

        // Loop while not finished
        while (!f_finished) {

          // Loop over all coarse rows that are found in the coloring for this coarse seed direction
          for (casadi_int k=D.colind(csd); k<D.colind(csd+1); ++k) {
            casadi_int cci = D.row(k);

            // The first and last rows of the fine block
            casadi_int fci_start = fine_lookup[coarse[cci]];
            casadi_int fci_end   = fine_lookup[coarse[cci+1]];

            // Local counter that modifies index into bvec
            casadi_int bvec_i_mod = 0;

            casadi_int value = -bvec_i + fci_offset + fci_start;

            //casadi_assert_dev(value>=0);

            // Loop over the rows of the fine block
            for (casadi_int fci = fci_offset;fci<min(fci_end-fci_start, fci_cap);++fci) {

              // Loop over the coarse block cols that appear in the
              // coloring for the current coarse seed direction
              for (casadi_int cri=r.colind(cci);cri<r.colind(cci+1);++cri) {
                lookup_col.push_back(r.row(cri));
                lookup_row.push_back(bvec_i+bvec_i_mod);
                lookup_value.push_back(value);
              }

              // Toggle on seeds
              bvec_toggle(get_ptr(seed), fine[fci+fci_start], fine[fci+fci_start+1],
                          bvec_i+bvec_i_mod);
              bvec_i_mod++;
            }
          }

          // Bump bvec_i for next major coarse direction
          bvec_i+= min(n_fine_blocks_max, fci_cap);

          // Check if bvec buffer is full
          if (bvec_i==bvec_size || csd==D.size2()-1) {
            // Calculate sparsity for bvec_size directions at once

            // Statistics
            nsweeps+=1;

            // Construct lookup table
            IM lookup = IM::triplet(lookup_row, lookup_col, lookup_value,
                                    bvec_size, coarse.size());

            std::reverse(lookup_col.begin(), lookup_col.end());
            std::reverse(lookup_row.begin(), lookup_row.end());
            std::reverse(lookup_value.begin(), lookup_value.end());
            IM duplicates =
              IM::triplet(lookup_row, lookup_col, lookup_value, bvec_size, coarse.size())
              - lookup;
            duplicates = sparsify(duplicates);
            lookup(duplicates.sparsity()) = -bvec_size;

            // Propagate the dependencies
            sp_forward(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), nullptr);

            // Temporary bit work vector
            bvec_t spsens;

            // Loop over the cols of coarse blocks
            for (casadi_int cri=0; cri<coarse.size()-1; ++cri) {

              // Loop over the cols of fine blocks within the current coarse block
              for (casadi_int fri=fine_lookup[coarse[cri]];fri<fine_lookup[coarse[cri+1]];++fri) {
                // Lump individual sensitivities together into fine block
                bvec_or(get_ptr(sens), spsens, fine[fri], fine[fri+1]);

                // Loop over all bvec_bits
                for (casadi_int bvec_i=0;bvec_i<bvec_size;++bvec_i) {
                  if (spsens & (bvec_t(1) << bvec_i)) {
                    // if dependency is found, add it to the new sparsity pattern
                    casadi_int ind = lookup.sparsity().get_nz(bvec_i, cri);
                    if (ind==-1) continue;
                    casadi_int lk = lookup->at(ind);
                    if (lk>-bvec_size) {
                      jrow.push_back(bvec_i+lk);
                      jcol.push_back(fri);
                      jrow.push_back(fri);
                      jcol.push_back(bvec_i+lk);
                    }
                  }
                }
              }
            }

            // Clear the forward seeds/adjoint sensitivities, ready for next bvec sweep
            fill(seed.begin(), seed.end(), 0);

            // Clean lookup table
            lookup_col.clear();
            lookup_row.clear();
            lookup_value.clear();
          }

          if (n_fine_blocks_max>fci_cap) {
            fci_offset += min(n_fine_blocks_max, fci_cap);
            bvec_i = 0;
            fci_cap = bvec_size;
          } else {
            f_finished = true;
          }
        }
      }

      // Construct fine sparsity pattern
      r = Sparsity::triplet(fine.size()-1, fine.size()-1, jrow, jcol);

      // There may be false positives here that are not present
      // in the reverse mode that precedes it.
      // This can lead to an assymetrical result
      //  cf. #1522
      r=r*r.T();

      coarse = fine;
      hasrun = true;
    }
    if (verbose_) {
      casadi_message("Number of sweeps: " + str(nsweeps));
      casadi_message("Formed Jacobian sparsity pattern (dimension " + str(r.size()) +
          ", " + str(r.nnz()) + " (" + str(r.density()) + " %) nonzeros.");
    }

    return r.T();
  }

  Sparsity FunctionInternal::
  getJacSparsityHierarchical(casadi_int iind, casadi_int oind) const {
    // Number of nonzero inputs
    casadi_int nz_in = nnz_in(iind);

    // Number of nonzero outputs
    casadi_int nz_out = nnz_out(oind);

    // Seeds and sensitivities
    vector<bvec_t> s_in(nz_in, 0);
    vector<bvec_t> s_out(nz_out, 0);

    // Evaluation buffers
    vector<const bvec_t*> arg_fwd(sz_arg(), nullptr);
    vector<bvec_t*> arg_adj(sz_arg(), nullptr);
    arg_fwd[iind] = arg_adj[iind] = get_ptr(s_in);
    vector<bvec_t*> res(sz_res(), nullptr);
    res[oind] = get_ptr(s_out);
    vector<casadi_int> iw(sz_iw());
    vector<bvec_t> w(sz_w());

    // Sparsity triplet accumulator
    std::vector<casadi_int> jcol, jrow;

    // Cols of the coarse blocks
    std::vector<casadi_int> coarse_col(2, 0); coarse_col[1] = nz_out;
    // Rows of the coarse blocks
    std::vector<casadi_int> coarse_row(2, 0); coarse_row[1] = nz_in;

    // Cols of the fine blocks
    std::vector<casadi_int> fine_col;

    // Rows of the fine blocks
    std::vector<casadi_int> fine_row;

    // In each iteration, subdivide each coarse block in this many fine blocks
    casadi_int subdivision = bvec_size;

    Sparsity r = Sparsity::dense(1, 1);

    // The size of a block
    casadi_int granularity_row = nz_in;
    casadi_int granularity_col = nz_out;

    bool use_fwd = true;

    casadi_int nsweeps = 0;

    bool hasrun = false;

    // Get weighting factor
    double sp_w = sp_weight();

    // Lookup table for bvec_t
    std::vector<bvec_t> bvec_lookup;
    bvec_lookup.reserve(bvec_size);
    for (casadi_int i=0;i<bvec_size;++i) {
      bvec_lookup.push_back(bvec_t(1) << i);
    }

    while (!hasrun || coarse_col.size()!=nz_out+1 || coarse_row.size()!=nz_in+1) {
      if (verbose_) {
        casadi_message("Block size: " + str(granularity_col) + " x " + str(granularity_row));
      }

      // Clear the sparsity triplet acccumulator
      jcol.clear();
      jrow.clear();

      // Clear the fine block structure
      fine_row.clear();
      fine_col.clear();

      // r transpose will be needed in the algorithm
      Sparsity rT = r.T();

      /**       Decide which ad_mode to take           */

      // Forward mode
      Sparsity D1 = rT.uni_coloring(r);
      // Adjoint mode
      Sparsity D2 = r.uni_coloring(rT);
      if (verbose_) {
        casadi_message("Coloring on " + str(r.dim()) + " (fwd seeps: " + str(D1.size2()) +
                 " , adj sweeps: " + str(D2.size1()) + ")");
      }

      // Use whatever required less colors if we tried both (with preference to forward mode)
      double fwd_cost = static_cast<double>(use_fwd ? granularity_row : granularity_col) *
        sp_w*static_cast<double>(D1.size2());
      double adj_cost = static_cast<double>(use_fwd ? granularity_col : granularity_row) *
        (1-sp_w)*static_cast<double>(D2.size2());
      use_fwd = fwd_cost <= adj_cost;
      if (verbose_) {
        casadi_message(string(use_fwd ? "Forward" : "Reverse") + " mode chosen "
            "(fwd cost: " + str(fwd_cost) + ", adj cost: " + str(adj_cost) + ")");
      }

      // Get seeds and sensitivities
      bvec_t* seed_v = use_fwd ? get_ptr(s_in) : get_ptr(s_out);
      bvec_t* sens_v = use_fwd ? get_ptr(s_out) : get_ptr(s_in);

      // The number of zeros in the seed and sensitivity directions
      casadi_int nz_seed = use_fwd ? nz_in  : nz_out;
      casadi_int nz_sens = use_fwd ? nz_out : nz_in;

      // Clear the seeds
      for (casadi_int i=0; i<nz_seed; ++i) seed_v[i]=0;

      // Choose the active jacobian coloring scheme
      Sparsity D = use_fwd ? D1 : D2;

      // Adjoint mode amounts to swapping
      if (!use_fwd) {
        std::swap(coarse_col, coarse_row);
        std::swap(granularity_col, granularity_row);
        std::swap(r, rT);
      }

      // Subdivide the coarse block cols
      for (casadi_int k=0;k<coarse_col.size()-1;++k) {
        casadi_int diff = coarse_col[k+1]-coarse_col[k];
        casadi_int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<casadi_int> temp = range(coarse_col[k], coarse_col[k+1], new_diff);
        fine_col.insert(fine_col.end(), temp.begin(), temp.end());
      }
      // Subdivide the coarse block rows
      for (casadi_int k=0;k<coarse_row.size()-1;++k) {
        casadi_int diff = coarse_row[k+1]-coarse_row[k];
        casadi_int new_diff = diff/subdivision;
        if (diff%subdivision>0) new_diff++;
        std::vector<casadi_int> temp = range(coarse_row[k], coarse_row[k+1], new_diff);
        fine_row.insert(fine_row.end(), temp.begin(), temp.end());
      }
      if (fine_row.back()!=coarse_row.back()) fine_row.push_back(coarse_row.back());
      if (fine_col.back()!=coarse_col.back()) fine_col.push_back(coarse_col.back());

      granularity_col = fine_col[1] - fine_col[0];
      granularity_row = fine_row[1] - fine_row[0];

      // The index into the bvec bit vector
      casadi_int bvec_i = 0;

      // Create lookup tables for the fine blocks
      std::vector<casadi_int> fine_col_lookup = lookupvector(fine_col, nz_sens+1);
      std::vector<casadi_int> fine_row_lookup = lookupvector(fine_row, nz_seed+1);

      // Triplet data used as a lookup table
      std::vector<casadi_int> lookup_col;
      std::vector<casadi_int> lookup_row;
      std::vector<casadi_int> lookup_value;


      // The maximum number of fine blocks contained in one coarse block
      casadi_int n_fine_blocks_max = 0;
      for (casadi_int i=0;i<coarse_row.size()-1;++i) {
        casadi_int del = fine_row_lookup[coarse_row[i+1]]-fine_row_lookup[coarse_row[i]];
        n_fine_blocks_max = std::max(n_fine_blocks_max, del);
      }

      // Loop over all coarse seed directions from the coloring
      for (casadi_int csd=0; csd<D.size2(); ++csd) {

        casadi_int fci_offset = 0;
        casadi_int fci_cap = bvec_size-bvec_i;

        // Flag to indicate if all fine blocks have been handled
        bool f_finished = false;

        // Loop while not finished
        while (!f_finished) {

          // Loop over all coarse rows that are found in the coloring for this coarse seed direction
          for (casadi_int k=D.colind(csd); k<D.colind(csd+1); ++k) {
            casadi_int cci = D.row(k);

            // The first and last rows of the fine block
            casadi_int fci_start = fine_row_lookup[coarse_row[cci]];
            casadi_int fci_end   = fine_row_lookup[coarse_row[cci+1]];

            // Local counter that modifies index into bvec
            casadi_int bvec_i_mod = 0;

            casadi_int value = -bvec_i + fci_offset + fci_start;

            // Loop over the rows of the fine block
            for (casadi_int fci = fci_offset; fci<min(fci_end-fci_start, fci_cap); ++fci) {

              // Loop over the coarse block cols that appear in the coloring
              // for the current coarse seed direction
              for (casadi_int cri=rT.colind(cci);cri<rT.colind(cci+1);++cri) {
                lookup_col.push_back(rT.row(cri));
                lookup_row.push_back(bvec_i+bvec_i_mod);
                lookup_value.push_back(value);
              }

              // Toggle on seeds
              bvec_toggle(seed_v, fine_row[fci+fci_start], fine_row[fci+fci_start+1],
                          bvec_i+bvec_i_mod);
              bvec_i_mod++;
            }
          }

          // Bump bvec_i for next major coarse direction
          bvec_i+= min(n_fine_blocks_max, fci_cap);

          // Check if bvec buffer is full
          if (bvec_i==bvec_size || csd==D.size2()-1) {
            // Calculate sparsity for bvec_size directions at once

            // Statistics
            nsweeps+=1;

            // Construct lookup table
            IM lookup = IM::triplet(lookup_row, lookup_col, lookup_value, bvec_size,
                                    coarse_col.size());

            // Propagate the dependencies
            if (use_fwd) {
              sp_forward(get_ptr(arg_fwd), get_ptr(res), get_ptr(iw), get_ptr(w), memory(0));
            } else {
              fill(w.begin(), w.end(), 0);
              sp_reverse(get_ptr(arg_adj), get_ptr(res), get_ptr(iw), get_ptr(w), memory(0));
            }

            // Temporary bit work vector
            bvec_t spsens;

            // Loop over the cols of coarse blocks
            for (casadi_int cri=0;cri<coarse_col.size()-1;++cri) {

              // Loop over the cols of fine blocks within the current coarse block
              for (casadi_int fri=fine_col_lookup[coarse_col[cri]];
                   fri<fine_col_lookup[coarse_col[cri+1]];++fri) {
                // Lump individual sensitivities together into fine block
                bvec_or(sens_v, spsens, fine_col[fri], fine_col[fri+1]);

                // Next iteration if no sparsity
                if (!spsens) continue;

                // Loop over all bvec_bits
                for (casadi_int bvec_i=0;bvec_i<bvec_size;++bvec_i) {
                  if (spsens & bvec_lookup[bvec_i]) {
                    // if dependency is found, add it to the new sparsity pattern
                    casadi_int ind = lookup.sparsity().get_nz(bvec_i, cri);
                    if (ind==-1) continue;
                    jrow.push_back(bvec_i+lookup->at(ind));
                    jcol.push_back(fri);
                  }
                }
              }
            }

            // Clear the forward seeds/adjoint sensitivities, ready for next bvec sweep
            fill(s_in.begin(), s_in.end(), 0);

            // Clear the adjoint seeds/forward sensitivities, ready for next bvec sweep
            fill(s_out.begin(), s_out.end(), 0);

            // Clean lookup table
            lookup_col.clear();
            lookup_row.clear();
            lookup_value.clear();
          }

          if (n_fine_blocks_max>fci_cap) {
            fci_offset += min(n_fine_blocks_max, fci_cap);
            bvec_i = 0;
            fci_cap = bvec_size;
          } else {
            f_finished = true;
          }

        }

      }

      // Swap results if adjoint mode was used
      if (use_fwd) {
        // Construct fine sparsity pattern
        r = Sparsity::triplet(fine_row.size()-1, fine_col.size()-1, jrow, jcol);
        coarse_col = fine_col;
        coarse_row = fine_row;
      } else {
        // Construct fine sparsity pattern
        r = Sparsity::triplet(fine_col.size()-1, fine_row.size()-1, jcol, jrow);
        coarse_col = fine_row;
        coarse_row = fine_col;
      }
      hasrun = true;
    }
    if (verbose_) {
      casadi_message("Number of sweeps: " + str(nsweeps));
      casadi_message("Formed Jacobian sparsity pattern (dimension " + str(r.size()) + ", " +
          str(r.nnz()) + " (" + str(r.density()) + " %) nonzeros.");
    }

    return r.T();
  }

  Sparsity FunctionInternal::getJacSparsity(casadi_int iind, casadi_int oind,
      bool symmetric) const {
    // Check if we are able to propagate dependencies through the function
    if (has_spfwd() || has_sprev()) {
      Sparsity sp;
      if (nnz_in(iind)>3*bvec_size && nnz_out(oind)>3*bvec_size &&
            GlobalOptions::hierarchical_sparsity) {
        if (symmetric) {
          sp = getJacSparsityHierarchicalSymm(iind, oind);
        } else {
          sp = getJacSparsityHierarchical(iind, oind);
        }
      } else {
        // Number of nonzero inputs and outputs
        casadi_int nz_in = nnz_in(iind);
        casadi_int nz_out = nnz_out(oind);

        // Number of forward sweeps we must make
        casadi_int nsweep_fwd = nz_in/bvec_size;
        if (nz_in%bvec_size) nsweep_fwd++;

        // Number of adjoint sweeps we must make
        casadi_int nsweep_adj = nz_out/bvec_size;
        if (nz_out%bvec_size) nsweep_adj++;

        // Get weighting factor
        double w = sp_weight();

        // Use forward mode?
        if (w*static_cast<double>(nsweep_fwd) <= (1-w)*static_cast<double>(nsweep_adj)) {
          sp = getJacSparsityGen<true>(iind, oind, false);
        } else {
          sp = getJacSparsityGen<false>(iind, oind, false);
        }
      }
      // There may be false positives here that are not present
      // in the reverse mode that precedes it.
      // This can lead to an assymetrical result
      //  cf. #1522
      if (symmetric) sp=sp*sp.T();
      return sp;
    } else {
      // Dense sparsity by default
      return Sparsity::dense(nnz_out(oind), nnz_in(iind));
    }
  }

  Sparsity& FunctionInternal::
  sparsity_jac(casadi_int iind, casadi_int oind, bool compact, bool symmetric) const {
    // Get an owning reference to the block
    Sparsity jsp = compact ? jac_sparsity_compact_.elem(oind, iind)
        : jac_sparsity_.elem(oind, iind);

    // Generate, if null
    if (jsp.is_null()) {
      if (compact) {

        // Use internal routine to determine sparsity
        jsp = getJacSparsity(iind, oind, symmetric);

      } else {

        // Get the compact sparsity pattern
        Sparsity sp = sparsity_jac(iind, oind, true, symmetric);

        // Enlarge if sparse output
        if (numel_out(oind)!=sp.size1()) {
          casadi_assert_dev(sp.size1()==nnz_out(oind));

          // New row for each old row
          vector<casadi_int> row_map = sparsity_out(oind).find();

          // Insert rows
          sp.enlargeRows(numel_out(oind), row_map);
        }

        // Enlarge if sparse input
        if (numel_in(iind)!=sp.size2()) {
          casadi_assert_dev(sp.size2()==nnz_in(iind));

          // New column for each old column
          vector<casadi_int> col_map = sparsity_in(iind).find();

          // Insert columns
          sp.enlargeColumns(numel_in(iind), col_map);
        }

        // Save
        jsp = sp;
      }
    }

    // If still null, not dependent
    if (jsp.is_null()) {
      jsp = Sparsity(nnz_out(oind), nnz_in(iind));
    }

    // Return a reference to the block
    Sparsity& jsp_ref = compact ? jac_sparsity_compact_.elem(oind, iind) :
        jac_sparsity_.elem(oind, iind);
    jsp_ref = jsp;
    return jsp_ref;
  }

  void FunctionInternal::get_partition(casadi_int iind, casadi_int oind, Sparsity& D1, Sparsity& D2,
                                       bool compact, bool symmetric,
                                       bool allow_forward, bool allow_reverse) const {
    if (verbose_) casadi_message(name_ + "::get_partition");
    casadi_assert(allow_forward || allow_reverse, "Inconsistent options");

    // Sparsity pattern with transpose
    Sparsity &AT = sparsity_jac(iind, oind, compact, symmetric);
    Sparsity A = symmetric ? AT : AT.T();

    // Get seed matrices by graph coloring
    if (symmetric) {
      casadi_assert_dev(enable_forward_ || enable_fd_);
      casadi_assert_dev(allow_forward);

      // Star coloring if symmetric
      if (verbose_) casadi_message("FunctionInternal::getPartition star_coloring");
      D1 = A.star_coloring();
      if (verbose_) {
        casadi_message("Star coloring completed: " + str(D1.size2())
          + " directional derivatives needed ("
                 + str(A.size1()) + " without coloring).");
      }

    } else {
      casadi_assert_dev(enable_forward_ || enable_fd_ || enable_reverse_);
      // Get weighting factor
      double w = ad_weight();

      // Which AD mode?
      if (w==1) allow_forward = false;
      if (w==0) allow_reverse = false;
      casadi_assert(allow_forward || allow_reverse, "Conflicting ad weights");

      // Best coloring encountered so far (relatively tight upper bound)
      double best_coloring = numeric_limits<double>::infinity();

      // Test forward mode first?
      bool test_fwd_first = allow_forward && w*static_cast<double>(A.size1()) <=
        (1-w)*static_cast<double>(A.size2());
      casadi_int mode_fwd = test_fwd_first ? 0 : 1;

      // Test both coloring modes
      for (casadi_int mode=0; mode<2; ++mode) {
        // Is this the forward mode?
        bool fwd = mode==mode_fwd;

        // Skip?
        if (!allow_forward && fwd) continue;
        if (!allow_reverse && !fwd) continue;

        // Perform the coloring
        if (fwd) {
          if (verbose_) casadi_message("Unidirectional coloring (forward mode)");
          bool d = best_coloring>=w*static_cast<double>(A.size1());
          casadi_int max_colorings_to_test =
            d ? A.size1() : static_cast<casadi_int>(floor(best_coloring/w));
          D1 = AT.uni_coloring(A, max_colorings_to_test);
          if (D1.is_null()) {
            if (verbose_) {
              casadi_message("Forward mode coloring interrupted (more than "
                             + str(max_colorings_to_test) + " needed).");
            }
          } else {
            if (verbose_) {
              casadi_message("Forward mode coloring completed: "
                             + str(D1.size2()) + " directional derivatives needed ("
                             + str(A.size1()) + " without coloring).");
            }
            D2 = Sparsity();
            best_coloring = w*static_cast<double>(D1.size2());
          }
        } else {
          if (verbose_) casadi_message("Unidirectional coloring (adjoint mode)");
          bool d = best_coloring>=(1-w)*static_cast<double>(A.size2());
          casadi_int max_colorings_to_test =
            d ? A.size2() : static_cast<casadi_int>(floor(best_coloring/(1-w)));

          D2 = A.uni_coloring(AT, max_colorings_to_test);
          if (D2.is_null()) {
            if (verbose_) {
              casadi_message("Adjoint mode coloring interrupted (more than "
                            + str(max_colorings_to_test) + " needed).");
            }
          } else {
            if (verbose_) {
              casadi_message("Adjoint mode coloring completed: "
                             + str(D2.size2()) + " directional derivatives needed ("
                             + str(A.size2()) + " without coloring).");
            }
            D1 = Sparsity();
            best_coloring = (1-w)*static_cast<double>(D2.size2());
          }
        }
      }

    }
  }

  std::vector<DM> FunctionInternal::eval_dm(const std::vector<DM>& arg) const {
    casadi_error("'eval', 'eval_dm' not defined for " + class_name());
  }

  int FunctionInternal::
  eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w, void* mem) const {
    casadi_error("'eval_sx' not defined for " + class_name());
  }

  Function FunctionInternal::forward(casadi_int nfwd) const {
    casadi_assert_dev(nfwd>=0);
    // Used wrapped function if forward not available
    if (!enable_forward_ && !enable_fd_) {
      // Derivative information must be available
      casadi_assert(has_derivative(), "Derivatives cannot be calculated for " + name_);
      return wrap().forward(nfwd);
    }
    // Retrieve/generate cached
    Function f;
    string fname = "fwd" + str(nfwd) + "_" + name_;
    if (!incache(fname, f)) {
      casadi_int i;
      // Names of inputs
      std::vector<std::string> inames;
      for (i=0; i<n_in_; ++i) inames.push_back(name_in_[i]);
      for (i=0; i<n_out_; ++i) inames.push_back("out_" + name_out_[i]);
      for (i=0; i<n_in_; ++i) inames.push_back("fwd_" + name_in_[i]);
      // Names of outputs
      std::vector<std::string> onames;
      for (i=0; i<n_out_; ++i) onames.push_back("fwd_" + name_out_[i]);
      // Options
      Dict opts;
      if (!enable_forward_) opts = fd_options_;
      opts["max_num_dir"] = max_num_dir_;
      opts["derivative_of"] = self();
      // Generate derivative function
      casadi_assert_dev(enable_forward_ || enable_fd_);
      if (enable_forward_) {
        f = get_forward(nfwd, fname, inames, onames, opts);
      } else {
        // Get FD method
        if (fd_method_.empty() || fd_method_=="central") {
          f = Function::create(new CentralDiff(fname, nfwd), opts);
        } else if (fd_method_=="forward") {
          f = Function::create(new ForwardDiff(fname, nfwd), opts);
        } else if (fd_method_=="backward") {
          f = Function::create(new BackwardDiff(fname, nfwd), opts);
        } else if (fd_method_=="smoothing") {
          f = Function::create(new Smoothing(fname, nfwd), opts);
        } else {
          casadi_error("Unknown 'fd_method': " + fd_method_);
        }
      }
      // Consistency check for inputs
      casadi_assert_dev(f.n_in()==n_in_ + n_out_ + n_in_);
      casadi_int ind=0;
      for (i=0; i<n_in_; ++i) f.assert_size_in(ind++, size1_in(i), size2_in(i));
      for (i=0; i<n_out_; ++i) f.assert_size_in(ind++, size1_out(i), size2_out(i));
      for (i=0; i<n_in_; ++i) f.assert_size_in(ind++, size1_in(i), nfwd*size2_in(i));
      // Consistency check for outputs
      casadi_assert_dev(f.n_out()==n_out_);
      for (i=0; i<n_out_; ++i) f.assert_size_out(i, size1_out(i), nfwd*size2_out(i));
      // Save to cache
      tocache(f);
    }
    return f;
  }

  Function FunctionInternal::reverse(casadi_int nadj) const {
    casadi_assert_dev(nadj>=0);
    // Used wrapped function if reverse not available
    if (!enable_reverse_) {
      // Derivative information must be available
      casadi_assert(has_derivative(), "Derivatives cannot be calculated for " + name_);
      return wrap().reverse(nadj);
    }
    // Retrieve/generate cached
    Function f;
    string fname = "adj" + str(nadj) + "_" + name_;
    if (!incache(fname, f)) {
      casadi_int i;
      // Names of inputs
      std::vector<std::string> inames;
      for (i=0; i<n_in_; ++i) inames.push_back(name_in_[i]);
      for (i=0; i<n_out_; ++i) inames.push_back("out_" + name_out_[i]);
      for (i=0; i<n_out_; ++i) inames.push_back("adj_" + name_out_[i]);
      // Names of outputs
      std::vector<std::string> onames;
      for (casadi_int i=0; i<n_in_; ++i) onames.push_back("adj_" + name_in_[i]);
      // Options
      Dict opts;
      opts["max_num_dir"] = max_num_dir_;
      opts["derivative_of"] = self();
      // Generate derivative function
      casadi_assert_dev(enable_reverse_);
      f = get_reverse(nadj, fname, inames, onames, opts);
      // Consistency check for inputs
      casadi_assert_dev(f.n_in()==n_in_ + n_out_ + n_out_);
      casadi_int ind=0;
      for (i=0; i<n_in_; ++i) f.assert_size_in(ind++, size1_in(i), size2_in(i));
      for (i=0; i<n_out_; ++i) f.assert_size_in(ind++, size1_out(i), size2_out(i));
      for (i=0; i<n_out_; ++i) f.assert_size_in(ind++, size1_out(i), nadj*size2_out(i));
      // Consistency check for outputs
      casadi_assert_dev(f.n_out()==n_in_);
      for (i=0; i<n_in_; ++i) f.assert_size_out(i, size1_in(i), nadj*size2_in(i));
      // Save to cache
      tocache(f);
    }
    return f;
  }

  Function FunctionInternal::
  get_forward(casadi_int nfwd, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    casadi_error("'get_forward' not defined for " + class_name());
  }

  Function FunctionInternal::
  get_reverse(casadi_int nadj, const std::string& name,
              const std::vector<std::string>& inames,
              const std::vector<std::string>& onames,
              const Dict& opts) const {
    casadi_error("'get_reverse' not defined for " + class_name());
  }

  void FunctionInternal::export_code(const std::string& lang, std::ostream &stream,
      const Dict& options) const {
    casadi_error("'export_code' not defined for " + class_name());
  }

  void FunctionInternal::serialize(std::ostream &stream) const {
    casadi_error("'serialize' not defined for " + class_name());
  }

  void assert_read(std::istream &stream, const std::string& s) {
    casadi_int n = s.size();
    char c;
    std::stringstream ss;
    for (casadi_int i=0;i<n;++i) {
      stream >> c;
      ss << c;
    }
    casadi_assert_dev(s==ss.str());
  }

  void FunctionInternal::serialize_header(std::ostream &stream) const {
    stream << "n" << name_.size() << name_;
    stream << "i" << n_in_ << ":";
    for (int i=0;i<n_in_;++i) stream << name_in_[i].size() << name_in_[i];
    for (int i=0;i<n_in_;++i) sparsity_in_[i].serialize(stream);
    stream << "o" << n_out_ << ":";
    for (int i=0;i<n_out_;++i) stream << name_out_[i].size() << name_out_[i];
    for (int i=0;i<n_out_;++i) sparsity_out_[i].serialize(stream);
    stream << "w" << sz_w();
    stream << "iw" << sz_iw();
  }

  void FunctionInternal::deserialize_header(std::istream& stream,
      std::string& name,
      std::vector<Sparsity>& sp_in, std::vector<Sparsity>& sp_out,
      std::vector<std::string>& names_in, std::vector<std::string>& names_out,
      casadi_int& sz_w, casadi_int& sz_iw) {
    assert_read(stream, "n");
    int name_size; stream >> name_size;
    char c;
    for (int i=0;i<name_size;++i) {
      stream >> c;
      name += c;
    }
    assert_read(stream, "i");
    casadi_int n_in;
    stream >> n_in; stream >> c;
    for (int i=0;i<n_in;++i) {
      std::string name;
      int n; stream >> n;
      for (int k=0;k<n;++k) {
        stream >> c;
        name+= c;
      }
      names_in.push_back(name);
    }
    for (int i=0;i<n_in;++i) sp_in.push_back(Sparsity::deserialize(stream));

    assert_read(stream, "o");
    casadi_int n_out;
    stream >> n_out; stream >> c;
    for (int i=0;i<n_out;++i) {
      std::string name;
      int n; stream >> n;
      for (int k=0;k<n;++k) {
        stream >> c;
        name+= c;
      }
      names_out.push_back(name);
    }
    for (int i=0;i<n_out;++i) sp_out.push_back(Sparsity::deserialize(stream));
    assert_read(stream, "w");
    stream >> sz_w;
    assert_read(stream, "iw");
    stream >> sz_iw;
  }

  casadi_int FunctionInternal::nnz_in() const {
    casadi_int ret=0;
    for (casadi_int iind=0; iind<n_in_; ++iind) ret += nnz_in(iind);
    return ret;
  }

  casadi_int FunctionInternal::nnz_out() const {
    casadi_int ret=0;
    for (casadi_int oind=0; oind<n_out_; ++oind) ret += nnz_out(oind);
    return ret;
  }

  casadi_int FunctionInternal::numel_in() const {
    casadi_int ret=0;
    for (casadi_int iind=0; iind<n_in_; ++iind) ret += numel_in(iind);
    return ret;
  }

  casadi_int FunctionInternal::numel_out() const {
    casadi_int ret=0;
    for (casadi_int oind=0; oind<n_out_; ++oind) ret += numel_out(oind);
    return ret;
  }

  void FunctionInternal::eval_mx(const MXVector& arg, MXVector& res,
                                 bool always_inline, bool never_inline) const {
    // The code below creates a call node, to inline, wrap in an MXFunction
    if (always_inline) {
      casadi_assert(!never_inline, "Inconsistent options");
      return wrap().call(arg, res, true);
    }

    // Create a call-node
    res = Call::create(self(), arg);
  }

  Function FunctionInternal::jac() const {
    // Used wrapped function if jacobian not available
    if (!has_jac()) {
      // Derivative information must be available
      casadi_assert(has_derivative(),
                    "Derivatives cannot be calculated for " + name_);
      return wrap().jac();
    }
    // Retrieve/generate cached
    Function f;
    string fname = "JAC_" + name_;
    if (!incache(fname, f)) {
      // Names of inputs
      std::vector<std::string> inames = name_in_;
      inames.insert(inames.end(), name_out_.begin(), name_out_.end());

      // Names of outputs
      std::vector<std::string> onames;
      onames.reserve(n_in_*n_out_);
      for (size_t oind=0; oind<n_out_; ++oind) {
        for (size_t iind=0; iind<n_in_; ++iind) {
          onames.push_back("D" + name_out_[oind] + "D" + name_in_[iind]);
        }
      }

      // Options
      Dict opts;
      opts["derivative_of"] = self();

      // Generate derivative function
      casadi_assert_dev(enable_jacobian_);
      f = get_jac(fname, inames, onames, opts);

      // Consistency check
      casadi_assert(f.n_in()==inames.size(),
                    "Return function has wrong number of inputs");
      casadi_assert(f.n_out()==onames.size(),
                    "Return function has wrong number of outputs");
      tocache(f);
    }
    return f;
  }

  Function FunctionInternal::jacobian() const {
    // Used wrapped function if jacobian not available
    if (!has_jacobian()) {
      // Derivative information must be available
      casadi_assert(has_derivative(),
                            "Derivatives cannot be calculated for " + name_);
      return wrap().jacobian();
    }

    // Quick return if cached
    if (jacobian_.alive()) {
      return shared_cast<Function>(jacobian_.shared());
    }

    // Give it a suitable name
    string name = "jac_" + name_;

    // Names of inputs
    std::vector<std::string> inames;
    for (casadi_int i=0; i<n_in_; ++i) inames.push_back(name_in_[i]);
    for (casadi_int i=0; i<n_out_; ++i) inames.push_back("out_" + name_out_[i]);

    // Names of outputs
    std::vector<std::string> onames = {"jac"};

    // Options
    Dict opts;
    opts["derivative_of"] = self();

    // Generate derivative function
    casadi_assert_dev(enable_jacobian_);
    Function ret = get_jacobian(name, inames, onames, opts);

    // Consistency check
    casadi_assert_dev(ret.n_in()==n_in_ + n_out_);
    casadi_assert_dev(ret.n_out()==1);

    // Cache it for reuse and return
    jacobian_ = ret;
    return ret;
  }

  Function FunctionInternal::
  get_jacobian(const std::string& name,
               const std::vector<std::string>& inames,
               const std::vector<std::string>& onames,
               const Dict& opts) const {
    casadi_error("'get_jacobian' not defined for " + class_name());
  }

  Function FunctionInternal::
  get_jac(const std::string& name,
               const std::vector<std::string>& inames,
               const std::vector<std::string>& onames,
               const Dict& opts) const {
    casadi_error("'get_jac' not defined for " + class_name());
  }

  Sparsity FunctionInternal::get_jacobian_sparsity() const {
    return wrap()->get_jacobian_sparsity();
  }

  void FunctionInternal::codegen(CodeGenerator& g, const std::string& fname) const {
    // Define function
    g << "/* " << definition() << " */\n";
    g << "static " << signature(fname) << " {\n";

    // Reset local variables, flush buffer
    g.flush(g.body);
    g.local_variables_.clear();
    g.local_default_.clear();

    // Generate function body (to buffer)
    codegen_body(g);

    // Order local variables
    std::map<string, set<pair<string, string>>> local_variables_by_type;
    for (auto&& e : g.local_variables_) {
      local_variables_by_type[e.second.first].insert(make_pair(e.first, e.second.second));
    }

    // Codegen local variables
    for (auto&& e : local_variables_by_type) {
      g.body << "  " << e.first;
      for (auto it=e.second.begin(); it!=e.second.end(); ++it) {
        g.body << (it==e.second.begin() ? " " : ", ") << it->second << it->first;
        // Insert definition, if any
        auto k=g.local_default_.find(it->first);
        if (k!=g.local_default_.end()) g.body << "=" << k->second;
      }
      g.body << ";\n";
    }

    // Finalize the function
    g << "return 0;\n";
    g << "}\n\n";

    // Flush to function body
    g.flush(g.body);
  }

  std::string FunctionInternal::signature(const std::string& fname) const {
    return "int " + fname + "(const casadi_real** arg, casadi_real** res, "
                            "casadi_int* iw, casadi_real* w, void* mem)";
  }

  void FunctionInternal::codegen_sparsities(CodeGenerator& g) const {
    g.add_io_sparsities(name_, sparsity_in_, sparsity_out_);
  }

  void FunctionInternal::codegen_meta(CodeGenerator& g) const {
    // Reference counter routines
    g << g.declare("void " + name_ + "_incref(void)") << " {\n";
    codegen_incref(g);
    g << "}\n\n"
      << g.declare("void " + name_ + "_decref(void)") << " {\n";
    codegen_decref(g);
    g << "}\n\n";

    // Number of inputs and outptus
    g << g.declare("casadi_int " + name_ + "_n_in(void)")
      << " { return " << n_in_ << ";}\n\n"
      << g.declare("casadi_int " + name_ + "_n_out(void)")
      << " { return " << n_out_ << ";}\n\n";

    // Input names
    g << g.declare("const char* " + name_ + "_name_in(casadi_int i)") << "{\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<n_in_; ++i) {
      g << "case " << i << ": return \"" << name_in_[i] << "\";\n";
    }
    g << "default: return 0;\n}\n"
      << "}\n\n";

    // Output names
    g << g.declare("const char* " + name_ + "_name_out(casadi_int i)") << "{\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<n_out_; ++i) {
      g << "case " << i << ": return \"" << name_out_[i] << "\";\n";
    }
    g << "default: return 0;\n}\n"
      << "}\n\n";

    // Codegen sparsities
    codegen_sparsities(g);

    // Determine work vector size
    casadi_int sz_w_codegen = sz_w();
    if (is_a("SXFunction", true) && !g.avoid_stack()) sz_w_codegen = 0;

    // Function that returns work vector lengths
    g << g.declare(
        "int " + name_ + "_work(casadi_int *sz_arg, casadi_int* sz_res, "
        "casadi_int *sz_iw, casadi_int *sz_w)")
      << " {\n"
      << "if (sz_arg) *sz_arg = " << sz_arg() << ";\n"
      << "if (sz_res) *sz_res = " << sz_res() << ";\n"
      << "if (sz_iw) *sz_iw = " << sz_iw() << ";\n"
      << "if (sz_w) *sz_w = " << sz_w_codegen << ";\n"
      << "return 0;\n"
      << "}\n\n";

    // Generate mex gateway for the function
    if (g.mex) {
      // Begin conditional compilation
      g << "#ifdef MATLAB_MEX_FILE\n";

      // Declare wrapper
      g << "void mex_" << name_
        << "(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {\n"
        << "casadi_int i, j;\n";

      // Work vectors, including input and output buffers
      casadi_int i_nnz = nnz_in(), o_nnz = nnz_out();
      size_t sz_w = this->sz_w();
      for (casadi_int i=0; i<n_in_; ++i) {
        const Sparsity& s = sparsity_in_[i];
        sz_w = max(sz_w, static_cast<size_t>(s.size1())); // To be able to copy a column
        sz_w = max(sz_w, static_cast<size_t>(s.size2())); // To be able to copy a row
      }
      sz_w += i_nnz + o_nnz;
      g << g.array("casadi_real", "w", sz_w);
      g << g.array("casadi_int", "iw", sz_iw());
      string fw = "w+" + str(i_nnz + o_nnz);

      // Copy inputs to buffers
      casadi_int offset=0;
      g << g.array("const casadi_real*", "arg", n_in_, "{0}");

      // Allocate output buffers
      g << "casadi_real* res[" << n_out_ << "] = {0};\n";

      // Check arguments
      g << "if (argc>" << n_in_ << ") mexErrMsgIdAndTxt(\"Casadi:RuntimeError\","
        << "\"Evaluation of \\\"" << name_ << "\\\" failed. Too many input arguments "
        << "(%d, max " << n_in_ << ")\", argc);\n";

      g << "if (resc>" << n_out_ << ") mexErrMsgIdAndTxt(\"Casadi:RuntimeError\","
        << "\"Evaluation of \\\"" << name_ << "\\\" failed. "
        << "Too many output arguments (%d, max " << n_out_ << ")\", resc);\n";

      for (casadi_int i=0; i<n_in_; ++i) {
        std::string p = "argv[" + str(i) + "]";
        g << "if (--argc>=0) arg[" << i << "] = "
          << g.from_mex(p, "w", offset, sparsity_in_[i], fw) << "\n";
        offset += nnz_in(i);
      }

      for (casadi_int i=0; i<n_out_; ++i) {
        if (i==0) {
          // if i==0, always store output (possibly ans output)
          g << "--resc;\n";
        } else {
          // Store output, if it exists
          g << "if (--resc>=0) ";
        }
        // Create and get pointer
        g << "res[" << i << "] = w+" << str(offset) << ";\n";
        offset += nnz_out(i);
      }

      // Call the function
      g << "i = " << name_ << "(arg, res, iw, " << fw << ", 0);\n"
        << "if (i) mexErrMsgIdAndTxt(\"Casadi:RuntimeError\",\"Evaluation of \\\"" << name_
        << "\\\" failed.\");\n";

      // Save results
      for (casadi_int i=0; i<n_out_; ++i) {
        string res_i = "res[" + str(i) + "]";
        g << "if (" << res_i << ") resv[" << i << "] = "
          << g.to_mex(sparsity_out_[i], res_i) << "\n";
      }

      // End conditional compilation and function
      g << "}\n"
        << "#endif\n\n";
    }

    if (g.main) {
      // Declare wrapper
      g << "casadi_int main_" << name_ << "(casadi_int argc, char* argv[]) {\n";

      // Work vectors and input and output buffers
      size_t nr = sz_w() + nnz_in() + nnz_out();
      g << g.array("casadi_int", "iw", sz_iw())
        << g.array("casadi_real", "w", nr);

      // Input buffers
      g << "const casadi_real* arg[" << sz_arg() << "] = {";
      casadi_int off=0;
      for (casadi_int i=0; i<n_in_; ++i) {
        if (i!=0) g << ", ";
        g << "w+" << off;
        off += nnz_in(i);
      }
      g << "};\n";

      // Output buffers
      g << "casadi_real* res[" << sz_res() << "] = {";
      for (casadi_int i=0; i<n_out_; ++i) {
        if (i!=0) g << ", ";
        g << "w+" << off;
        off += nnz_out(i);
      }
      g << "};\n";

      // TODO(@jaeandersson): Read inputs from file. For now; read from stdin
      g << "casadi_int j;\n"
        << "casadi_real* a = w;\n"
        << "for (j=0; j<" << nnz_in() << "; ++j) "
        << "scanf(\"%lf\", a++);\n";

      // Call the function
      g << "casadi_int flag = " << name_ << "(arg, res, iw, w+" << off << ", 0);\n"
        << "if (flag) return flag;\n";

      // TODO(@jaeandersson): Write outputs to file. For now: print to stdout
      g << "const casadi_real* r = w+" << nnz_in() << ";\n"
        << "for (j=0; j<" << nnz_out() << "; ++j) "
        << g.printf("%g ", "*r++") << "\n";

      // End with newline
      g << g.printf("\\n") << "\n";

      // Finalize function
      g << "return 0;\n"
        << "}\n\n";
    }

    if (g.with_mem) {
      // Allocate memory
      g << g.declare("casadi_functions* " + name_ + "_functions(void)") << " {\n"
        << "static casadi_functions fun = {\n"
        << name_ << "_incref,\n"
        << name_ << "_decref,\n"
        << name_ << "_n_in,\n"
        << name_ << "_n_out,\n"
        << name_ << "_name_in,\n"
        << name_ << "_name_out,\n"
        << name_ << "_sparsity_in,\n"
        << name_ << "_sparsity_out,\n"
        << name_ << "_work,\n"
        << name_ << "\n"
        << "};\n"
        << "return &fun;\n"
        << "}\n";
    }
    // Flush
    g.flush(g.body);
  }

  std::string FunctionInternal::codegen_name(const CodeGenerator& g) const {
    // Get the index of the function
    for (auto&& e : g.added_functions_) {
      if (e.f.get()==this) return e.codegen_name;
    }
    casadi_error("Function '" + name_ + "' not found");
  }

  void FunctionInternal::codegen_declarations(CodeGenerator& g) const {
    // Nothing to declare
  }

  void FunctionInternal::codegen_body(CodeGenerator& g) const {
    casadi_warning("The function \"" + name_ + "\", which is of type \""
                   + class_name() + "\" cannot be code generated. The generation "
                   "will proceed, but compilation of the code will not be possible.");
    g << "#error Code generation not supported for " << class_name() << "\n";
  }

  std::string FunctionInternal::
  generate_dependencies(const std::string& fname, const Dict& opts) const {
    casadi_error("'generate_dependencies' not defined for " + class_name());
  }

  int FunctionInternal::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    // Loop over outputs
    for (casadi_int oind=0; oind<n_out_; ++oind) {
      // Skip if nothing to assign
      if (res[oind]==nullptr || nnz_out(oind)==0) continue;

      // Clear result
      casadi_fill(res[oind], nnz_out(oind), bvec_t(0));

      // Loop over inputs
      for (casadi_int iind=0; iind<n_in_; ++iind) {
        // Skip if no seeds
        if (arg[iind]==nullptr || nnz_in(iind)==0) continue;

        // Get the sparsity of the Jacobian block
        Sparsity sp = sparsity_jac(iind, oind, true, false);
        if (sp.is_null() || sp.nnz() == 0) continue; // Skip if zero

        // Carry out the sparse matrix-vector multiplication
        casadi_int d1 = sp.size2();
        const casadi_int *colind = sp.colind(), *row = sp.row();
        for (casadi_int cc=0; cc<d1; ++cc) {
          for (casadi_int el = colind[cc]; el < colind[cc+1]; ++el) {
            res[oind][row[el]] |= arg[iind][cc];
          }
        }
      }
    }
    return 0;
  }

  int FunctionInternal::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
    // Loop over outputs
    for (casadi_int oind=0; oind<n_out_; ++oind) {
      // Skip if nothing to assign
      if (res[oind]==nullptr || nnz_out(oind)==0) continue;

      // Loop over inputs
      for (casadi_int iind=0; iind<n_in_; ++iind) {
        // Skip if no seeds
        if (arg[iind]==nullptr || nnz_in(iind)==0) continue;

        // Get the sparsity of the Jacobian block
        Sparsity sp = sparsity_jac(iind, oind, true, false);
        if (sp.is_null() || sp.nnz() == 0) continue; // Skip if zero

        // Carry out the sparse matrix-vector multiplication
        casadi_int d1 = sp.size2();
        const casadi_int *colind = sp.colind(), *row = sp.row();
        for (casadi_int cc=0; cc<d1; ++cc) {
          for (casadi_int el = colind[cc]; el < colind[cc+1]; ++el) {
            arg[iind][cc] |= res[oind][row[el]];
          }
        }
      }

      // Clear seeds
      casadi_fill(res[oind], nnz_out(oind), bvec_t(0));
    }
    return 0;
  }

  void FunctionInternal::sz_work(size_t& sz_arg, size_t& sz_res,
                                 size_t& sz_iw, size_t& sz_w) const {
    sz_arg = this->sz_arg();
    sz_res = this->sz_res();
    sz_iw = this->sz_iw();
    sz_w = this->sz_w();
  }

  void FunctionInternal::alloc_arg(size_t sz_arg, bool persistent) {
    if (persistent) {
      sz_arg_per_ += sz_arg;
    } else {
      sz_arg_tmp_ = max(sz_arg_tmp_, sz_arg);
    }
  }

  void FunctionInternal::alloc_res(size_t sz_res, bool persistent) {
    if (persistent) {
      sz_res_per_ += sz_res;
    } else {
      sz_res_tmp_ = max(sz_res_tmp_, sz_res);
    }
  }

  void FunctionInternal::alloc_iw(size_t sz_iw, bool persistent) {
    if (persistent) {
      sz_iw_per_ += sz_iw;
    } else {
      sz_iw_tmp_ = max(sz_iw_tmp_, sz_iw);
    }
  }

  void FunctionInternal::alloc_w(size_t sz_w, bool persistent) {
    if (persistent) {
      sz_w_per_ += sz_w;
    } else {
      sz_w_tmp_ = max(sz_w_tmp_, sz_w);
    }
  }

  void FunctionInternal::alloc(const Function& f, bool persistent) {
    if (f.is_null()) return;
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f.sz_work(sz_arg, sz_res, sz_iw, sz_w);
    alloc_arg(sz_arg, persistent);
    alloc_res(sz_res, persistent);
    alloc_iw(sz_iw, persistent);
    alloc_w(sz_w, persistent);
  }

  bool FunctionInternal::has_derivative() const {
    return enable_forward_ || enable_reverse_ || enable_jacobian_ || enable_fd_;
  }

  bool FunctionInternal::fwdViaJac(casadi_int nfwd) const {
    if (!enable_forward_ && !enable_fd_) return true;
    if (jac_penalty_==-1) return false;

    // Heuristic 1: Jac calculated via forward mode likely cheaper
    if (jac_penalty_*static_cast<double>(nnz_in())<nfwd) return true;

    // Heuristic 2: Jac calculated via reverse mode likely cheaper
    double w = ad_weight();
    if (enable_reverse_ &&
        jac_penalty_*(1-w)*static_cast<double>(nnz_out())<w*static_cast<double>(nfwd))
      return true;

    return false;
  }

  bool FunctionInternal::adjViaJac(casadi_int nadj) const {
    if (!enable_reverse_) return true;
    if (jac_penalty_==-1) return false;

    // Heuristic 1: Jac calculated via reverse mode likely cheaper
    if (jac_penalty_*static_cast<double>(nnz_out())<nadj) return true;

    // Heuristic 2: Jac calculated via forward mode likely cheaper
    double w = ad_weight();
    if ((enable_forward_ || enable_fd_) &&
        jac_penalty_*w*static_cast<double>(nnz_in())<(1-w)*static_cast<double>(nadj))
      return true;

    return false;
  }

  Dict FunctionInternal::info() const {
    return Dict();
  }

  void FunctionInternal::
  call_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
             const std::vector<std::vector<MX> >& fseed,
             std::vector<std::vector<MX> >& fsens,
             bool always_inline, bool never_inline) const {
    casadi_assert(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert(!always_inline, "Class " + class_name() +
                          " cannot be inlined in an MX expression");

    // Derivative information must be available
    casadi_assert(has_derivative(),
                          "Derivatives cannot be calculated for " + name_);

    // Number of directional derivatives
    casadi_int nfwd = fseed.size();
    fsens.resize(nfwd);

    // Quick return if no seeds
    if (nfwd==0) return;

    // Check if seeds need to have dimensions corrected
    casadi_int npar = 1;
    for (auto&& r : fseed) {
      if (!matching_arg(r, npar)) {
        return FunctionInternal::call_forward(arg, res, replace_fseed(fseed, npar),
                                            fsens, always_inline, never_inline);
      }
    }

    // Calculating full Jacobian and then multiplying
    if (fwdViaJac(nfwd)) {
      // Join forward seeds
      vector<MX> v(nfwd);
      for (casadi_int d=0; d<nfwd; ++d) {
        v[d] = veccat(fseed[d]);
      }

      // Multiply the Jacobian from the right
      vector<MX> darg = arg;
      darg.insert(darg.end(), res.begin(), res.end());
      MX J = jacobian()(darg).at(0);
      v = horzsplit(mtimes(J, horzcat(v)));

      // Vertical offsets
      vector<casadi_int> offset(n_out_+1, 0);
      for (casadi_int i=0; i<n_out_; ++i) {
        offset[i+1] = offset[i]+numel_out(i);
      }

      // Collect forward sensitivities
      for (casadi_int d=0; d<nfwd; ++d) {
        fsens[d] = vertsplit(v[d], offset);
        for (casadi_int i=0; i<n_out_; ++i) {
          fsens[d][i] = reshape(fsens[d][i], size_out(i));
        }
      }

    } else {
      // Evaluate in batches
      casadi_assert_dev(enable_forward_ || enable_fd_);
      casadi_int max_nfwd = max_num_dir_;
      if (!enable_fd_) {
        while (!has_forward(max_nfwd)) max_nfwd/=2;
      }
      casadi_int offset = 0;
      while (offset<nfwd) {
        // Number of derivatives, in this batch
        casadi_int nfwd_batch = min(nfwd-offset, max_nfwd);

        // All inputs and seeds
        vector<MX> darg;
        darg.reserve(n_in_ + n_out_ + n_in_);
        darg.insert(darg.end(), arg.begin(), arg.end());
        darg.insert(darg.end(), res.begin(), res.end());
        vector<MX> v(nfwd_batch);
        for (casadi_int i=0; i<n_in_; ++i) {
          for (casadi_int d=0; d<nfwd_batch; ++d) v[d] = fseed[offset+d][i];
          darg.push_back(horzcat(v));
        }

        // Create the evaluation node
        Function dfcn = self().forward(nfwd_batch);
        vector<MX> x = Call::create(dfcn, darg);

        casadi_assert_dev(x.size()==n_out_);

        // Retrieve sensitivities
        for (casadi_int d=0; d<nfwd_batch; ++d) fsens[offset+d].resize(n_out_);
        for (casadi_int i=0; i<n_out_; ++i) {
          if (size2_out(i)>0) {
            v = horzsplit(x[i], size2_out(i));
            casadi_assert_dev(v.size()==nfwd_batch);
          } else {
            v = vector<MX>(nfwd_batch, MX(size_out(i)));
          }
          for (casadi_int d=0; d<nfwd_batch; ++d) fsens[offset+d][i] = v[d];
        }

        // Update offset
        offset += nfwd_batch;
      }
    }
  }

  void FunctionInternal::
  call_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
             const std::vector<std::vector<MX> >& aseed,
             std::vector<std::vector<MX> >& asens,
             bool always_inline, bool never_inline) const {
    casadi_assert(!(always_inline && never_inline), "Inconsistent options");
    casadi_assert(!always_inline, "Class " + class_name() +
                          " cannot be inlined in an MX expression");

    // Derivative information must be available
    casadi_assert(has_derivative(),
                          "Derivatives cannot be calculated for " + name_);

    // Number of directional derivatives
    casadi_int nadj = aseed.size();
    asens.resize(nadj);

    // Quick return if no seeds
    if (nadj==0) return;

    // Check if seeds need to have dimensions corrected
    casadi_int npar = 1;
    for (auto&& r : aseed) {
      if (!matching_res(r, npar)) {
        return FunctionInternal::call_reverse(arg, res, replace_aseed(aseed, npar),
                                            asens, always_inline, never_inline);
      }
    }

    // Calculating full Jacobian and then multiplying likely cheaper
    if (adjViaJac(nadj)) {
      // Join adjoint seeds
      vector<MX> v(nadj);
      for (casadi_int d=0; d<nadj; ++d) {
        v[d] = veccat(aseed[d]);
      }

      // Multiply the transposed Jacobian from the right
      vector<MX> darg = arg;
      darg.insert(darg.end(), res.begin(), res.end());
      MX J = jacobian()(darg).at(0);
      v = horzsplit(mtimes(J.T(), horzcat(v)));

      // Vertical offsets
      vector<casadi_int> offset(n_in_+1, 0);
      for (casadi_int i=0; i<n_in_; ++i) {
        offset[i+1] = offset[i]+numel_in(i);
      }

      // Collect adjoint sensitivities
      for (casadi_int d=0; d<nadj; ++d) {
        asens[d].resize(n_in_);
        vector<MX> a = vertsplit(v[d], offset);
        for (casadi_int i=0; i<n_in_; ++i) {
          if (asens[d][i].is_empty(true)) {
            asens[d][i] = reshape(a[i], size_in(i));
          } else {
            asens[d][i] += reshape(a[i], size_in(i));
          }
        }
      }
    } else {
      // Evaluate in batches
      casadi_assert_dev(enable_reverse_);
      casadi_int max_nadj = max_num_dir_;

      while (!has_reverse(max_nadj)) max_nadj/=2;
      casadi_int offset = 0;
      while (offset<nadj) {
        // Number of derivatives, in this batch
        casadi_int nadj_batch = min(nadj-offset, max_nadj);

        // All inputs and seeds
        vector<MX> darg;
        darg.reserve(n_in_ + n_out_ + n_out_);
        darg.insert(darg.end(), arg.begin(), arg.end());
        darg.insert(darg.end(), res.begin(), res.end());
        vector<MX> v(nadj_batch);
        for (casadi_int i=0; i<n_out_; ++i) {
          for (casadi_int d=0; d<nadj_batch; ++d) v[d] = aseed[offset+d][i];
          darg.push_back(horzcat(v));
        }

        // Create the evaluation node
        Function dfcn = self().reverse(nadj_batch);
        vector<MX> x = Call::create(dfcn, darg);
        casadi_assert_dev(x.size()==n_in_);

        // Retrieve sensitivities
        for (casadi_int d=0; d<nadj_batch; ++d) asens[offset+d].resize(n_in_);
        for (casadi_int i=0; i<n_in_; ++i) {
          if (size2_in(i)>0) {
            v = horzsplit(x[i], size2_in(i));
            casadi_assert_dev(v.size()==nadj_batch);
          } else {
            v = vector<MX>(nadj_batch, MX(size_in(i)));
          }
          for (casadi_int d=0; d<nadj_batch; ++d) {
            if (asens[offset+d][i].is_empty(true)) {
              asens[offset+d][i] = v[d];
            } else {
              asens[offset+d][i] += v[d];
            }
          }
        }
        // Update offset
        offset += nadj_batch;
      }
    }
  }

  void FunctionInternal::
  call_forward(const std::vector<SX>& arg, const std::vector<SX>& res,
             const std::vector<std::vector<SX> >& fseed,
             std::vector<std::vector<SX> >& fsens,
             bool always_inline, bool never_inline) const {
    casadi_assert(!(always_inline && never_inline), "Inconsistent options");
    if (fseed.empty()) { // Quick return if no seeds
      fsens.clear();
      return;
    }
    casadi_error("'forward' (SX) not defined for " + class_name());
  }

  void FunctionInternal::
  call_reverse(const std::vector<SX>& arg, const std::vector<SX>& res,
             const std::vector<std::vector<SX> >& aseed,
             std::vector<std::vector<SX> >& asens,
             bool always_inline, bool never_inline) const {
    casadi_assert(!(always_inline && never_inline), "Inconsistent options");
    if (aseed.empty()) { // Quick return if no seeds
      asens.clear();
      return;
    }
    casadi_error("'reverse' (SX) not defined for " + class_name());
  }

  double FunctionInternal::ad_weight() const {
    // If reverse mode derivatives unavailable, use forward
    if (!enable_reverse_) return 0;

    // If forward mode derivatives unavailable, use reverse
    if (!enable_forward_ && !enable_fd_) return 1;

    // Use the (potentially user set) option
    return ad_weight_;
  }

  double FunctionInternal::sp_weight() const {
    // If reverse mode propagation unavailable, use forward
    if (!has_sprev()) return 0;

    // If forward mode propagation unavailable, use reverse
    if (!has_spfwd()) return 1;

    // Use the (potentially user set) option
    return ad_weight_sp_;
  }

  const SX FunctionInternal::sx_in(casadi_int ind) const {
    return SX::sym("x_" + str(ind), sparsity_in(ind));
  }

  const SX FunctionInternal::sx_out(casadi_int ind) const {
    return SX::sym("r_" + str(ind), sparsity_out(ind));
  }

  const std::vector<SX> FunctionInternal::sx_in() const {
    vector<SX> ret(n_in_);
    for (casadi_int i=0; i<ret.size(); ++i) {
      ret[i] = sx_in(i);
    }
    return ret;
  }

  const std::vector<SX> FunctionInternal::sx_out() const {
    vector<SX> ret(n_out_);
    for (casadi_int i=0; i<ret.size(); ++i) {
      ret[i] = sx_out(i);
    }
    return ret;
  }

  const MX FunctionInternal::mx_in(casadi_int ind) const {
    return MX::sym("x_" + str(ind), sparsity_in(ind));
  }

  const MX FunctionInternal::mx_out(casadi_int ind) const {
    return MX::sym("r_" + str(ind), sparsity_out(ind));
  }

  const std::vector<MX> FunctionInternal::mx_in() const {
    vector<MX> ret(n_in_);
    for (casadi_int i=0; i<ret.size(); ++i) {
      ret[i] = mx_in(i);
    }
    return ret;
  }

  const std::vector<MX> FunctionInternal::mx_out() const {
    vector<MX> ret(n_out_);
    for (casadi_int i=0; i<ret.size(); ++i) {
      ret[i] = mx_out(i);
    }
    return ret;
  }

  bool FunctionInternal::is_a(const std::string& type, bool recursive) const {
    return type == "FunctionInternal";
  }

  std::vector<MX> FunctionInternal::free_mx() const {
    casadi_error("'free_mx' only defined for 'MXFunction'");
  }

  std::vector<SX> FunctionInternal::free_sx() const {
    casadi_error("'free_sx' only defined for 'SXFunction'");
  }

  void FunctionInternal::generate_lifted(Function& vdef_fcn,
                                         Function& vinit_fcn) const {
    casadi_error("'generate_lifted' only defined for 'MXFunction'");
  }

  casadi_int FunctionInternal::n_instructions() const {
    casadi_error("'n_instructions' not defined for " + class_name());
  }

  casadi_int FunctionInternal::instruction_id(casadi_int k) const {
    casadi_error("'instruction_id' not defined for " + class_name());
  }

  std::vector<casadi_int> FunctionInternal::instruction_input(casadi_int k) const {
    casadi_error("'instruction_input' not defined for " + class_name());
  }

  double FunctionInternal::instruction_constant(casadi_int k) const {
    casadi_error("'instruction_constant' not defined for " + class_name());
  }

  std::vector<casadi_int> FunctionInternal::instruction_output(casadi_int k) const {
    casadi_error("'instruction_output' not defined for " + class_name());
  }

  MX FunctionInternal::instruction_MX(casadi_int k) const {
    casadi_error("'instruction_MX' not defined for " + class_name());
  }

  casadi_int FunctionInternal::n_nodes() const {
    casadi_error("'n_nodes' not defined for " + class_name());
  }

  std::vector<MX>
  FunctionInternal::mapsum_mx(const std::vector<MX > &x,
                              const std::string& parallelization) {
    if (x.empty()) return x;
    // Check number of arguments
    casadi_assert(x.size()==n_in_, "mapsum_mx: Wrong number of arguments");
    // Number of parallel calls
    casadi_int npar = 1;
    // Check/replace arguments
    std::vector<MX> x_mod(x.size());
    for (casadi_int i=0; i<n_in_; ++i) {
      if (check_mat(x[i].sparsity(), sparsity_in_[i], npar)) {
        x_mod[i] = replace_mat(x[i], sparsity_in_[i], npar);
      } else {
        // Mismatching sparsity: The following will throw an error message
        npar = 0;
        check_arg(x, npar);
      }
    }

    casadi_int n = 1;
    for (casadi_int i=0; i<x_mod.size(); ++i) {
      n = max(x_mod[i].size2()/size2_in(i), n);
    }

    vector<casadi_int> reduce_in;
    for (casadi_int i=0; i<x_mod.size(); ++i) {
      if (x_mod[i].size2()/size2_in(i)!=n) {
        reduce_in.push_back(i);
      }
    }

    Function ms = self().map("mapsum", parallelization, n, reduce_in, range(n_out_));

    // Call the internal function
    return ms(x_mod);
  }

  bool FunctionInternal::check_mat(const Sparsity& arg, const Sparsity& inp, casadi_int& npar) {
    // Matching dimensions
    if (arg.size()==inp.size()) return true;
    // Calling with empty matrix - set all to zero
    if (arg.is_empty()) return true;
    // Calling with a scalar - set all
    if (arg.is_scalar()) return true;
    // Vectors that are transposes of each other
    if (arg.is_vector() && inp.size()==make_pair(arg.size2(), arg.size1())) return true;
    // Horizontal repmat
    if (arg.size1()==inp.size1() && arg.size2()>0 && inp.size2()>0
        && inp.size2()%arg.size2()==0) return true;
    // Evaluate with multiple arguments
    if (arg.size1()==inp.size1() && arg.size2()>0 && inp.size2()>0
        && arg.size2()%(npar*inp.size2())==0) {
      npar *= arg.size2()/(npar*inp.size2());
      return true;
    }
    // No match
    return false;
  }

  void FunctionInternal::setup(void* mem, const double** arg, double** res,
                               casadi_int* iw, double* w) const {
    set_work(mem, arg, res, iw, w);
    set_temp(mem, arg, res, iw, w);
  }

  void ProtoFunction::free_mem(void *mem) const {
    casadi_warning("'free_mem' not defined for " + class_name());
  }

  void ProtoFunction::clear_mem() {
    for (auto&& i : mem_) {
      if (i!=nullptr) free_mem(i);
    }
    mem_.clear();
  }

  size_t FunctionInternal::get_n_in() {
    if (!derivative_of_.is_null()) {
      string n = derivative_of_.name();
      if (name_ == "jac_" + n) {
        return derivative_of_.n_in() + derivative_of_.n_out();
      }
    }
    // One by default
    return 1;
  }

  size_t FunctionInternal::get_n_out() {
    if (!derivative_of_.is_null()) {
      string n = derivative_of_.name();
      if (name_ == "jac_" + n) {
        return 1;
      }
    }
    // One by default
    return 1;
  }

  Sparsity FunctionInternal::get_sparsity_in(casadi_int i) {
    if (!derivative_of_.is_null()) {
      string n = derivative_of_.name();
      if (name_ == "jac_" + n) {
        if (i < derivative_of_.n_in()) {
          // Same as nondifferentiated function
          return derivative_of_.sparsity_in(i);
        } else {
          // Dummy output
          return Sparsity(derivative_of_.size_out(i-derivative_of_.n_in()));
        }
      }
    }
    // Scalar by default
    return Sparsity::scalar();
  }

  Sparsity FunctionInternal::get_sparsity_out(casadi_int i) {
    if (!derivative_of_.is_null()) {
      string n = derivative_of_.name();
      if (name_ == "jac_" + n) {
        // Dense Jacobian by default
        return Sparsity::dense(derivative_of_.nnz_out(), derivative_of_.nnz_in());
      }
    }
    // Scalar by default
    return Sparsity::scalar();
  }

  void* ProtoFunction::memory(casadi_int ind) const {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mtx_);
#endif //CASADI_WITH_THREAD
    return mem_.at(ind);
  }

  casadi_int ProtoFunction::checkout() const {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mtx_);
#endif //CASADI_WITH_THREAD
    if (unused_.empty()) {
      // Allocate a new memory object
      void* m = alloc_mem();
      mem_.push_back(m);
      if (init_mem(m)) {
        casadi_error("Failed to create or initialize memory object");
      }
      return mem_.size()-1;
    } else {
      // Use an unused memory object
      casadi_int m = unused_.top();
      unused_.pop();
      return m;
    }
  }

  void ProtoFunction::release(casadi_int mem) const {
#ifdef CASADI_WITH_THREAD
    std::lock_guard<std::mutex> lock(mtx_);
#endif //CASADI_WITH_THREAD
    unused_.push(mem);
  }

  Function FunctionInternal::
  factory(const std::string& name,
          const std::vector<std::string>& s_in,
          const std::vector<std::string>& s_out,
          const Function::AuxOut& aux,
          const Dict& opts) const {
    return wrap().factory(name, s_in, s_out, aux, opts);
  }

  std::vector<std::string> FunctionInternal::get_function() const {
    // No functions
    return std::vector<std::string>();
  }

  const Function& FunctionInternal::get_function(const std::string &name) const {
    casadi_error("'get_function' not defined for " + class_name());
    static Function singleton;
    return singleton;
  }

  vector<bool> FunctionInternal::
  which_depends(const string& s_in, const vector<string>& s_out, casadi_int order, bool tr) const {
    casadi_error("'which_depends' not defined for " + class_name());
    return vector<bool>();
  }

  const Function& FunctionInternal::oracle() const {
    casadi_error("'oracle' not defined for " + class_name());
    static Function singleton;
    return singleton;
  }

  Function FunctionInternal::slice(const std::string& name,
        const std::vector<casadi_int>& order_in,
        const std::vector<casadi_int>& order_out, const Dict& opts) const {
    return wrap().slice(name, order_in, order_out, opts);
  }

  bool FunctionInternal::all_scalar() const {
    // Check inputs
    for (casadi_int i=0; i<n_in_; ++i) {
      if (!sparsity_in_[i].is_scalar()) return false;
    }
    // Check outputs
    for (casadi_int i=0; i<n_out_; ++i) {
      if (!sparsity_out_[i].is_scalar()) return false;
    }
    // All are scalar
    return true;
  }

  void FunctionInternal::set_jac_sparsity(const Sparsity& sp) {
    // Make sure that it's of the right size
    casadi_assert_dev(sp.size1()==numel_out());
    casadi_assert_dev(sp.size2()==numel_in());
    // Split up into the individual patterns
    std::vector<casadi_int> v_offset(n_out_+1, 0);
    for (casadi_int i=0; i<n_out_; ++i) v_offset[i+1] = v_offset[i] + numel_out(i);
    std::vector<casadi_int> h_offset(n_in_+1, 0);
    for (casadi_int i=0; i<n_in_; ++i) h_offset[i+1] = h_offset[i] + numel_in(i);
    vector<vector<Sparsity>> blocks = blocksplit(sp, v_offset, h_offset);
    // Save to jac_sparsity_ and jac_sparsity_compact_
    for (casadi_int oind=0; oind<n_out_; ++oind) {
      vector<casadi_int> row_nz = sparsity_out(oind).find();
      for (casadi_int iind=0; iind<n_in_; ++iind) {
        vector<casadi_int> col_nz = sparsity_in(iind).find();
        const Sparsity& sp = blocks.at(oind).at(iind);
        jac_sparsity_.elem(oind, iind) = sp;
        vector<casadi_int> mapping;
        jac_sparsity_compact_.elem(oind, iind) = sp.sub(row_nz, col_nz, mapping);
      }
    }
  }

  int FunctionInternal::
  eval_fallback(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // As a fallback, redirect to (the less efficient) eval_dm

    // Allocate input matrices
    std::vector<DM> argv(n_in_);
    for (casadi_int i=0; i<n_in_; ++i) {
      argv[i] = DM(sparsity_in_[i]);
      casadi_copy(arg[i], argv[i].nnz(), argv[i].ptr());
    }

    // Try to evaluate using eval_dm
    try {
      std::vector<DM> resv = eval_dm(argv);

      // Check number of outputs
      casadi_assert(resv.size()==n_out_,
        "Expected " + str(n_out_) + " outputs, got " + str(resv.size()) + ".");

      // Get outputs
      for (casadi_int i=0; i<n_out_; ++i) {
        if (resv[i].sparsity()!=sparsity_out_[i]) {
          if (resv[i].size()==size_out(i)) {
            resv[i] = project(resv[i], sparsity_out_[i]);
          } else {
            casadi_error("Shape mismatch for output " + str(i) + ": got " + resv[i].dim() + ", "
                         "expected " + sparsity_out_[i].dim() + ".");
          }
        }
        if (res[i]) casadi_copy(resv[i].ptr(), resv[i].nnz(), res[i]);
      }
    } catch (KeyboardInterruptException&) {
      throw;
    } catch (exception& e) {
      casadi_error("Failed to evaluate 'eval_dm' for " + name_ + ":\n" + e.what());
      return 1;
    }

    // Successful return
    return 0;
  }

  int FunctionInternal::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
#ifdef WITH_EXTRA_CHECKS
    casadi_int call_depth = Function::call_depth_;
    Function::call_depth_ = 0;
#endif // WITH_EXTRA_CHECKS
    try {
      int ret = eval_fallback(arg, res, iw, w, mem);
#ifdef WITH_EXTRA_CHECKS
      Function::call_depth_ = call_depth;
#endif // WITH_EXTRA_CHECKS
      return ret;
    } catch (exception&) {
#ifdef WITH_EXTRA_CHECKS
      Function::call_depth_ = call_depth;
#endif
      throw;
    }
  }

  void FunctionInternal::sprint(char* buf, size_t buf_sz, const char* fmt, ...) const {
    // Variable number of arguments
    va_list args;
    va_start(args, fmt);
    // Print to buffer
    casadi_int n = vsnprintf(buf, buf_sz, fmt, args);
    // Cleanup
    va_end(args);
    // Throw error if failure
    casadi_assert(n>=0 && n<buf_sz, "Print failure while processing '" + string(fmt) + "'");
  }

  void FunctionInternal::print(const char* fmt, ...) const {
    // Variable number of arguments
    va_list args;
    va_start(args, fmt);
    // Static & dynamic buffers
    char buf[256];
    size_t buf_sz = sizeof(buf);
    char* buf_dyn = nullptr;
    // Try to print with a small buffer
    casadi_int n = vsnprintf(buf, buf_sz, fmt, args);
    // Need a larger buffer?
    if (n>static_cast<casadi_int>(buf_sz)) {
      buf_sz = static_cast<size_t>(n+1);
      buf_dyn = new char[buf_sz];
      n = vsnprintf(buf_dyn, buf_sz, fmt, args);
    }
    // Print buffer content
    if (n>=0) uout() << (buf_dyn ? buf_dyn : buf) << std::flush;
    // Cleanup
    if (buf_dyn) delete[] buf_dyn;
    va_end(args);
    // Throw error if failure
    casadi_assert(n>=0, "Print failure while processing '" + string(fmt) + "'");
  }

  void FunctionInternal::
  call_gen(const MXVector& arg, MXVector& res, casadi_int npar,
           bool always_inline, bool never_inline) const {
    if (npar==1) {
      eval_mx(arg, res, always_inline, never_inline);
    } else {
      // Split it up arguments
      std::vector<std::vector<MX>> v(npar, arg);
      std::vector<MX> t;
      for (int i=0; i<n_in_; ++i) {
        if (arg[i].size2()!=size2_in(i)) {
          t = horzsplit(arg[i], size2_in(i));
          casadi_assert_dev(t.size()==npar);
          for (int p=0; p<npar; ++p) v[p][i] = t[p];
        }
      }
      // Unroll the loop
      for (int p=0; p<npar; ++p) {
        eval_mx(v[p], t, always_inline, never_inline);
        v[p] = t;
      }
      // Concatenate results
      t.resize(npar);
      res.resize(n_out_);
      for (int i=0; i<n_out_; ++i) {
        for (int p=0; p<npar; ++p) t[p] = v[p][i];
        res[i] = horzcat(t);
      }
    }
  }


} // namespace casadi
