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


#include "external_impl.hpp"
#include "casadi_misc.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

  Function external(const string& name, const Importer& li,
                    const Dict& opts) {
    return Function::create(new GenericExternal(name, li), opts);
  }

  Function external(const string& name, const Dict& opts) {
    return external(name, "./" + name + SHARED_LIBRARY_SUFFIX, opts);
  }

  Function external(const string& name, const string& bin_name,
                    const Dict& opts) {
    return external(name, Importer(bin_name, "dll"), opts);
  }

  External::External(const std::string& name, const Importer& li)
    : FunctionInternal(name), li_(li) {

    // Increasing/decreasing reference counter
    incref_ = (signal_t)li_.get_function(name_ + "_incref");
    decref_ = (signal_t)li_.get_function(name_ + "_decref");

    // Getting number of inputs and outputs
    get_n_in_ = (getint_t)li_.get_function(name + "_n_in");
    get_n_out_ = (getint_t)li_.get_function(name + "_n_out");

    // Getting names of inputs and outputs
    get_name_in_ = (name_t)li_.get_function(name + "_name_in");
    get_name_out_ = (name_t)li_.get_function(name + "_name_out");

    // Work vector sizes
    work_ = (work_t)li_.get_function(name_ + "_work");

    // Increase reference counter - external function memory initialized at this point
    if (incref_) incref_();
  }

  GenericExternal::GenericExternal(const std::string& name, const Importer& li)
    : External(name, li) {

    // Functions for retrieving sparsities of inputs and outputs
    get_sparsity_in_ = (sparsity_t)li_.get_function(name + "_sparsity_in");
    get_sparsity_out_ = (sparsity_t)li_.get_function(name + "_sparsity_out");

    // Memory allocation functions
    alloc_mem_ = (alloc_mem_t)li_.get_function(name_ + "_alloc_mem");
    init_mem_ = (init_mem_t)li_.get_function(name_ + "_init_mem");
    free_mem_ = (free_mem_t)li_.get_function(name_ + "_free_mem");

    // Function for numerical evaluation
    eval_ = (eval_t)li_.get_function(name_);
  }

  External::~External() {
    if (decref_) decref_();
    clear_mem();
  }

  size_t External::get_n_in() {
    if (get_n_in_) {
      return get_n_in_();
    } else if (li_.has_meta(name_ + "_N_IN")) {
      return li_.meta_int(name_ + "_N_IN");
    } else {
      // Fall back to base class
      return FunctionInternal::get_n_in();
    }
  }

  size_t External::get_n_out() {
    if (get_n_out_) {
      return get_n_out_();
    } else if (li_.has_meta(name_ + "_N_OUT")) {
      return li_.meta_int(name_ + "_N_OUT");
    } else {
      // Fall back to base class
      return FunctionInternal::get_n_out();
    }
  }

  string External::get_name_in(casadi_int i) {
    if (get_name_in_) {
      // Use function pointer
      const char* n = get_name_in_(i);
      casadi_assert(n!=nullptr, "Error querying input name");
      return n;
    } else if (li_.has_meta(name_ + "_NAME_IN", i)) {
      // Read meta
      return li_.meta_string(name_ + "_NAME_IN", i);
    } else {
      // Default name
      return FunctionInternal::get_name_in(i);
    }
  }

  string External::get_name_out(casadi_int i) {
    if (get_name_out_) {
      // Use function pointer
      const char* n = get_name_out_(i);
      casadi_assert(n!=nullptr, "Error querying output name");
      return n;
    } else if (li_.has_meta(name_ + "_NAME_OUT", i)) {
      // Read meta
      return li_.meta_string(name_ + "_NAME_OUT", i);
    } else {
      // Default name
      return FunctionInternal::get_name_out(i);
    }
  }

  Sparsity GenericExternal::get_sparsity_in(casadi_int i) {
    // Use sparsity retrieval function, if present
    if (get_sparsity_in_) {
      return Sparsity::compressed(get_sparsity_in_(i));
    } else if (li_.has_meta(name_ + "_SPARSITY_IN", i)) {
      return Sparsity::compressed(li_.meta_vector<casadi_int>(name_ + "_SPARSITY_IN", i));
    } else {
      // Fall back to base class
      return FunctionInternal::get_sparsity_in(i);
    }
  }

  Sparsity GenericExternal::get_sparsity_out(casadi_int i) {
    // Use sparsity retrieval function, if present
    if (get_sparsity_out_) {
      return Sparsity::compressed(get_sparsity_out_(i));
    } else if (li_.has_meta(name_ + "_SPARSITY_OUT", i)) {
      return Sparsity::compressed(li_.meta_vector<casadi_int>(name_ + "_SPARSITY_OUT", i));
    } else {
      // Fall back to base class
      return FunctionInternal::get_sparsity_out(i);
    }
  }

  void* GenericExternal::alloc_mem() const {
    if (alloc_mem_) {
      return alloc_mem_();
    } else {
      return FunctionInternal::alloc_mem();
    }
  }

  int GenericExternal::init_mem(void* mem) const {
    if (init_mem_) {
      return init_mem_(mem);
    } else {
      return FunctionInternal::init_mem(mem);
    }
  }

  void GenericExternal::free_mem(void *mem) const {
    if (free_mem_) {
      return free_mem_(mem);
    } else {
      return FunctionInternal::free_mem(mem);
    }
  }

  void External::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Reference counting?
    has_refcount_ = li_.has_function(name_ + "_incref");
    casadi_assert(has_refcount_==li_.has_function(name_ + "_decref"),
                          "External functions must provide functions for both increasing "
                          "and decreasing the reference count, or neither.");

    // Allocate work vectors
    casadi_int sz_arg=0, sz_res=0, sz_iw=0, sz_w=0;
    if (work_) {
      casadi_int flag = work_(&sz_arg, &sz_res, &sz_iw, &sz_w);
      casadi_assert(flag==0, "External: \"work\" failed");
    } else if (li_.has_meta(name_ + "_WORK")) {
      vector<casadi_int> v = li_.meta_vector<casadi_int>(name_ + "_WORK");
      casadi_assert_dev(v.size()==4);
      sz_arg = v[0];
      sz_res = v[1];
      sz_iw = v[2];
      sz_w = v[3];
    }

    // Get information about Jacobian sparsity, if any
    sparsity_t jac_sparsity_fcn = (sparsity_t)li_.get_function("jac_" + name_ + "_sparsity_out");
    if (jac_sparsity_fcn) {
      set_jac_sparsity(Sparsity::compressed(jac_sparsity_fcn(0)));
    } else if (li_.has_meta("JAC_" + name_ + "_SPARSITY_OUT", 0)) {
      vector<casadi_int> sp = li_.meta_vector<casadi_int>("jac_" + name_ + "_SPARSITY_OUT", 0);
      set_jac_sparsity(Sparsity::compressed(sp));
    }

    alloc_arg(sz_arg);
    alloc_res(sz_res);
    alloc_iw(sz_iw);
    alloc_w(sz_w);
  }

  void GenericExternal::init(const Dict& opts) {
    // Call recursively
    External::init(opts);
  }

  void External::codegen_declarations(CodeGenerator& g) const {
    if (!li_.inlined(name_)) {
      g.add_external(signature(name_) + ";");
    }
  }

  void External::codegen_body(CodeGenerator& g) const {
    if (li_.inlined(name_)) {
      // Function body is inlined
      g << li_.body(name_) << "\n";
    } else {
      g << "if (" << name_ << "(arg, res, iw, w, mem)) return 1;\n";
    }
  }

  bool External::has_jacobian() const {
    if (FunctionInternal::has_jacobian()) return true;
    return li_.has_function("jac_" + name_);
  }

  Function External
  ::get_jacobian(const std::string& name,
                    const std::vector<std::string>& inames,
                    const std::vector<std::string>& onames,
                    const Dict& opts) const {
    if (has_jacobian()) {
      return external(name, li_, opts);
    } else {
      return FunctionInternal::get_jacobian(name, inames, onames, opts);
    }
  }

  Function External
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Consistency check
    casadi_int n=1;
    while (n<nfwd) n*=2;
    if (n!=nfwd || !has_forward(nfwd)) {
      // Inefficient code to be replaced later
      Function fwd1 = forward(1);
      return fwd1.map(name, "serial", nfwd, range(n_in_+n_out_), std::vector<casadi_int>(), opts);
      //casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  bool External::has_forward(casadi_int nfwd) const {
    return li_.has_function("fwd" + str(nfwd) + "_" + name_);
  }

  Function External
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Consistency check
    casadi_int n=1;
    while (n<nadj) n*=2;
    if (n!=nadj || !has_reverse(nadj)) {
      // Inefficient code to be replaced later
      Function adj1 = reverse(1);
      return adj1.map(name, "serial", nadj, range(n_in_+n_out_), std::vector<casadi_int>(), opts);
      //casadi_error("Internal error: Refactoring needed, cf. #1055");
    }
    return external(name, li_, opts);
  }

  bool External::has_reverse(casadi_int nadj) const {
    return li_.has_function("adj" + str(nadj) + "_" + name_);
  }

  Function External::factory(const std::string& name,
                             const std::vector<std::string>& s_in,
                             const std::vector<std::string>& s_out,
                             const Function::AuxOut& aux,
                             const Dict& opts) const {
    // If not available, call base class function
    if (!li_.has_function(name)) {
      return FunctionInternal::factory(name, s_in, s_out, aux, opts);
    }

    // Retrieve function
    Function ret = external(name, li_, opts);

    // Inputs consistency checks
    casadi_assert(s_in.size() == ret.n_in(),
      "Inconsistent number of inputs. Expected " + str(s_in.size())+ "  "
      "(" + str(s_in) + "), got " + str(ret.n_in()) + ".");
    for (casadi_int i=0; i<s_in.size(); ++i) {
      string s = s_in[i];
      replace(s.begin(), s.end(), ':', '_');
      casadi_assert(s == ret.name_in(i),
        "Inconsistent input name. Expected: " + str(s_in) + ", "
        "got " + ret.name_in(i) + " for input " + str(i));
    }

    // Outputs consistency checks
    casadi_assert(s_out.size() == ret.n_out(),
      "Inconsistent number of outputs. Expected " + str(s_out.size()) + " "
      "(" + str(s_out) + "), got " + str(ret.n_out()) + ".");
    for (casadi_int i=0; i<s_out.size(); ++i) {
      string s = s_out[i];
      replace(s.begin(), s.end(), ':', '_');
      casadi_assert(s == ret.name_out(i),
        "Inconsistent output name. Expected: " + str(s_out) + ", "
        "got " + ret.name_out(i) + " for output " + str(i));
    }

    return ret;
  }

} // namespace casadi
