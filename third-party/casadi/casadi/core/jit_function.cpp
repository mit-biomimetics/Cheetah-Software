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


#include "jit_function.hpp"
#include "casadi_misc.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace casadi {
  using namespace std;

  JitFunction::JitFunction(const std::string& name, const std::string& body,
                      const std::vector<std::string>& name_in,
                      const std::vector<std::string>& name_out,
                      const std::vector<Sparsity>& sparsity_in,
                      const std::vector<Sparsity>& sparsity_out)
                      : FunctionInternal(name), body_(body) {
    // Set sparsity
    sparsity_in_ = sparsity_in;
    sparsity_out_ = sparsity_out;
    name_in_ = name_in;
    name_out_ = name_out;

    // Default options
    jit_ = true; // override default
    buffered_ = true;
    enable_fd_ = true; // override default
  }

  Options JitFunction::options_
  = {{&FunctionInternal::options_},
     {{"buffered",
      {OT_BOOL,
        "Buffer the calls, user does not need to "}},
       {"jac",
      {OT_STRING,
        "Function body for Jacobian"}},
      {"hess",
       {OT_STRING,
        "Function body for Hessian"}}
     }
  };

  void JitFunction::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="buffered") {
        buffered_ = op.second;
      } else if (op.first=="jac") {
        jac_body_ = op.second.to_string();
      } else if (op.first=="hess") {
        hess_body_ = op.second.to_string();
      }
    }

    // Arrays for holding inputs and outputs
    if (buffered_) {
      alloc_w(nnz_in() + nnz_out());
    }
  }

  JitFunction::~JitFunction() {
  }

  void JitFunction::codegen_body(CodeGenerator& g) const {
    // Add all input arguments as local variables
    for (casadi_int i=0; i<n_in_; ++i) {
      g.local(name_in_[i], "const casadi_real", "*");
      if (buffered_) {
        g << g.copy("*arg++", nnz_in(i), "w") << "\n"
          << name_in_[i] << " = w; w += " << nnz_in(i) << ";\n";
      } else {
        g << name_in_[i] << " = *arg++;\n";
      }
    }

    // Add all output arguments as local variables
    for (casadi_int i=0; i<n_out_; ++i) {
      g.local(name_out_[i], "casadi_real", "*");
      if (buffered_) {
        g << name_out_[i] << " = w; w += " << nnz_out(i) << ";\n";
      } else {
        g << name_out_[i] << " = *res++;\n";
      }
    }

    // Codegen function body
    g << body_;

    // Get results
    for (casadi_int i=0; i<n_out_; ++i) {
      if (buffered_) {
        g << g.copy(name_out_[i], nnz_out(i), "*res++") << "\n";
      }
    }
  }

  bool JitFunction::has_jacobian() const {
    return !jac_body_.empty();
  }

  Function JitFunction::get_jacobian(const std::string& name,
                                   const std::vector<std::string>& inames,
                                   const std::vector<std::string>& onames,
                                   const Dict& opts) const {
    // Create a JIT-function for the Jacobian
    Dict jac_opts;
    if (!hess_body_.empty()) jac_opts["jac"] = hess_body_;
    return Function::jit(name, jac_body_, inames, onames, jac_opts);
  }

} // namespace casadi
