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


#include "switch.hpp"

using namespace std;

namespace casadi {

  Switch::Switch(const std::string& name,
                 const std::vector<Function>& f, const Function& f_def)
    : FunctionInternal(name), f_(f), f_def_(f_def) {

    // Consitency check
    casadi_assert_dev(!f_.empty());
  }

  Switch::~Switch() {
  }

  size_t Switch::get_n_in() {
    for (auto&& i : f_) if (!i.is_null()) return 1+i.n_in();
    casadi_assert_dev(!f_def_.is_null());
    return 1+f_def_.n_in();
  }

  size_t Switch::get_n_out() {
    for (auto&& i : f_) if (!i.is_null()) return i.n_out();
    casadi_assert_dev(!f_def_.is_null());
    return f_def_.n_out();
  }

  Sparsity Switch::get_sparsity_in(casadi_int i) {
    if (i==0) {
      return Sparsity::scalar();
    } else {
      Sparsity ret;
      for (auto&& fk : f_) {
        if (!fk.is_null()) {
          const Sparsity& s = fk.sparsity_in(i-1);
          ret = ret.is_null() ? s : ret.unite(s);
        }
      }
      casadi_assert_dev(!f_def_.is_null());
      const Sparsity& s = f_def_.sparsity_in(i-1);
      ret = ret.is_null() ? s : ret.unite(s);
      return ret;
    }
  }

  Sparsity Switch::get_sparsity_out(casadi_int i) {
    Sparsity ret;
    for (auto&& fk : f_) {
      if (!fk.is_null()) {
        const Sparsity& s = fk.sparsity_out(i);
        ret = ret.is_null() ? s : ret.unite(s);
      }
    }
    casadi_assert_dev(!f_def_.is_null());
    const Sparsity& s = f_def_.sparsity_out(i);
    ret = ret.is_null() ? s : ret.unite(s);
    return ret;
  }

  void Switch::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Buffer for mismatching sparsities
    size_t sz_buf=0;

    // Keep track of sparsity projections
    project_in_ = project_out_ = false;

    // Get required work
    for (casadi_int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      if (fk.is_null()) continue;

      // Memory for evaluation
      alloc(fk);

      // Required work vectors
      size_t sz_buf_k=0;

      // Add size for input buffers
      for (casadi_int i=1; i<n_in_; ++i) {
        const Sparsity& s = fk.sparsity_in(i-1);
        if (s!=sparsity_in_[i]) {
          project_in_ = true;
          alloc_w(s.size1()); // for casadi_project
          sz_buf_k += s.nnz();
        }
      }

      // Add size for output buffers
      for (casadi_int i=0; i<n_out_; ++i) {
        const Sparsity& s = fk.sparsity_out(i);
        if (s!=sparsity_out_[i]) {
          project_out_ = true;
          alloc_w(s.size1()); // for casadi_project
          sz_buf_k += s.nnz();
        }
      }

      // Only need the largest of these work vectors
      sz_buf = max(sz_buf, sz_buf_k);
    }

    // Memory for the work vectors
    alloc_w(sz_buf, true);
  }

  int Switch::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // Get the function to be evaluated
    casadi_int k = arg[0] ? static_cast<casadi_int>(*arg[0]) : 0;
    const Function& fk = k>=0 && k<f_.size() ? f_[k] : f_def_;

    // Project arguments with different sparsity
    const double** arg1;
    if (project_in_) {
      // Project one or more argument
      arg1 = arg + n_in_;
      for (casadi_int i=0; i<n_in_-1; ++i) {
        const Sparsity& f_sp = fk.sparsity_in(i);
        const Sparsity& sp = sparsity_in_[i+1];
        arg1[i] = arg[i+1];
        if (arg1[i] && f_sp!=sp) {
          casadi_project(arg1[i], sp, w, f_sp, w + f_sp.nnz());
          arg1[i] = w; w += f_sp.nnz();
        }
      }
    } else {
      // No inputs projected
      arg1 = arg + 1;
    }

    // Temporary memory for results with different sparsity
    double** res1;
    if (project_out_) {
      // Project one or more results
      res1 = res + n_out_;
      for (casadi_int i=0; i<n_out_; ++i) {
        const Sparsity& f_sp = fk.sparsity_out(i);
        const Sparsity& sp = sparsity_out_[i];
        res1[i] = res[i];
        if (res1[i] && f_sp!=sp) {
          res1[i] = w;
          w += f_sp.nnz();
        }
      }
    } else {
      // No outputs projected
      res1 = res;
    }

    // Evaluate the corresponding function
    if (fk(arg1, res1, iw, w, 0)) return 1;

    // Project results with different sparsity
    if (project_out_) {
      for (casadi_int i=0; i<n_out_; ++i) {
        const Sparsity& f_sp = fk.sparsity_out(i);
        const Sparsity& sp = sparsity_out_[i];
        if (res[i] && f_sp!=sp) {
          casadi_project(res1[i], f_sp, res[i], sp, w);
        }
      }
    }
    return 0;
  }

  Function Switch
  ::get_forward(casadi_int nfwd, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (casadi_int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].forward(nfwd);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.forward(nfwd);

    // New Switch for derivatives
    Function sw = Function::conditional("switch_" + name, der, der_def);

    // Get expressions for the derivative switch
    vector<MX> arg = sw.mx_in();
    vector<MX> res = sw(arg);

    // Ignore seed for ind
    arg.insert(arg.begin() + n_in_ + n_out_, MX(1, nfwd));

    // Create wrapper
    return Function(name, arg, res, inames, onames, opts);
  }

  Function Switch
  ::get_reverse(casadi_int nadj, const std::string& name,
                const std::vector<std::string>& inames,
                const std::vector<std::string>& onames,
                const Dict& opts) const {
    // Derivative of each case
    vector<Function> der(f_.size());
    for (casadi_int k=0; k<f_.size(); ++k) {
      if (!f_[k].is_null()) der[k] = f_[k].reverse(nadj);
    }

    // Default case
    Function der_def;
    if (!f_def_.is_null()) der_def = f_def_.reverse(nadj);

    // New Switch for derivatives
    Function sw = Function::conditional("switch_" + name, der, der_def);

    // Get expressions for the derivative switch
    vector<MX> arg = sw.mx_in();
    vector<MX> res = sw(arg);

    // No derivatives with respect to index
    res.insert(res.begin(), MX(1, nadj));

    // Create wrapper
    return Function(name, arg, res, inames, onames, opts);
  }

  void Switch::disp_more(ostream &stream) const {
    // Print more
    if (f_.size()==1) {
      // Print as if-then-else
      stream << f_def_.name() << ", " << f_[0].name();
    } else {
      // Print generic
      stream << "[";
      for (casadi_int k=0; k<f_.size(); ++k) {
        if (k!=0) stream << ", ";
        stream << f_[k].name();
      }
      stream << "], " << f_def_.name();
    }
  }

  void Switch::codegen_declarations(CodeGenerator& g) const {
    for (casadi_int k=0; k<=f_.size(); ++k) {
      const Function& fk = k<f_.size() ? f_[k] : f_def_;
      g.add_dependency(fk);
    }
  }

  int Switch::eval_sx(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, void* mem) const {
    // Input and output buffers
    const SXElem** arg1 = arg + n_in_;
    SXElem** res1 = res + n_out_;

    // Extra memory needed for chaining if_else calls
    std::vector<SXElem> w_extra(nnz_out());
    std::vector<SXElem*> res_tempv(n_out_);
    SXElem** res_temp = get_ptr(res_tempv);

    for (casadi_int k=0; k<f_.size()+1; ++k) {

      // Local work vector
      SXElem* wl = w;

      // Local work vector
      SXElem* wll = get_ptr(w_extra);

      if (k==0) {
        // For the default case, redirect the temporary results to res
        copy_n(res, n_out_, res_temp);
      } else {
        // For the other cases, store the temporary results
        for (casadi_int i=0; i<n_out_; ++i) {
          res_temp[i] = wll;
          wll += nnz_out(i);
        }
      }

      copy_n(arg+1, n_in_-1, arg1);
      copy_n(res_temp, n_out_, res1);

      const Function& fk = k==0 ? f_def_ : f_[k-1];

      // Project arguments with different sparsity
      for (casadi_int i=0; i<n_in_-1; ++i) {
        if (arg1[i]) {
          const Sparsity& f_sp = fk.sparsity_in(i);
          const Sparsity& sp = sparsity_in_[i+1];
          if (f_sp!=sp) {
            SXElem *t = wl; wl += f_sp.nnz(); // t is non-const
            casadi_project(arg1[i], sp, t, f_sp, wl);
            arg1[i] = t;
          }
        }
      }

      // Temporary memory for results with different sparsity
      for (casadi_int i=0; i<n_out_; ++i) {
        if (res1[i]) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out_[i];
          if (f_sp!=sp) { res1[i] = wl; wl += f_sp.nnz();}
        }
      }

      // Evaluate the corresponding function
      if (fk(arg1, res1, iw, wl, 0)) return 1;

      // Project results with different sparsity
      for (casadi_int i=0; i<n_out_; ++i) {
        if (res1[i]) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out_[i];
          if (f_sp!=sp) casadi_project(res1[i], f_sp, res_temp[i], sp, wl);
        }
      }

      if (k>0) { // output the temporary results via an if_else
        SXElem cond = k-1==arg[0][0];
        for (casadi_int i=0; i<n_out_; ++i) {
          if (res[i]) {
            for (casadi_int j=0; j<nnz_out(i); ++j) {
              res[i][j] = if_else(cond, res_temp[i][j], res[i][j]);
            }
          }
        }
      }

    }
    return 0;
  }

  void Switch::codegen_body(CodeGenerator& g) const {
    // Project arguments with different sparsity
    if (project_in_) {
      // Project one or more argument
      g.local("i", "casadi_int");
      g << "const casadi_real** arg1 = arg + " << n_in_ << ";\n";
    }

    // Temporary memory for results with different sparsity
    if (project_out_) {
      // Project one or more results
      g.local("i", "casadi_int");
      g << "casadi_real** res1 = res + " << n_out_ << ";\n";
    }

    if (project_in_)
      g  << "for (i=0; i<" << n_in_-1 << "; ++i) arg1[i]=arg[i+1];\n";

    if (project_out_)
      g  << "for (i=0; i<" << n_out_ << "; ++i) res1[i]=res[i];\n";

    // Codegen condition
    bool if_else = f_.size()==1;
    g.add_auxiliary(CodeGenerator::AUX_TO_INT);
    g << (if_else ? "if" : "switch")  << " (arg[0] ? casadi_to_int(*arg[0]) : 0) {\n";

    // Loop over cases/functions
    for (casadi_int k=0; k<=f_.size(); ++k) {

      // For if,  reverse order
      casadi_int k1 = if_else ? 1-k : k;

      if (!if_else) {
        // Codegen cases
        if (k1<f_.size()) {
          g << "case " << k1 << ":\n";
        } else {
          g << "default:\n";
        }
      } else if (k1==0) {
        // Else
        g << "} else {\n";
      }

      // Get the function:
      const Function& fk = k1<f_.size() ? f_[k1] : f_def_;
      if (fk.is_null()) {
        g << "return 1;\n";
      } else {
        // Project arguments with different sparsity
        for (casadi_int i=0; i<n_in_-1; ++i) {
          const Sparsity& f_sp = fk.sparsity_in(i);
          const Sparsity& sp = sparsity_in_[i+1];
          if (f_sp!=sp) {
            if (f_sp.nnz()==0) {
              g << "arg1[" << i << "]=0;\n";
            } else {
              g.local("t", "casadi_real", "*");
              g << "t=w, w+=" << f_sp.nnz() << ";\n"
                << g.project("arg1[" + str(i) + "]", sp, "t", f_sp, "w") << "\n"
                << "arg1[" << i << "]=t;\n";
            }
          }
        }

        // Temporary memory for results with different sparsity
        for (casadi_int i=0; i<n_out_; ++i) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out_[i];
          if (f_sp!=sp) {
            if (f_sp.nnz()==0) {
              g << "res1[" << i << "]=0;\n";
            } else {
              g << "res1[" << i << "]=w, w+=" << f_sp.nnz() << ";\n";
            }
          }
        }

        // Function call
        g << "if (" << g(fk, project_in_ ? "arg1" : "arg+1",
                         project_out_ ? "res1" : "res",
                         "iw", "w") << ") return 1;\n";

        // Project results with different sparsity
        for (casadi_int i=0; i<n_out_; ++i) {
          const Sparsity& f_sp = fk.sparsity_out(i);
          const Sparsity& sp = sparsity_out_[i];
          if (f_sp!=sp) {
            g << g.project("res1[" + str(i) + "]", f_sp,
                           "res[" + str(i) + "]", sp, "w") << "\n";
          }
        }

        // Break (if switch)
        if (!if_else) g << "break;\n";
      }
    }

    // End switch/else
    g << "}\n";
  }

  Dict Switch::info() const {
    return {{"project_in", project_in_}, {"project_out", project_out_},
            {"f_def", f_def_}, {"f", f_}};
  }

} // namespace casadi
