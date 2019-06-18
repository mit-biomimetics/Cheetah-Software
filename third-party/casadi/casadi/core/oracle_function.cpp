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


#include "oracle_function.hpp"
#include "external.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

namespace casadi {

  OracleFunction::OracleFunction(const std::string& name, const Function& oracle)
  : FunctionInternal(name), oracle_(oracle) {
  }

  OracleFunction::~OracleFunction() {
  }

  Options OracleFunction::options_
  = {{&FunctionInternal::options_},
     {{"monitor",
       {OT_STRINGVECTOR,
        "Set of user problem functions to be monitored"}},
      {"common_options",
       {OT_DICT,
        "Options for auto-generated functions"}},
      {"specific_options",
       {OT_DICT,
        "Options for specific auto-generated functions,"
        " overwriting the defaults from common_options. Nested dictionary."}}
    }
  };

  void OracleFunction::init(const Dict& opts) {

    FunctionInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="common_options") {
        common_options_ = op.second;
      } else if (op.first=="specific_options") {
        specific_options_ = op.second;
        for (auto&& i : specific_options_) {
          casadi_assert(i.second.is_dict(),
            "specific_option must be a nested dictionary."
            " Type mismatch for entry '" + i.first+ "': "
            " got type " + i.second.get_description() + ".");
        }
      }
    }
  }

  void OracleFunction::finalize(const Dict& opts) {
    // Default options
    vector<string> monitor;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="monitor") {
        monitor = op.second;
      }
    }

    // Set corresponding monitors
    for (const string& fname : monitor) {
      auto it = all_functions_.find(fname);
      if (it==all_functions_.end()) {
        casadi_warning("Ignoring monitor '" + fname + "'."
                       " Available functions: " + join(get_function()) + ".");
      } else {
        if (it->second.monitored) casadi_warning("Duplicate monitor " + fname);
        it->second.monitored = true;
      }
    }

    // Check specific options
    for (auto&& i : specific_options_) {
      if (all_functions_.find(i.first)==all_functions_.end())
        casadi_warning("Ignoring specific_options entry '" + i.first+"'."
                       " Available functions: " + join(get_function()) + ".");
    }

    // Recursive call
    FunctionInternal::finalize(opts);
  }

  Function OracleFunction::create_function(const std::string& fname,
                                   const std::vector<std::string>& s_in,
                                   const std::vector<std::string>& s_out,
                                   const Function::AuxOut& aux) {
    // Print progress
    if (verbose_) {
      casadi_message(name_ + "::create_function " + fname + ":" + str(s_in) + "->" + str(s_out));
    }

    // Retrieve specific set of options if available
    Dict specific_options;
    auto it = specific_options_.find(fname);
    if (it!=specific_options_.end()) specific_options = it->second;

    // Combine specific and common options
    Dict opt = combine(specific_options, common_options_);

    // Generate the function
    Function ret = oracle_.factory(fname, s_in, s_out, aux, opt);

    // Make sure that it's sound
    if (ret.has_free()) {
      casadi_error("Cannot create '" + fname + "' since " + str(ret.get_free()) + " are free.");
    }

    // Save and return
    set_function(ret, fname, true);
    return ret;
  }

  void OracleFunction::
  set_function(const Function& fcn, const std::string& fname, bool jit) {
    casadi_assert(!has_function(fname), "Duplicate function " + fname);
    RegFun& r = all_functions_[fname];
    r.f = fcn;
    r.jit = jit;
    alloc(fcn);
  }

  casadi_int OracleFunction::
  calc_function(OracleMemory* m, const std::string& fcn,
                const double* const* arg) const {
    // Is the function monitored?
    bool monitored = this->monitored(fcn);

    // Print progress
    if (monitored) casadi_message("Calling \"" + fcn + "\"");

    // Respond to a possible Crl+C signals
    InterruptHandler::check();

    // Get function
    const Function& f = get_function(fcn);

    // Get statistics structure
    FStats& fstats = m->fstats.at(fcn);

    // Number of inputs and outputs
    casadi_int n_in = f.n_in(), n_out = f.n_out();

    // Prepare stats, start timer
    fstats.tic();

    // Input buffers
    if (arg) {
      fill_n(m->arg, n_in, nullptr);
      for (casadi_int i=0; i<n_in; ++i) m->arg[i] = *arg++;
    }

    // Print inputs nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " input nonzeros:\n";
      for (casadi_int i=0; i<n_in; ++i) {
        s << " " << i << " (" << f.name_in(i) << "): ";
        if (m->arg[i]) {
          // Print nonzeros
          s << "[";
          for (casadi_int k=0; k<f.nnz_in(i); ++k) {
            if (k!=0) s << ", ";
            DM::print_scalar(s, m->arg[i][k]);
          }
          s << "]\n";
        } else {
          // All zero input
          s << "0\n";
        }
      }
      casadi_message(s.str());
    }

    // Evaluate memory-less
    try {
      f(m->arg, m->res, m->iw, m->w);
    } catch(exception& ex) {
      // Fatal error
      casadi_warning(name_ + ":" + fcn + " failed:" + std::string(ex.what()));
      return 1;
    }

    // Print output nonzeros
    if (monitored) {
      std::stringstream s;
      s << fcn << " output nonzeros:\n";
      for (casadi_int i=0; i<n_out; ++i) {
        s << " " << i << " (" << f.name_out(i) << "): ";
        if (m->res[i]) {
          // Print nonzeros
          s << "[";
          for (casadi_int k=0; k<f.nnz_out(i); ++k) {
            if (k!=0) s << ", ";
            DM::print_scalar(s, m->res[i][k]);
          }
          s << "]\n";
        } else {
          // Ignored output
          s << " N/A\n";
        }
      }
      casadi_message(s.str());
    }

    // Make sure not NaN or Inf
    for (casadi_int i=0; i<n_out; ++i) {
      if (!m->res[i]) continue;
      if (!all_of(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return isfinite(v);})) {
        std::stringstream ss;

        auto it = find_if(m->res[i], m->res[i]+f.nnz_out(i), [](double v) { return !isfinite(v);});
        casadi_int k = distance(m->res[i], it);
        bool is_nan = isnan(m->res[i][k]);
        ss << name_ << ":" << fcn << " failed: " << (is_nan? "NaN" : "Inf") <<
        " detected for output " << f.name_out(i) << ", at " << f.sparsity_out(i).repr_el(k) << ".";

        if (regularity_check_) {
          casadi_error(ss.str());
        } else {
          casadi_warning(ss.str());
        }
        return -1;
      }
    }

    // Update stats
    fstats.toc();

    // Success
    return 0;
  }

  std::string OracleFunction::
  generate_dependencies(const std::string& fname, const Dict& opts) const {
    CodeGenerator gen(fname, opts);
    gen.add(oracle_);
    for (auto&& e : all_functions_) {
      if (e.second.jit) gen.add(e.second.f);
    }
    return gen.generate();
  }

  void OracleFunction::jit_dependencies(const std::string& fname) {
    if (verbose_)
      if (verbose_) casadi_message("compiling to "+ fname+"'.");
    // JIT dependent functions
    compiler_ = Importer(generate_dependencies(fname, Dict()),
                         compilerplugin_, jit_options_);

    // Replace the Oracle functions with generated functions
    for (auto&& e : all_functions_) {
      if (verbose_)
        if (verbose_) casadi_message("loading '" + e.second.f.name() + "' from '" + fname + "'.");
      if (e.second.jit) e.second.f = external(e.second.f.name(), compiler_);
    }
  }

  void OracleFunction::expand() {
    oracle_ = oracle_.expand();
  }

  void OracleFunction::print_fstats(const OracleMemory* m) const {
    // Length of the name being printed
    size_t name_len=0;
    for (auto &&s : m->fstats) {
      name_len = max(s.first.size(), name_len);
    }

    // Print name with a given length. Format: "%NNs "
    char namefmt[10];
    sprint(namefmt, sizeof(namefmt), "%%%ds ", static_cast<casadi_int>(name_len));

    // Print header
    print(namefmt, "");
    print("%12s %12s %9s\n", "t_proc [s]", "t_wall [s]", "n_eval");

    // Print keys
    for (auto &&s : m->fstats) {
      const FStats& fs = m->fstats.at(s.first);
      if (fs.n_call!=0) {
        print(namefmt, s.first.c_str());
        print("%12.3g %12.3g %9d\n", fs.t_proc, fs.t_wall, fs.n_call);
      }
    }
  }

  Dict OracleFunction::get_stats(void *mem) const {
    auto m = static_cast<OracleMemory*>(mem);

    // Add timing statistics
    Dict stats;
    for (auto&& s : m->fstats) {
      stats["n_call_" +s.first] = s.second.n_call;
      stats["t_wall_" +s.first] = s.second.t_wall;
      stats["t_proc_" +s.first] = s.second.t_proc;
    }
    return stats;
  }

  int OracleFunction::init_mem(void* mem) const {
    if (!mem) return 1;
    auto m = static_cast<OracleMemory*>(mem);

    // Create statistics
    for (auto&& e : all_functions_) {
      m->fstats[e.first] = FStats();
    }
    return 0;
  }

  void OracleFunction::set_temp(void* mem, const double** arg, double** res,
                            casadi_int* iw, double* w) const {
    auto m = static_cast<OracleMemory*>(mem);
    m->arg = arg;
    m->res = res;
    m->iw = iw;
    m->w = w;
  }

  std::vector<std::string> OracleFunction::get_function() const {
    std::vector<std::string> ret;
    ret.reserve(all_functions_.size());
    for (auto&& e : all_functions_) {
      ret.push_back(e.first);
    }
    return ret;
  }

  const Function& OracleFunction::get_function(const std::string &name) const {
    auto it = all_functions_.find(name);
    casadi_assert(it!=all_functions_.end(),
      "No function \"" + name + "\" in " + name_ + ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.f;
  }

  bool OracleFunction::monitored(const std::string &name) const {
    auto it = all_functions_.find(name);
    casadi_assert(it!=all_functions_.end(),
      "No function \"" + name + "\" in " + name_+ ". " +
      "Available functions: " + join(get_function()) + ".");
    return it->second.monitored;
  }

  bool OracleFunction::has_function(const std::string& fname) const {
    return all_functions_.find(fname) != all_functions_.end();
  }


} // namespace casadi
