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



#include "code_generator.hpp"
#include "function_internal.hpp"
#include <iomanip>
#include <casadi_runtime_str.h>

using namespace std;
namespace casadi {

  CodeGenerator::CodeGenerator(const string& name, const Dict& opts) {
    // Default options
    this->verbose = true;
    this->mex = false;
    this->cpp = false;
    this->main = false;
    this->casadi_real = "double";
    this->casadi_int_type = CASADI_INT_TYPE_STR;
    this->codegen_scalars = false;
    this->with_header = false;
    this->with_mem = false;
    this->with_export = true;
    this->with_import = false;
    this->include_math = true;
    avoid_stack_ = false;
    indent_ = 2;

    // Read options
    for (auto&& e : opts) {
      if (e.first=="verbose") {
        this->verbose = e.second;
      } else if (e.first=="mex") {
        this->mex = e.second;
      } else if (e.first=="cpp") {
        this->cpp = e.second;
      } else if (e.first=="main") {
        this->main = e.second;
      } else if (e.first=="casadi_real") {
        this->casadi_real = e.second.to_string();
      }  else if (e.first=="casadi_int") {
        this->casadi_int_type = e.second.to_string();
      } else if (e.first=="codegen_scalars") {
        this->codegen_scalars = e.second;
      } else if (e.first=="with_header") {
        this->with_header = e.second;
      } else if (e.first=="with_mem") {
        this->with_mem = e.second;
      } else if (e.first=="with_export") {
        this->with_export = e.second;
      } else if (e.first=="with_import") {
        this->with_import = e.second;
      } else if (e.first=="include_math") {
        this->include_math = e.second;
      } else if (e.first=="indent") {
        indent_ = e.second;
        casadi_assert_dev(indent_>=0);
      } else if (e.first=="avoid_stack") {
        avoid_stack_ = e.second;
      } else {
        casadi_error("Unrecongnized option: " + str(e.first));
      }
    }

    // Start at new line with no indentation
    newline_ = true;
    current_indent_ = 0;

    // Divide name into base and suffix (if any)
    string::size_type dotpos = name.rfind('.');
    if (dotpos==string::npos) {
      this->name = name;
      this->suffix = this->cpp ? ".cpp" : ".c";
    } else {
      this->name = name.substr(0, dotpos);
      this->suffix = name.substr(dotpos);
    }

    // Symbol prefix
    if (this->with_export) dll_export = "CASADI_SYMBOL_EXPORT ";
    if (this->with_import) dll_import = "CASADI_SYMBOL_IMPORT ";

    // Make sure that the base name is sane
    casadi_assert_dev(Function::check_name(this->name));

    // Includes needed
    if (this->include_math) add_include("math.h");
    if (this->main) add_include("stdio.h");

    // Mex and main need string.h
    if (this->mex || this->main) {
      add_include("string.h");
    }

    // Memory struct entry point
    if (this->with_mem) {
      add_include("casadi/mem.h", false);
      this->header << "#include <casadi/mem.h>\n";
    }

    // Mex
    if (this->mex) {
      add_include("mex.h", false, "MATLAB_MEX_FILE");
    }
  }

  string CodeGenerator::add_dependency(const Function& f) {
    // Quick return if it already exists
    for (auto&& e : added_functions_) if (e.f==f) return e.codegen_name;

    // Give it a name
    string fname = shorthand("f" + str(added_functions_.size()));

    // Add to list of functions
    added_functions_.push_back({f, fname});

    // Generate declarations
    f->codegen_declarations(*this);

    // Print to file
    f->codegen(*this, fname);

    // Codegen reference count functions, if needed
    if (f->has_refcount_) {
      // Increase reference counter
      *this << "void " << fname << "_incref(void) {\n";
      f->codegen_incref(*this);
      *this << "}\n\n";

      // Decrease reference counter
      *this << "void " << fname << "_decref(void) {\n";
      f->codegen_decref(*this);
      *this << "}\n\n";
    }

    // Flush to body
    flush(this->body);

    return fname;
  }

    void CodeGenerator::add(const Function& f, bool with_jac_sparsity) {
    // Add if not already added
    string codegen_name = add_dependency(f);

    // Define function
    *this << declare(f->signature(f.name())) << "{\n"
          << "return " << codegen_name <<  "(arg, res, iw, w, mem);\n"
          << "}\n\n";

    // Generate meta information
    f->codegen_meta(*this);

    // Generate Jacobian sparsity information
    if (with_jac_sparsity) {
      // Generate/get Jacobian sparsity
      Sparsity jac = f->get_jacobian_sparsity();
      // Code generate the sparsity pattern
      add_io_sparsities("jac_" + f.name(), f->sparsity_in_, {jac});

      // Flush buffers
      flush(this->body);
    }

    // Add to list of exposed symbols
    this->exposed_fname.push_back(f.name());
  }

  string CodeGenerator::dump() const {
    stringstream s;
    dump(s);
    return s.str();
  }

  void CodeGenerator::file_open(std::ofstream& f, const string& name) const {
    // Open a file for writing
    f.open(name);

    // Print header
    f << "/* This file was automatically generated by CasADi.\n"
      << "   The CasADi copyright holders make no ownership claim of its contents. */\n";

    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus\n"
        << "extern \"C\" {\n"
        << "#endif\n\n";
    }
  }

  void CodeGenerator::file_close(std::ofstream& f) const {
    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus\n"
        << "} /* extern \"C\" */\n"
        << "#endif\n";
    }

    // Close file(s)
    f.close();
  }

  void CodeGenerator::generate_casadi_real(std::ostream &s) const {
    s << "#ifndef casadi_real\n"
      << "#define casadi_real " << this->casadi_real << endl
      << "#endif\n\n";
  }

  void  CodeGenerator::generate_export_symbol(std::ostream &s) const {
      s << "/* Symbol visibility in DLLs */\n"
      << "#ifndef CASADI_SYMBOL_EXPORT\n"
      << "  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)\n"
      << "    #if defined(STATIC_LINKED)\n"
      << "      #define CASADI_SYMBOL_EXPORT\n"
      << "    #else\n"
      << "      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)\n"
      << "    #endif\n"
      << "  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)\n"
      << "    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility (\"default\")))\n"
      << "  #else"  << endl
      << "    #define CASADI_SYMBOL_EXPORT\n"
      << "  #endif\n"
      << "#endif\n\n";
  }

  void  CodeGenerator::generate_import_symbol(std::ostream &s) const {
      s << "/* Symbol visibility in DLLs */\n"
      << "#ifndef CASADI_SYMBOL_IMPORT\n"
      << "  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)\n"
      << "    #if defined(STATIC_LINKED)\n"
      << "      #define CASADI_SYMBOL_IMPORT\n"
      << "    #else\n"
      << "      #define CASADI_SYMBOL_IMPORT __declspec(dllimport)\n"
      << "    #endif\n"
      << "  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)\n"
      << "    #define CASADI_SYMBOL_IMPORT __attribute__ ((visibility (\"default\")))\n"
      << "  #else"  << endl
      << "    #define CASADI_SYMBOL_IMPORT\n"
      << "  #endif\n"
      << "#endif\n\n";
  }

  void CodeGenerator::generate_casadi_int(std::ostream &s) const {
    s << "#ifndef casadi_int\n"
      << "#define casadi_int " << this->casadi_int_type << endl
      << "#endif\n\n";
  }

  string CodeGenerator::generate(const string& prefix) const {
    // Throw an error if the prefix contains the filename, since since syntax
    // has changed
    casadi_assert(prefix.find(this->name + this->suffix)==string::npos,
       "The signature of CodeGenerator::generate has changed. "
       "Instead of providing the filename, only provide the prefix.");

    // Create c file
    ofstream s;
    string fullname = prefix + this->name + this->suffix;
    file_open(s, fullname);

    // Dump code to file
    dump(s);

    // Mex entry point
    if (this->mex) generate_mex(s);

    // Main entry point
    if (this->main) generate_main(s);

    // Finalize file
    file_close(s);

    // Generate header
    if (this->with_header) {
      // Create a header file
      file_open(s, prefix + this->name + ".h");

      // Define the casadi_real type (typically double)
      generate_casadi_real(s);

      // Define the casadi_int type
      generate_casadi_int(s);

      // Generate export symbol macros
      if (this->with_import) generate_import_symbol(s);

      // Add declarations
      s << this->header.str();

      // Finalize file
      file_close(s);
    }
    return fullname;
  }

  void CodeGenerator::generate_mex(std::ostream &s) const {
    // Begin conditional compilation
    s << "#ifdef MATLAB_MEX_FILE\n";

    // Function prototype
    if (this->cpp) s << "extern \"C\"\n"; // C linkage
    s << "void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {"
      << endl;

    // Create a buffer
    size_t buf_len = 0;
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      buf_len = std::max(buf_len, exposed_fname[i].size());
    }
    s << "  char buf[" << (buf_len+1) << "];\n";

    // Read string argument
    s << "  int buf_ok = --argc >= 0 && !mxGetString(*argv++, buf, sizeof(buf));\n";

    // Create switch
    s << "  if (!buf_ok) {\n"
      << "    /* name error */\n";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(buf, \"" << exposed_fname[i] << "\")==0) {\n"
        << "    mex_" << exposed_fname[i] << "(resc, resv, argc, argv);\n"
        << "    return;\n";
    }
    s << "  }\n";

    // Error
    s << "  mexErrMsgTxt(\"First input should be a command string. Possible values:";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\");\n";

    // End conditional compilation and function
    s << "}\n"
         << "#endif\n";
  }

  void CodeGenerator::generate_main(std::ostream &s) const {
    s << "int main(int argc, char* argv[]) {\n";

    // Create switch
    s << "  if (argc<2) {\n"
      << "    /* name error */\n";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(argv[1], \"" << exposed_fname[i] << "\")==0) {\n"
        << "    return main_" << exposed_fname[i] << "(argc-2, argv+2);\n";
    }
    s << "  }\n";

    // Error
    s << "  fprintf(stderr, \"First input should be a command string. Possible values:";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\\n\");\n";

    // End main
    s << "  return 1;\n"
      << "}\n";
  }

  void CodeGenerator::dump(std::ostream& s) const {
    // Consistency check
    casadi_assert_dev(current_indent_ == 0);

    // Prefix internal symbols to avoid symbol collisions
    s << "/* How to prefix internal symbols */\n"
      << "#ifdef CODEGEN_PREFIX\n"
      << "  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)\n"
      << "  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID\n"
      << "  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)\n"
      << "#else\n"
      << "  #define CASADI_PREFIX(ID) " << this->name << "_ ## ID\n"
      << "#endif\n\n";

    s << this->includes.str();
    s << endl;

    // Real type (usually double)
    generate_casadi_real(s);

    // Integer type (usually long long)
    generate_casadi_int(s);

    // Macros
    if (!added_shorthands_.empty()) {
      s << "/* Add prefix to internal symbols */\n";
      for (auto&& i : added_shorthands_) {
        s << "#define " << "casadi_" << i <<  " CASADI_PREFIX(" << i <<  ")\n";
      }
      s << endl;
    }

    if (this->with_export) generate_export_symbol(s);

    // Print integer constants
    if (!integer_constants_.empty()) {
      for (casadi_int i=0; i<integer_constants_.size(); ++i) {
        print_vector(s, "casadi_s" + str(i), integer_constants_[i]);
      }
      s << endl;
    }

    // Print double constants
    if (!double_constants_.empty()) {
      for (casadi_int i=0; i<double_constants_.size(); ++i) {
        print_vector(s, "casadi_c" + str(i), double_constants_[i]);
      }
      s << endl;
    }

    // External function declarations
    if (!added_externals_.empty()) {
      s << "/* External functions */\n";
      for (auto&& i : added_externals_) {
        s << i << endl;
      }
      s << endl << endl;
    }

    // Codegen auxiliary functions
    s << this->auxiliaries.str();

    // Codegen body
    s << this->body.str();

    // End with new line
    s << endl;
  }

  string CodeGenerator::work(casadi_int n, casadi_int sz) const {
    if (n<0 || sz==0) {
      return "0";
    } else if (sz==1 && !this->codegen_scalars) {
      return "(&w" + str(n) + ")";
    } else {
      return "w" + str(n);
    }
  }

  string CodeGenerator::workel(casadi_int n) const {
    if (n<0) return "0";
    stringstream s;
    if (this->codegen_scalars) s << "*";
    s << "w" << n;
    return s.str();
  }

  string CodeGenerator::array(const string& type, const string& name, casadi_int len,
                                   const string& def) {
    stringstream s;
    s << type << " ";
    if (len==0) {
      s << "*" << name << " = 0";
    } else {
      s << name << "[" << len << "]";
      if (!def.empty()) s << " = " << def;
    }
    s << ";\n";
    return s.str();
  }

  void CodeGenerator::print_vector(std::ostream &s, const string& name,
      const vector<casadi_int>& v) {
    s << array("static const casadi_int", name, v.size(), initializer(v));
  }

  void CodeGenerator::print_vector(std::ostream &s, const string& name,
                                  const vector<double>& v) {
    s << array("static const casadi_real", name, v.size(), initializer(v));
  }

  std::string CodeGenerator::print_op(casadi_int op, const std::string& a0) {
    switch (op) {
      case OP_SQ:
        add_auxiliary(AUX_SQ);
        return "casadi_sq("+a0+")";
      case OP_SIGN:
        add_auxiliary(AUX_SIGN);
        return "casadi_sign("+a0+")";
      default:
        return casadi_math<double>::print(op, a0);
    }
  }
  std::string CodeGenerator::print_op(casadi_int op, const std::string& a0, const std::string& a1) {
    switch (op) {
      case OP_FMIN:
        add_auxiliary(AUX_FMIN);
        return "casadi_fmin("+a0+","+a1+")";
      case OP_FMAX:
        add_auxiliary(AUX_FMAX);
        return "casadi_fmax("+a0+","+a1+")";
      default:
        return casadi_math<double>::print(op, a0, a1);
    }
  }

  void CodeGenerator::add_include(const string& new_include, bool relative_path,
                                 const string& use_ifdef) {
    // Register the new element
    bool added = added_includes_.insert(new_include).second;

    // Quick return if it already exists
    if (!added) return;

    // Ifdef opening
    if (!use_ifdef.empty()) this->includes << "#ifdef " << use_ifdef << endl;

    // Print to the header section
    if (relative_path) {
      this->includes << "#include \"" << new_include << "\"\n";
    } else {
      this->includes << "#include <" << new_include << ">\n";
    }

    // Ifdef closing
    if (!use_ifdef.empty()) this->includes << "#endif\n";
  }

  string CodeGenerator::
  operator()(const Function& f, const string& arg,
             const string& res, const string& iw,
             const string& w, const string& mem) const {
    return f->codegen_name(*this) + "(" + arg + ", " + res + ", "
      + iw + ", " + w + ", " + mem + ")";
  }

  void CodeGenerator::add_external(const string& new_external) {
    added_externals_.insert(new_external);
  }

  string CodeGenerator::shorthand(const string& name) const {
    casadi_assert(added_shorthands_.count(name), "No such macro: " + name);
    return "casadi_" + name;
  }

  string CodeGenerator::shorthand(const string& name, bool allow_adding) {
    bool added = added_shorthands_.insert(name).second;
    if (!allow_adding) {
      casadi_assert(added, "Duplicate macro: " + name);
    }
    return "casadi_" + name;
  }

  casadi_int CodeGenerator::add_sparsity(const Sparsity& sp) {
    return get_constant(sp, true);
  }

  string CodeGenerator::sparsity(const Sparsity& sp) {
    return shorthand("s" + str(add_sparsity(sp)));
  }

  casadi_int CodeGenerator::get_sparsity(const Sparsity& sp) const {
    return const_cast<CodeGenerator&>(*this).get_constant(sp, false);
  }

  size_t CodeGenerator::hash(const vector<double>& v) {
    // Calculate a hash value for the vector
    std::size_t seed=0;
    if (!v.empty()) {
      casadi_assert_dev(sizeof(double) % sizeof(size_t)==0);
      const casadi_int int_len = v.size()*(sizeof(double)/sizeof(size_t));
      const size_t* int_v = reinterpret_cast<const size_t*>(&v.front());
      for (size_t i=0; i<int_len; ++i) {
        hash_combine(seed, int_v[i]);
      }
    }
    return seed;
  }

  size_t CodeGenerator::hash(const vector<casadi_int>& v) {
    size_t seed=0;
    hash_combine(seed, v);
    return seed;
  }

  casadi_int CodeGenerator::get_constant(const vector<double>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    auto eq = added_double_constants_.equal_range(h);
    for (auto i=eq.first; i!=eq.second; ++i) {
      if (equal(v, double_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      casadi_int ind = double_constants_.size();
      double_constants_.push_back(v);
      added_double_constants_.insert(make_pair(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  casadi_int CodeGenerator::get_constant(const vector<casadi_int>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator> eq =
      added_integer_constants_.equal_range(h);
    for (multimap<size_t, size_t>::iterator i=eq.first; i!=eq.second; ++i) {
      if (equal(v, integer_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      casadi_int ind = integer_constants_.size();
      integer_constants_.push_back(v);
      added_integer_constants_.insert(pair<size_t, size_t>(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  string CodeGenerator::constant(const vector<casadi_int>& v) {
    return shorthand("s" + str(get_constant(v, true)));
  }

  string CodeGenerator::constant(const vector<double>& v) {
    return shorthand("c" + str(get_constant(v, true)));
  }

  void CodeGenerator::add_auxiliary(Auxiliary f, const vector<string>& inst) {
    // Look for existing instantiations
    auto f_match = added_auxiliaries_.equal_range(f);
    // Look for duplicates
    for (auto it=f_match.first; it!=f_match.second; ++it) {
      if (it->second==inst) return;
    }
    added_auxiliaries_.insert(make_pair(f, inst));

    // Add the appropriate function
    switch (f) {
    case AUX_COPY:
      this->auxiliaries << sanitize_source(casadi_copy_str, inst);
      break;
    case AUX_SWAP:
      this->auxiliaries << sanitize_source(casadi_swap_str, inst);
      break;
    case AUX_SCAL:
      this->auxiliaries << sanitize_source(casadi_scal_str, inst);
      break;
    case AUX_AXPY:
      this->auxiliaries << sanitize_source(casadi_axpy_str, inst);
      break;
    case AUX_DOT:
      this->auxiliaries << sanitize_source(casadi_dot_str, inst);
      break;
    case AUX_BILIN:
      this->auxiliaries << sanitize_source(casadi_bilin_str, inst);
      break;
    case AUX_RANK1:
      this->auxiliaries << sanitize_source(casadi_rank1_str, inst);
      break;
    case AUX_IAMAX:
      this->auxiliaries << sanitize_source(casadi_iamax_str, inst);
      break;
    case AUX_INTERPN:
      add_auxiliary(AUX_INTERPN_WEIGHTS);
      add_auxiliary(AUX_INTERPN_INTERPOLATE);
      add_auxiliary(AUX_FLIP, {});
      add_auxiliary(AUX_FILL);
      add_auxiliary(AUX_FILL, {"casadi_int"});
      this->auxiliaries << sanitize_source(casadi_interpn_str, inst);
      break;
    case AUX_INTERPN_GRAD:
      add_auxiliary(AUX_INTERPN);
      this->auxiliaries << sanitize_source(casadi_interpn_grad_str, inst);
      break;
    case AUX_DE_BOOR:
      this->auxiliaries << sanitize_source(casadi_de_boor_str, inst);
      break;
    case AUX_ND_BOOR_EVAL:
      add_auxiliary(AUX_DE_BOOR);
      add_auxiliary(AUX_FILL);
      add_auxiliary(AUX_FILL, {"casadi_int"});
      add_auxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_nd_boor_eval_str, inst);
      break;
    case AUX_FLIP:
      this->auxiliaries << sanitize_source(casadi_flip_str, inst);
      break;
    case AUX_LOW:
      this->auxiliaries << sanitize_source(casadi_low_str, inst);
      break;
    case AUX_INTERPN_WEIGHTS:
      add_auxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_interpn_weights_str, inst);
      break;
    case AUX_INTERPN_INTERPOLATE:
      this->auxiliaries << sanitize_source(casadi_interpn_interpolate_str, inst);
      break;
    case AUX_NORM_1:
      this->auxiliaries << sanitize_source(casadi_norm_1_str, inst);
      break;
    case AUX_NORM_2:
      this->auxiliaries << sanitize_source(casadi_norm_2_str, inst);
      break;
    case AUX_NORM_INF:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_norm_inf_str, inst);
      break;
    case AUX_FILL:
      this->auxiliaries << sanitize_source(casadi_fill_str, inst);
      break;
    case AUX_MV:
      this->auxiliaries << sanitize_source(casadi_mv_str, inst);
      break;
    case AUX_MV_DENSE:
      this->auxiliaries << sanitize_source(casadi_mv_dense_str, inst);
      break;
    case AUX_MTIMES:
      this->auxiliaries << sanitize_source(casadi_mtimes_str, inst);
      break;
    case AUX_PROJECT:
      this->auxiliaries << sanitize_source(casadi_project_str, inst);
      break;
    case AUX_DENSIFY:
      add_auxiliary(AUX_FILL);
      add_auxiliary(AUX_CAST);
      {
        vector<string> inst2 = inst;
        if (inst.size()==1) inst2.push_back(inst[0]);
        this->auxiliaries << sanitize_source(casadi_densify_str, inst2);
      }
      break;
    case AUX_TRANS:
      this->auxiliaries << sanitize_source(casadi_trans_str, inst);
      break;
    case AUX_TO_MEX:
      add_auxiliary(AUX_TO_DOUBLE);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_to_mex_str, inst)
                        << "#endif\n\n";
      break;
    case AUX_FROM_MEX:
      add_auxiliary(AUX_FILL);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_from_mex_str, inst)
                        << "#endif\n\n";
      break;
    case AUX_FINITE_DIFF:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_finite_diff_str, inst);
      break;
    case AUX_QR:
      add_auxiliary(AUX_IF_ELSE);
      add_auxiliary(AUX_SCAL);
      add_auxiliary(AUX_DOT);
      add_auxiliary(AUX_FILL);
      this->auxiliaries << sanitize_source(casadi_qr_str, inst);
      break;
    case AUX_LDL:
      this->auxiliaries << sanitize_source(casadi_ldl_str, inst);
      break;
    case AUX_NEWTON:
      add_auxiliary(AUX_COPY);
      add_auxiliary(AUX_AXPY);
      add_auxiliary(AUX_NORM_INF);
      add_auxiliary(AUX_QR);
      this->auxiliaries << sanitize_source(casadi_newton_str, inst);
      break;
    case AUX_TO_DOUBLE:
      this->auxiliaries << "#define casadi_to_double(x) "
                        << "(" << (this->cpp ? "static_cast<double>(x)" : "(double) x") << ")\n\n";
      break;
    case AUX_TO_INT:
      this->auxiliaries << "#define casadi_to_int(x) "
                        << "(" << (this->cpp ? "static_cast<casadi_int>(x)" : "(casadi_int) x")
                        << ")\n\n";
      break;
    case AUX_CAST:
      this->auxiliaries << "#define CASADI_CAST(x,y) "
                        << "(" << (this->cpp ? "static_cast<x>(y)" : "(x) y") << ")\n\n";
      break;
    case AUX_SQ:
      shorthand("sq");
      this->auxiliaries << "casadi_real casadi_sq(casadi_real x) { return x*x;}\n\n";
      break;
    case AUX_SIGN:
      shorthand("sign");
      this->auxiliaries << "casadi_real casadi_sign(casadi_real x) "
                        << "{ return x<0 ? -1 : x>0 ? 1 : x;}\n\n";
      break;
    case AUX_IF_ELSE:
      shorthand("if_else");
      this->auxiliaries << "casadi_real casadi_if_else"
                        << "(casadi_real c, casadi_real x, casadi_real y) "
                        << "{ return c!=0 ? x : y;}\n\n";
      break;
    case AUX_PRINTF:
      if (this->mex) {
        this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                          << "  #define CASADI_PRINTF mexPrintf\n"
                          << "#else\n"
                          << "  #define CASADI_PRINTF printf\n"
                          << "#endif\n\n";
      } else {
        this->auxiliaries << "#define CASADI_PRINTF printf\n\n";
      }
      break;
    case AUX_FMIN:
      shorthand("fmin");
      this->auxiliaries << "casadi_real casadi_fmin(casadi_real x, casadi_real y) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x<y ? x : y;\n"
                        << "#else\n"
                        << "  return fmin(x, y);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_FMAX:
      shorthand("fmax");
      this->auxiliaries << "casadi_real casadi_fmax(casadi_real x, casadi_real y) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x>y ? x : y;\n"
                        << "#else\n"
                        << "  return fmax(x, y);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    }
  }

  string CodeGenerator::to_mex(const Sparsity& sp, const string& arg) {
    add_auxiliary(AUX_TO_MEX);
    stringstream s;
    s << "casadi_to_mex(" << sparsity(sp) << ", " << arg << ");";
    return s.str();
  }

  string CodeGenerator::from_mex(string& arg,
                                      const string& res, std::size_t res_off,
                                      const Sparsity& sp_res, const string& w) {
    // Handle offset with recursion
    if (res_off!=0) return from_mex(arg, res+"+"+str(res_off), 0, sp_res, w);

    add_auxiliary(AUX_FROM_MEX);
    stringstream s;
    s << "casadi_from_mex(" << arg
      << ", " << res << ", " << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  string CodeGenerator::constant(casadi_int v) {
    return constant(static_cast<double>(v));
  }
  string CodeGenerator::constant(double v) {
    stringstream s;
    if (isnan(v)) {
      s << "NAN";
    } else if (isinf(v)) {
      if (v<0) s << "-";
      s << "INFINITY";
    } else {
      casadi_int v_int = static_cast<casadi_int>(v);
      if (static_cast<double>(v_int)==v) {
        // Print integer
        s << v_int << ".";
      } else {
        // Print real
        std::ios_base::fmtflags fmtfl = s.flags(); // get current format flags
        s << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << v;
        s.flags(fmtfl); // reset current format flags
      }
    }
    return s.str();
  }

  string CodeGenerator::initializer(const vector<double>& v) {
    stringstream s;
    s << "{";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << constant(v[i]);
    }
    s << "}";
    return s.str();
  }

  string CodeGenerator::initializer(const vector<casadi_int>& v) {
    stringstream s;
    s << "{";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << v[i];
    }
    s << "}";
    return s.str();
  }

  string CodeGenerator::copy(const string& arg,
                                  std::size_t n, const string& res) {
    stringstream s;
    // Perform operation
    add_auxiliary(AUX_COPY);
    s << "casadi_copy(" << arg << ", " << n << ", " << res << ");";
    return s.str();
  }

  string CodeGenerator::fill(const string& res,
                                  std::size_t n, const string& v) {
    stringstream s;
    // Perform operation
    add_auxiliary(AUX_FILL);
    s << "casadi_fill(" << res << ", " << n << ", " << v << ");";
    return s.str();
  }

  string CodeGenerator::dot(casadi_int n, const string& x,
                                 const string& y) {
    add_auxiliary(AUX_DOT);
    stringstream s;
    s << "casadi_dot(" << n << ", " << x << ", " << y << ")";
    return s.str();
  }

  string CodeGenerator::bilin(const string& A, const Sparsity& sp_A,
                                   const string& x, const string& y) {
    add_auxiliary(AUX_BILIN);
    stringstream s;
    s << "casadi_bilin(" << A << ", " << sparsity(sp_A) << ", " << x << ", " << y << ")";
    return s.str();
  }

  string CodeGenerator::rank1(const string& A, const Sparsity& sp_A,
                                   const string& alpha, const string& x,
                                   const string& y) {
    add_auxiliary(AUX_RANK1);
    stringstream s;
    s << "casadi_rank1(" << A << ", " << sparsity(sp_A) << ", "
      << alpha << ", " << x << ", " << y << ");";
    return s.str();
  }

  string CodeGenerator::interpn(const std::string& res, casadi_int ndim, const string& grid,
                                   const string& offset,
                                   const string& values, const string& x,
                                   const string& lookup_mode, casadi_int m,
                                   const string& iw, const string& w) {
    add_auxiliary(AUX_INTERPN);
    stringstream s;
    s << "casadi_interpn(" << res << ", " << ndim << ", " << grid << ", "  << offset << ", "
      << values << ", " << x << ", " << lookup_mode << ", " << m << ", " << iw << ", " << w << ");";
    return s.str();
  }

  string CodeGenerator::interpn_grad(const string& grad,
                                   casadi_int ndim, const string& grid, const string& offset,
                                   const string& values, const string& x,
                                   const string& lookup_mode, casadi_int m,
                                   const string& iw, const string& w) {
    add_auxiliary(AUX_INTERPN_GRAD);
    stringstream s;
    s << "casadi_interpn_grad(" << grad << ", " << ndim << ", " << grid << ", " << offset << ", "
      << values << ", " << x << ", " << lookup_mode << "," << m << ", " << iw << ", " << w << ");";
    return s.str();
  }

  string CodeGenerator::trans(const string& x, const Sparsity& sp_x,
                                   const string& y, const Sparsity& sp_y,
                                   const string& iw) {
    add_auxiliary(CodeGenerator::AUX_TRANS);
    return "casadi_trans(" + x + "," + sparsity(sp_x) + ", "
            + y + ", " + sparsity(sp_y) + ", " + iw + ")";
  }

  string CodeGenerator::declare(string s) {
    // Add c linkage
    string cpp_prefix = this->cpp ? "extern \"C\" " : "";

    // To header file
    if (this->with_header) {
      this->header << cpp_prefix << this->dll_import << s << ";\n";
    }

    // Return name with declarations
    return cpp_prefix + this->dll_export + s;
  }

  string
  CodeGenerator::project(const string& arg, const Sparsity& sp_arg,
                         const string& res, const Sparsity& sp_res,
                         const string& w) {
    // If sparsity match, simple copy
    if (sp_arg==sp_res) return copy(arg, sp_arg.nnz(), res);

    // Create call
    add_auxiliary(AUX_PROJECT);
    stringstream s;
    s << "casadi_project(" << arg << ", " << sparsity(sp_arg) << ", " << res << ", "
      << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  string CodeGenerator::printf(const string& str, const vector<string>& arg) {
    add_include("stdio.h");
    add_auxiliary(AUX_PRINTF);
    stringstream s;
    s << "CASADI_PRINTF(\"" << str << "\"";
    for (casadi_int i=0; i<arg.size(); ++i) s << ", " << arg[i];
    s << ");";
    return s.str();
  }

  string CodeGenerator::printf(const string& str, const string& arg1) {
    vector<string> arg;
    arg.push_back(arg1);
    return printf(str, arg);
  }

  string CodeGenerator::printf(const string& str, const string& arg1,
                                    const string& arg2) {
    vector<string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    return printf(str, arg);
  }

  string CodeGenerator::printf(const string& str, const string& arg1,
                                    const string& arg2, const string& arg3) {
    vector<string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    arg.push_back(arg3);
    return printf(str, arg);
  }

  string CodeGenerator::axpy(casadi_int n, const string& a,
                                  const string& x, const string& y) {
    add_auxiliary(AUX_AXPY);
    return "casadi_axpy(" + str(n) + ", " + a + ", " + x + ", " + y + ");";
  }

  string CodeGenerator::scal(casadi_int n, const string& alpha, const string& x) {
    add_auxiliary(AUX_SCAL);
    return "casadi_scal(" + str(n) + ", " + alpha + ", " + x + ");";
  }

  string CodeGenerator::mv(const string& x, const Sparsity& sp_x,
                                const string& y, const string& z, bool tr) {
    add_auxiliary(AUX_MV);
    return "casadi_mv(" + x + ", " + sparsity(sp_x) + ", " + y + ", "
           + z + ", " +  (tr ? "1" : "0") + ");";
  }

  string CodeGenerator::mv(const string& x, casadi_int nrow_x, casadi_int ncol_x,
                                const string& y, const string& z, bool tr) {
    add_auxiliary(AUX_MV_DENSE);
    return "casadi_mv_dense(" + x + ", " + str(nrow_x) + ", " + str(ncol_x) + ", "
           + y + ", " + z + ", " +  (tr ? "1" : "0") + ");";
  }

  string CodeGenerator::mtimes(const string& x, const Sparsity& sp_x,
                                    const string& y, const Sparsity& sp_y,
                                    const string& z, const Sparsity& sp_z,
                                    const string& w, bool tr) {
    add_auxiliary(AUX_MTIMES);
    return "casadi_mtimes(" + x + ", " + sparsity(sp_x) + ", " + y + ", " + sparsity(sp_y) + ", "
      + z + ", " + sparsity(sp_z) + ", " + w + ", " +  (tr ? "1" : "0") + ");";
  }

  void CodeGenerator::print_formatted(const string& s) {
    // Quick return if empty
    if (s.empty()) return;

    // If new line, add indentation
    if (newline_) {
      casadi_int shift = s.front()=='}' ? -1 : 0;
      casadi_assert_dev(current_indent_+shift>=0);
      this->buffer << string(indent_*(current_indent_+shift), ' ');
      newline_ = false;
    }

    // Print to body
    this->buffer << s;

    // Brackets change indentation for next row
    // NOTE(@jaeandersson): Should ignore strings, comments
    for (char c : s) {
      if (c=='{') {
        indent();
      } else if (c=='}') {
        unindent();
      }
    }
  }

  CodeGenerator& CodeGenerator::operator<<(const string& s) {
    // Loop over newline characters
    size_t off=0;
    while (true) {
      size_t pos = s.find('\n', off);
      if (pos==string::npos) {
        // No more newline characters
        print_formatted(s.substr(off));
        break;
      } else {
        // Ends with newline
        print_formatted(s.substr(off, pos-off));
        this->buffer << '\n';
        newline_ = true;
        off = pos+1;
      }
    }

    return *this;
  }

  void CodeGenerator::flush(std::ostream &s) {
    s << this->buffer.str();
    this->buffer.str(string());
  }

  void CodeGenerator::local(const string& name, const string& type,
                            const string& ref) {
    // Check if the variable already exists
    auto it = local_variables_.find(name);
    if (it==local_variables_.end()) {
      // Add it
      local_variables_[name] = make_pair(type, ref);
    } else {
      // Consistency check
      casadi_assert(it->second.first==type, "Type mismatch for " + name);
      casadi_assert(it->second.second==ref, "Type mismatch for " + name);
    }
  }

  std::string CodeGenerator::sx_work(casadi_int i) {
    if (avoid_stack_) {
      return "w[" + str(i) + "]";
    } else {
      std::string name = "a"+str(i);

      // Make sure work vector element has been declared
      local(name, "casadi_real");

      return name;
    }
  }

  void CodeGenerator::init_local(const string& name, const string& def) {
    bool inserted = local_default_.insert(make_pair(name, def)).second;
    casadi_assert(inserted, name + " already defined");
  }

  string CodeGenerator::
  sanitize_source(const string& src,
                  const vector<string>& inst, bool add_shorthand) {
    // Create suffix if templates type are not all "casadi_real"
    string suffix;
    for (const string& s : inst) {
      if (s!="casadi_real") {
        for (const string& s : inst) suffix += "_" + s;
        break;
      }
    }

    // Construct map of name replacements
    vector<std::pair<string, string> > rep;
    for (casadi_int i=0; i<inst.size(); ++i) {
      rep.push_back(make_pair("T" + str(i+1), inst[i]));
    }

    // Return object
    stringstream ret;
    // Process C++ source
    string line;
    istringstream stream(src);
    while (std::getline(stream, line)) {
      size_t n1, n2;

      // C++ template declarations are ignored
      if (line.find("template")==0) continue;

      // Macro definitions are ignored
      if (line.find("#define")==0) continue;
      if (line.find("#undef")==0) continue;

      // Inline declaration
      if (line == "inline") continue;

      // If line starts with "// SYMBOL", add shorthand
      if (line.find("// SYMBOL") != string::npos) {
        n1 = line.find("\"");
        n2 = line.find("\"", n1+1);
        string sym = line.substr(n1+1, n2-n1-1);
        if (add_shorthand) shorthand(sym + suffix);
        if (!suffix.empty()) {
          rep.push_back(make_pair(sym, sym + suffix));
        }
        continue;
      }

      // If line starts with "// C-REPLACE", add to list of replacements
      if (line.find("// C-REPLACE") != string::npos) {
        // Get C++ string
        n1 = line.find("\"");
        n2 = line.find("\"", n1+1);
        string key = line.substr(n1+1, n2-n1-1);
        // Get C string
        n1 = line.find("\"", n2+1);
        n2 = line.find("\"", n1+1);
        string sub = line.substr(n1+1, n2-n1-1);
        // Add to replacements
        rep.push_back(make_pair(key, sub));
        continue;
      }

      // Ignore other C++ style comment
      if ((n1 = line.find("//")) != string::npos) line.erase(n1);

      // Remove trailing spaces
      if ((n1 = line.find_last_not_of(' ')) != string::npos) {
        line.erase(n1 + 1);
      } else {
        continue;
      }

      // Perform string replacements
      for (auto&& it = rep.rbegin(); it!=rep.rend(); ++it) {
        string::size_type n = 0;
        while ((n = line.find(it->first, n)) != string::npos) {
          line.replace(n, it->first.size(), it->second);
          n += it->second.size();
        }
      }

      // Append to return
      ret << line << "\n";
    }

    // Trailing newline
    ret << "\n";
    return ret.str();
  }

  void CodeGenerator::comment(const string& s) {
    if (verbose) {
      *this << "/* " << s << " */\n";
    }
  }

  void CodeGenerator::
  add_io_sparsities(const string& name,
                    const vector<Sparsity>& sp_in,
                    const vector<Sparsity>& sp_out) {
    // Insert element, quick return if it already exists
    if (!sparsity_meta.insert(name).second) return;

    // Input sparsities
    *this << declare("const casadi_int* " + name + "_sparsity_in(casadi_int i)") << " {\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<sp_in.size(); ++i) {
      *this << "case " << i << ": return " << sparsity(sp_in[i]) << ";\n";
    }
    *this << "default: return 0;\n}\n"
      << "}\n\n";

    // Output sparsities
    *this << declare("const casadi_int* " + name + "_sparsity_out(casadi_int i)") << " {\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<sp_out.size(); ++i) {
      *this << "case " << i << ": return " << sparsity(sp_out[i]) << ";\n";
    }
    *this << "default: return 0;\n}\n"
      << "}\n\n";
  }

  string CodeGenerator::
  qr(const string& sp, const string& A, const string& w,
     const string& sp_v, const string& v, const string& sp_r,
     const string& r, const string& beta, const string& prinv, const string& pc) {
    add_auxiliary(CodeGenerator::AUX_QR);
    return "casadi_qr(" + sp + ", " + A + ", " + w + ", "
           + sp_v + ", " + v + ", " + sp_r + ", " + r + ", "
           + beta + ", " + prinv + ", " + pc + ");";
  }

  string CodeGenerator::
  qr_solve(const string& x, casadi_int nrhs, bool tr,
           const string& sp_v, const string& v,
           const string& sp_r, const string& r,
           const string& beta, const string& prinv,
           const string& pc, const string& w) {
    add_auxiliary(CodeGenerator::AUX_QR);
    return "casadi_qr_solve(" + x + ", " + str(nrhs) + ", " + (tr ? "1" : "0") + ", "
           + sp_v + ", " + v + ", " + sp_r + ", " + r + ", "
           + beta + ", " + prinv + ", " + pc + ", " + w + ");";
  }

  std::string CodeGenerator::
  ldl(const std::string& sp_a, const std::string& a,
      const std::string& sp_lt, const std::string& lt, const std::string& d,
      const std::string& p, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_LDL);
    return "casadi_ldl(" + sp_a + ", " + a + ", " + sp_lt + ", " + lt + ", "
           + d + ", " + p + ", " + w + ");";
  }

  std::string CodeGenerator::
  ldl_solve(const std::string& x, casadi_int nrhs,
    const std::string& sp_lt, const std::string& lt, const std::string& d,
    const std::string& p, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_LDL);
    return "casadi_ldl_solve(" + x + ", " + str(nrhs) + ", " + sp_lt + ", "
           + lt + ", " + d + ", " + p + ", " + w + ");";
  }

} // namespace casadi
