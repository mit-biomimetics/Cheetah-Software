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


#include "linsol_internal.hpp"

using namespace std;
namespace casadi {

  LinsolInternal::LinsolInternal(const std::string& name, const Sparsity& sp)
   : ProtoFunction(name), sp_(sp) {
  }

  LinsolInternal::~LinsolInternal() {
  }

  void LinsolInternal::init(const Dict& opts) {
    // Call the base class initializer
    ProtoFunction::init(opts);

  }

  void LinsolInternal::disp(ostream &stream, bool more) const {
    stream << "Linear solver " << class_name();
    if (more) {
      stream << endl;
      disp_more(stream);
    }
  }

  int LinsolInternal::init_mem(void* mem) const {
    if (!mem) return 1;
    //auto m = static_cast<LinsolMemory*>(mem);
    return 0;
  }

  void LinsolInternal::linsol_eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w,
                                      void* mem, bool tr, casadi_int nrhs) const {
    casadi_error("eval_sx not defined for " + class_name());
  }

  int LinsolInternal::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    casadi_error("'solve' not defined for " + class_name());
  }

#if 0
  casadi_int LinsolInternal::factorize(void* mem, const double* A) const {
    // Symbolic factorization, if needed
    if (needs_sfact(mem, A)) {
      if (sfact(mem, A)) return 1;
    }

    // Numeric factorization, if needed
    if (needs_nfact(mem, A)) {
      if (nfact(mem, A)) return 1;
    }

    return 0;
  }

  bool LinsolInternal::needs_sfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolMemory*>(mem);
    return m->is_sfact;
  }

  bool LinsolInternal::needs_nfact(void* mem, const double* A) const {
    auto m = static_cast<LinsolMemory*>(mem);
    return m->is_nfact;
  }
#endif

  int LinsolInternal::nfact(void* mem, const double* A) const {
    casadi_error("'nfact' not defined for " + class_name());
  }

  casadi_int LinsolInternal::neig(void* mem, const double* A) const {
    casadi_error("'neig' not defined for " + class_name());
  }

  casadi_int LinsolInternal::rank(void* mem, const double* A) const {
    casadi_error("'rank' not defined for " + class_name());
  }

  void LinsolInternal::generate(CodeGenerator& g, const std::string& A, const std::string& x,
                                casadi_int nrhs, bool tr) const {
    g << "#error " <<  class_name() << " does not support code generation\n";
  }

  std::map<std::string, LinsolInternal::Plugin> LinsolInternal::solvers_;

  const std::string LinsolInternal::infix_ = "linsol";

} // namespace casadi
