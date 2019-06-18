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


#include "lapack_lu.hpp"
#include "../../core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LAPACKLU_EXPORT
  casadi_register_linsol_lapacklu(LinsolInternal::Plugin* plugin) {
    plugin->creator = LapackLu::creator;
    plugin->name = "lapacklu";
    plugin->doc = LapackLu::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LapackLu::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKLU_EXPORT casadi_load_linsol_lapacklu() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lapacklu);
  }

  LapackLu::LapackLu(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {

    // Default options
    equilibriate_ = true;
    allow_equilibration_failure_ = false;
  }

  LapackLu::~LapackLu() {
    clear_mem();
  }

  Options LapackLu::options_
  = {{&FunctionInternal::options_},
     {{"equilibration",
       {OT_BOOL,
        "Equilibrate the matrix"}},
      {"allow_equilibration_failure",
       {OT_BOOL,
        "Non-fatal error when equilibration fails"}}
     }
  };

  void LapackLu::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="equilibration") {
        equilibriate_ = op.second;
      } else if (op.first=="allow_equilibration_failure") {
        allow_equilibration_failure_ = op.second;
      }
    }
  }

  int LapackLu::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 0;
    auto m = static_cast<LapackLuMemory*>(mem);

    // Allocate matrix
    m->mat.resize(nrow() * ncol());
    m->ipiv.resize(ncol());

    // Equilibration
    if (equilibriate_) {
      m->r.resize(nrow());
      m->c.resize(ncol());
    }
    m->equed = 'N'; // No equilibration
    return 0;
  }

  int LapackLu::nfact(void* mem, const double* A) const {
    auto m = static_cast<LapackLuMemory*>(mem);

    // Dimensions
    int nrow = this->nrow();
    int ncol = this->ncol();

    // Get the elements of the matrix, dense format
    casadi_densify(A, sp_, get_ptr(m->mat), false);

    if (equilibriate_) {
      // Calculate the col and row scaling factors
      double colcnd, rowcnd; // ratio of the smallest to the largest col/row scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&ncol, &nrow, get_ptr(m->mat), &ncol, get_ptr(m->r),
              get_ptr(m->c), &colcnd, &rowcnd, &amax, &info);
      if (info < 0) return 1;
      if (info > 0) {
        stringstream ss;
        ss << "LapackLu::prepare: ";
        if (info<=ncol)  ss << (info-1) << "-th row (zero-based) is exactly zero";
        else             ss << (info-1-ncol) << "-th col (zero-based) is exactly zero";
        uout() << "Warning: " << ss.str() << endl;
        if (allow_equilibration_failure_)  uout() << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }

      // Equilibrate the matrix if scaling was successful
      if (info!=0)
        dlaqge_(&ncol, &nrow, get_ptr(m->mat), &ncol, get_ptr(m->r), get_ptr(m->c),
                &colcnd, &rowcnd, &amax, &m->equed);
      else
        m->equed = 'N';
    }

    // Factorize the matrix
    int info = -100;
    dgetrf_(&ncol, &ncol, get_ptr(m->mat), &ncol, get_ptr(m->ipiv), &info);
    if (info) {
      if (verbose_) casadi_warning("dgetrf_ failed: Info: " + str(info));
      return 1;
    }
    return 0;
  }

  int LapackLu::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LapackLuMemory*>(mem);

    // Dimensions
    int nrow = this->nrow();
    int ncol = this->ncol();

    int n_rhs = nrhs;

    // Scale the right hand side
    if (tr) {
      if (m->equed=='C' || m->equed=='B')
        for (casadi_int rhs=0; rhs<nrhs; ++rhs)
          for (casadi_int i=0; i<nrow; ++i)
            x[i+rhs*nrow] *= m->c[i];
    } else {
      if (m->equed=='R' || m->equed=='B')
        for (casadi_int rhs=0; rhs<nrhs; ++rhs)
          for (casadi_int i=0; i<ncol; ++i)
            x[i+rhs*nrow] *= m->r[i];
    }

    // Solve the system of equations
    int info = 100;
    char trans = tr ? 'T' : 'N';
    dgetrs_(&trans, &ncol, &n_rhs, get_ptr(m->mat), &ncol, get_ptr(m->ipiv), x, &ncol, &info);
    if (info) return 1;

    // Scale the solution
    if (tr) {
      if (m->equed=='R' || m->equed=='B')
        for (casadi_int rhs=0; rhs<nrhs; ++rhs)
          for (casadi_int i=0; i<ncol; ++i)
            x[i+rhs*nrow] *= m->r[i];
    } else {
      if (m->equed=='C' || m->equed=='B')
        for (casadi_int rhs=0; rhs<nrhs; ++rhs)
          for (casadi_int i=0; i<nrow; ++i)
            x[i+rhs*nrow] *= m->c[i];
    }
    return 0;
  }

} // namespace casadi
