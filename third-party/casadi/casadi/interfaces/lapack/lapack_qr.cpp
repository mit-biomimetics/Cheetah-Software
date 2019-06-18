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


#include "lapack_qr.hpp"
#include "../../core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LAPACKQR_EXPORT
  casadi_register_linsol_lapackqr(LinsolInternal::Plugin* plugin) {
    plugin->creator = LapackQr::creator;
    plugin->name = "lapackqr";
    plugin->doc = LapackQr::meta_doc.c_str();;
    plugin->version = CASADI_VERSION;
    plugin->options = &LapackQr::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LAPACKQR_EXPORT casadi_load_linsol_lapackqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lapackqr);
  }

  LapackQr::LapackQr(const std::string& name, const Sparsity& sp) :
    LinsolInternal(name, sp) {
  }

  LapackQr::~LapackQr() {
    clear_mem();
  }

  Options LapackQr::options_
  = {{&FunctionInternal::options_},
     {{"max_nrhs",
       {OT_INT,
        "Maximum number of right-hand-sides that get processed in a single pass [default:10]."}}
     }
  };

  void LapackQr::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);

    max_nrhs_ = 10;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="max_nrhs") {
        max_nrhs_ = op.second;
      }
    }
  }

  int LapackQr::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<LapackQrMemory*>(mem);
    m->mat.resize(ncol() * ncol());
    m->tau.resize(ncol());
    m->work.resize(max(max_nrhs_, ncol())*10);
    return 0;
  }

  int LapackQr::nfact(void* mem, const double* A) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    // Dimensions
    //casadi_int nrow = this->nrow();
    int ncol = this->ncol();

    // Get the elements of the matrix, dense format
    casadi_densify(A, sp_, get_ptr(m->mat), false);

    // Factorize the matrix
    int info = -100;
    int lwork = m->work.size();
    dgeqrf_(&ncol, &ncol, get_ptr(m->mat), &ncol, get_ptr(m->tau),
            get_ptr(m->work), &lwork, &info);
    if (info) {
      if (verbose_) casadi_warning("dgeqrf_ failed: Info: " + str(info));
      return 1;
    }
    return 0;
  }

  int LapackQr::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    // Solve up to max_nrhs rhs at a time
    casadi_int offset = 0;
    while (nrhs>0) {
      if (solve_batch(m, A, x+offset, min(max_nrhs_, nrhs), tr)) return 1;
      nrhs-= max_nrhs_;
      offset+= max_nrhs_*nrow();
    }
    return 0;
  }

  int LapackQr::solve_batch(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<LapackQrMemory*>(mem);

    // Dimensions
    //casadi_int nrow = this->nrow();
    int ncol = this->ncol();

    // Properties of R
    char uploR = 'U';
    char diagR = 'N';
    char sideR = 'L';
    double alphaR = 1.;
    char transR = tr ? 'T' : 'N';

    // Properties of Q
    char transQ = tr ? 'N' : 'T';
    char sideQ = 'L';
    int k = m->tau.size(); // minimum of ncol and nrow
    int lwork = m->work.size();

    int n_rhs = nrhs;

    if (tr) {

      // Solve for transpose(R)
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &n_rhs, &alphaR,
             get_ptr(m->mat), &ncol, x, &ncol);

      // Multiply by Q
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &n_rhs, &k, get_ptr(m->mat), &ncol, get_ptr(m->tau), x,
              &ncol, get_ptr(m->work), &lwork, &info);
      casadi_assert(info == 0,
        "LapackQr::solve: dormqr_ A failed to solve the linear system. Info: " + str(info) + ".");
    } else {

      // Multiply by transpose(Q)
      int info = 100;
      dormqr_(&sideQ, &transQ, &ncol, &n_rhs, &k, get_ptr(m->mat), &ncol, get_ptr(m->tau), x,
              &ncol, get_ptr(m->work), &lwork, &info);
      casadi_assert(info == 0,
        "LapackQr::solve: dormqr_ B failed to solve the linear system. Info: " + str(info) + ".");

      // Solve for R
      dtrsm_(&sideR, &uploR, &transR, &diagR, &ncol, &n_rhs, &alphaR,
             get_ptr(m->mat), &ncol, x, &ncol);
    }
    return 0;
  }

} // namespace casadi
