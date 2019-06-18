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


#include "slicot_expm.hpp"
#include "slicot_layer.hpp"
#include "slicot_la.hpp"

#include "../../core/casadi_misc.hpp"
#include "../../core/mx_function.hpp"
#include "../../core/sx_function.hpp"

#include <cassert>
#include <ctime>
#include <numeric>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_EXPM_SLICOT_EXPORT
  casadi_register_expm_slicot(Expm::Plugin* plugin) {
    plugin->creator = SlicotExpm::creator;
    plugin->name = "slicot";
    plugin->doc = SlicotExpm::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &SlicotExpm::options_;
    return 0;
  }

  extern "C"
  void CASADI_EXPM_SLICOT_EXPORT casadi_load_expm_slicot() {
    Expm::registerPlugin(casadi_register_expm_slicot);
  }

  SlicotExpm::SlicotExpm(const std::string& name, const Sparsity& A) : Expm(name, A) {

  }

  SlicotExpm::~SlicotExpm() {
    clear_mem();
  }

  bool SlicotExpm::has_loaded_ = false;

  void SlicotExpm::init(const Dict& opts) {

    if (!has_loaded_) {
      has_loaded_ = true;
      casadi_warning("Loaded plugin with GPL license.");
    }

    Expm::init(opts);

    n_ = A_.size1();

    alloc_w(n_*n_, true);   // A
    alloc_w(n_*n_, true);   // H
    alloc_w(2*n_*n_, false); // dwork
    alloc_iw(2*n_, false);   // iwork
  }


  void SlicotExpm::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<SlicotExpmMemory*>(mem);

    // Set work in base classes
    Expm::set_work(mem, arg, res, iw, w);

    m->A = w; w += n_*n_;
    m->H = w; w += n_*n_;
    m->dwork = w;
    m->iwork = reinterpret_cast<int*>(iw);
  }


  /** \brief Initalize memory block */
  int SlicotExpm::init_mem(void* mem) const {
    return Expm::init_mem(mem);
  }

  int SlicotExpm::eval(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<SlicotExpmMemory*>(mem);

    setup(mem, arg, res, iw, w);

    double tol = 1e-8;
    int ret = slicot_mb05nd(n_, arg[1][0], arg[0], n_, m->A, n_, m->H, n_,
      tol, m->iwork, m->dwork, 2*n_*n_);
    casadi_assert(ret==0, "Slicot mb05nd failed with status " + str(ret) + ".");
    if (res[0]) std::copy(m->A, m->A+n_*n_, res[0]);
    return 0;
  }




} // namespace casadi
