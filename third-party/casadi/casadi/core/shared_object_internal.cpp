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

#include "shared_object_internal.hpp"

using namespace std;
namespace casadi {

  SharedObjectInternal::SharedObjectInternal(const SharedObjectInternal& node) {
    count = 0; // reference counter is _not_ copied
    weak_ref_ = nullptr; // nor will they have the same weak references
  }

  SharedObjectInternal& SharedObjectInternal::operator=(const SharedObjectInternal& node) {
    // do _not_ copy the reference counter
    return *this;
  }

  SharedObjectInternal::SharedObjectInternal() {
    count = 0;
    weak_ref_ = nullptr;
  }

  SharedObjectInternal::~SharedObjectInternal() {
    #ifdef WITH_REFCOUNT_WARNINGS
    if (count!=0) {
      // Note that casadi_assert_warning cannot be used in destructors
      std::cerr << "Reference counting failure." <<
                   "Possible cause: Circular dependency in user code." << std::endl;
    }
    #endif // WITH_REFCOUNT_WARNINGS
    if (weak_ref_!=nullptr) {
      weak_ref_->kill();
      delete weak_ref_;
    }
  }

  casadi_int SharedObjectInternal::getCount() const {
    return count;
  }

  WeakRef* SharedObjectInternal::weak() {
    if (weak_ref_==nullptr) {
      weak_ref_ = new WeakRef(this);
    }
    return weak_ref_;
  }

  WeakRefInternal::WeakRefInternal(SharedObjectInternal* raw) : raw_(raw) {
  }

  WeakRefInternal::~WeakRefInternal() {
  }

  void WeakRefInternal::disp(std::ostream& stream, bool more) const {
    if (raw_==nullptr) {
      stream << "NULL";
    } else {
      raw_->disp(stream, more);
    }
  }


} // namespace casadi
