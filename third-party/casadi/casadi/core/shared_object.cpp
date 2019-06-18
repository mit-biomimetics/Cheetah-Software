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
#include "sparse_storage_impl.hpp"
#ifdef WITH_EXTRA_CHECKS
#include "function.hpp"
#endif // WITH_EXTRA_CHECKS
#include <typeinfo>

using namespace std;
namespace casadi {

  // Instantiate templates
  template class SparseStorage<WeakRef>;

  SharedObject::SharedObject() {
    node = nullptr;
  }

  SharedObject::SharedObject(const SharedObject& ref) {
    node = ref.node;
    count_up();
  }

  SharedObject::~SharedObject() {
    count_down();
  }

  void SharedObject::own(SharedObjectInternal* node_) {
    count_down();
    node = node_;
    count_up();
  }

  void SharedObject::assign(SharedObjectInternal* node_) {
    node = node_;
  }

  SharedObject& SharedObject::operator=(const SharedObject& ref) {
    // quick return if the old and new pointers point to the same object
    if (node == ref.node) return *this;

    // decrease the counter and delete if this was the last pointer
    count_down();

    // save the new pointer
    node = ref.node;
    count_up();
    return *this;
  }

  SharedObjectInternal* SharedObject::get() const {
    return node;
  }

  bool SharedObject::is_null() const {
    return node==nullptr;
  }

  void SharedObject::count_up() {
#ifdef WITH_EXTRA_CHECKS
    casadi_assert_dev(Function::call_depth_==0);
#endif // WITH_EXTRA_CHECKS
    if (node) node->count++;
  }

  void SharedObject::count_down() {
#ifdef WITH_EXTRA_CHECKS
    casadi_assert_dev(Function::call_depth_==0);
#endif // WITH_EXTRA_CHECKS
    if (node && --node->count == 0) {
      delete node;
      node = nullptr;
    }
  }

  SharedObjectInternal* SharedObject::operator->() const {
    casadi_assert_dev(!is_null());
    return node;
  }

  std::string SharedObject::class_name() const {
    return (*this)->class_name();
  }

  void SharedObject::disp(std::ostream& stream, bool more) const {
    if (is_null()) {
      stream << "NULL";
    } else {
      (*this)->disp(stream, more);
    }
  }

  void SharedObject::print_ptr(std::ostream &stream) const {
    stream << node;
  }

  void SharedObject::swap(SharedObject& other) {
    SharedObject temp = *this;
    *this = other;
    other = temp;
  }

  casadi_int SharedObject::getCount() const {
    return (*this)->getCount();
  }

  WeakRef* SharedObject::weak() {
    return (*this)->weak();
  }

  casadi_int SharedObject::__hash__() const {
    return reinterpret_cast<casadi_int>(get());
  }

  WeakRef::WeakRef(int dummy) {
    casadi_assert_dev(dummy==0);
  }

  bool WeakRef::alive() const {
    return !is_null() && (*this)->raw_ != nullptr;
  }

  SharedObject WeakRef::shared() {
    SharedObject ret;
    if (alive()) {
      ret.own((*this)->raw_);
    }
    return ret;
  }

  const WeakRefInternal* WeakRef::operator->() const {
    return static_cast<const WeakRefInternal*>(SharedObject::operator->());
  }

  WeakRefInternal* WeakRef::operator->() {
    return static_cast<WeakRefInternal*>(SharedObject::operator->());
  }

  WeakRef::WeakRef(SharedObject shared) {
    own(shared.weak()->get());
  }

  WeakRef::WeakRef(SharedObjectInternal* raw) {
    own(new WeakRefInternal(raw));
  }

  void WeakRef::kill() {
    (*this)->raw_ = nullptr;
  }

} // namespace casadi
