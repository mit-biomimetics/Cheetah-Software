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


#include "generic_type_internal.hpp"
#include "casadi_misc.hpp"
#include "exception.hpp"
#include <cmath>

#include "function.hpp"

using namespace std;

namespace casadi {

  /// \cond INTERNAL
  typedef GenericTypeInternal<OT_STRING, std::string> StringType;
  typedef GenericTypeInternal<OT_DOUBLE, double> DoubleType;
  typedef GenericTypeInternal<OT_INT, casadi_int> IntType;
  typedef GenericTypeInternal<OT_BOOL, bool> BoolType;
  typedef GenericTypeInternal<OT_DOUBLEVECTOR, std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<OT_DOUBLEVECTORVECTOR,
                              std::vector< std::vector<double> > > DoubleVectorVectorType;
  typedef GenericTypeInternal<OT_INTVECTOR, std::vector<casadi_int> > IntVectorType;
  typedef GenericTypeInternal<OT_INTVECTORVECTOR,
                              std::vector< std::vector<casadi_int> > > IntVectorVectorType;
  typedef GenericTypeInternal<OT_STRINGVECTOR, std::vector<std::string> > StringVectorType;
  typedef GenericTypeInternal<OT_FUNCTION, Function> FunctionType;
  typedef GenericTypeInternal<OT_FUNCTIONVECTOR, std::vector<Function> > FunctionVectorType;
  typedef GenericTypeInternal<OT_DICT, Dict> DictType;
  typedef GenericTypeInternal<OT_VOIDPTR, void*> VoidPointerType;
  /// \endcond

  bool GenericType::can_cast_to(TypeID other) const {
    switch (other) {
    case OT_BOOL:
      return is_bool() || is_int() || is_double();
    case OT_BOOLVECTOR:
      return is_int_vector() || is_double_vector();
    case OT_INT:
    case OT_DOUBLE:
      return is_int() || is_double();
    case OT_INTVECTOR:
    case OT_DOUBLEVECTOR:
      return is_double_vector() || is_int_vector();
    case OT_INTVECTORVECTOR:
    case OT_DOUBLEVECTORVECTOR:
      return is_double_vector_vector() || is_int_vector_vector();
    case OT_STRINGVECTOR:
      return is_string_vector() || is_string();
    default:
      return getType() == other;
    }
  }

  GenericType GenericType::from_type(TypeID type) {
    switch (type) {
    case OT_INTVECTOR:
      return std::vector<casadi_int>();
    case OT_INTVECTORVECTOR:
      return std::vector< std::vector<casadi_int> >();
    case OT_BOOLVECTOR:
      return std::vector<bool>();
    case OT_DOUBLEVECTOR:
      return std::vector<double>();
    case OT_DOUBLEVECTORVECTOR:
      return std::vector< std::vector<double> >();
    case OT_STRINGVECTOR:
      return std::vector<std::string>();
    default:
      casadi_error("empty_from_type. Unsupported type " + str(type));
    }
  }

  std::string GenericType::get_type_description(TypeID type) {
    switch (type) {
    case OT_BOOL:
      return "OT_BOOL";
    case OT_INT:
      return "OT_INT";
    case OT_DOUBLE:
      return "OT_DOUBLE";
    case OT_STRING:
      return "OT_STRING";
    case OT_INTVECTOR:
      return "OT_INTVECTOR";
    case OT_INTVECTORVECTOR:
      return "OT_INTVECTORVECTOR";
    case OT_BOOLVECTOR:
      return "OT_BOOLVECTOR";
    case OT_DOUBLEVECTOR:
      return "OT_DOUBLEVECTOR";
    case OT_DOUBLEVECTORVECTOR:
      return "OT_DOUBLEVECTORVECTOR";
    case OT_STRINGVECTOR:
      return "OT_STRINGVECTOR";
    case OT_DICT:
      return "OT_DICT";
    case OT_FUNCTION:
      return "OT_FUNCTION";
    case OT_FUNCTIONVECTOR:
      return "OT_FUNCTIONVECTOR";
    case OT_VOIDPTR:
      return "OT_VOIDPTR";
    default:
      return "OT_UNKNOWN";

    }
  }


  bool GenericType::is_bool() const {
    return getType()==OT_BOOL;
  }

  bool GenericType::is_int() const {
    return getType()==OT_INT;
  }

  bool GenericType::is_double() const {
    return getType()==OT_DOUBLE;
  }

  bool GenericType::is_string() const {
    return getType()==OT_STRING;
  }

  bool GenericType::is_empty_vector() const {
    return (is_int_vector() && to_int_vector().size()==0) ||
      (is_int_vector_vector() && to_int_vector_vector().size()==0) ||
      (is_double_vector_vector() && to_double_vector_vector().size()==0) ||
      (is_double_vector() && to_double_vector().size()==0) ||
      (is_string_vector() && to_string_vector().size()==0) ||
      (is_bool_vector() && to_bool_vector().size()==0);
  }

  bool GenericType::is_int_vector() const {
    return getType()==OT_INTVECTOR;
  }

  bool GenericType::is_bool_vector() const {
    return getType()==OT_BOOLVECTOR;
  }

  bool GenericType::is_int_vector_vector() const {
    return getType()==OT_INTVECTORVECTOR;
  }

  bool GenericType::is_double_vector() const {
    return getType()==OT_DOUBLEVECTOR;
  }

  bool GenericType::is_double_vector_vector() const {
    return getType()==OT_DOUBLEVECTORVECTOR;
  }

  bool GenericType::is_string_vector() const {
    return getType()==OT_STRINGVECTOR;
  }

  bool GenericType::is_function() const {
    return getType()==OT_FUNCTION;
  }

  bool GenericType::is_function_vector() const {
    return getType()==OT_FUNCTIONVECTOR;
  }

  bool GenericType::is_void_pointer() const {
    return getType()==OT_VOIDPTR;
  }

  bool GenericType::is_dict() const {
    return getType()==OT_DICT;
  }

  GenericType::GenericType() {
  }

  GenericType::GenericType(bool b) {
    own(new BoolType(b));
  }

  GenericType::GenericType(casadi_int i) {
    own(new IntType(i));
  }

  GenericType::GenericType(double d) {
    own(new DoubleType(d));
  }

  GenericType::GenericType(const vector<casadi_int>& iv) {
    own(new IntVectorType(iv));
  }

  GenericType::GenericType(const vector<int>& iv) {
    std::vector<casadi_int> temp(iv.size());
    std::copy(iv.begin(), iv.end(), temp.begin());
    own(new IntVectorType(temp));
  }

  GenericType::GenericType(const vector<vector<casadi_int> >& ivv) {
    own(new IntVectorVectorType(ivv));
  }

  GenericType::GenericType(const vector<bool>& b_vec) {
    vector<casadi_int> i_vec(b_vec.size());
    copy(b_vec.begin(), b_vec.end(), i_vec.begin());
    own(new IntVectorType(i_vec));
  }

  GenericType::GenericType(const vector<double>& dv) {
    own(new DoubleVectorType(dv));
  }

  GenericType::GenericType(const vector< vector<double> >& dv) {
    own(new DoubleVectorVectorType(dv));
  }

  GenericType::GenericType(const vector<string>& sv) {
    own(new StringVectorType(sv));
  }

  GenericType::GenericType(const string& s) {
    own(new StringType(s));
  }

  GenericType::GenericType(const char s[]) {
    own(new StringType(s));
  }

  GenericType::GenericType(const Function& f) {
    own(new FunctionType(f));
  }

  GenericType::GenericType(const std::vector<Function>& f) {
    own(new FunctionVectorType(f));
  }

  const bool& GenericType::as_bool() const {
    casadi_assert_dev(is_bool());
    return static_cast<const BoolType*>(get())->d_;
  }

  const casadi_int& GenericType::as_int() const {
    casadi_assert_dev(is_int());
    return static_cast<const IntType*>(get())->d_;
  }

  const double& GenericType::as_double() const {
    casadi_assert_dev(is_double());
    return static_cast<const DoubleType*>(get())->d_;
  }

  const std::string& GenericType::as_string() const {
    casadi_assert_dev(is_string());
    return static_cast<const StringType*>(get())->d_;
  }

  const std::vector<casadi_int>& GenericType::as_int_vector() const {
    casadi_assert_dev(is_int_vector());
    return static_cast<const IntVectorType*>(get())->d_;
  }

  const std::vector<casadi_int>& GenericType::as_bool_vector() const {
    casadi_assert_dev(is_bool_vector());
    return static_cast<const IntVectorType*>(get())->d_;
  }

  const std::vector<std::vector<casadi_int> >& GenericType::as_int_vector_vector() const {
    casadi_assert_dev(is_int_vector_vector());
    return static_cast<const IntVectorVectorType*>(get())->d_;
  }

  const std::vector<double>& GenericType::as_double_vector() const {
    casadi_assert_dev(is_double_vector());
    return static_cast<const DoubleVectorType*>(get())->d_;
  }

  const std::vector< std::vector<double> >& GenericType::as_double_vector_vector() const {
    casadi_assert_dev(is_double_vector_vector());
    return static_cast<const DoubleVectorVectorType*>(get())->d_;
  }

  const std::vector<std::string>& GenericType::as_string_vector() const {
    casadi_assert_dev(is_string_vector());
    return static_cast<const StringVectorType*>(get())->d_;
  }

  const GenericType::Dict& GenericType::as_dict() const {
    casadi_assert_dev(is_dict());
    return static_cast<const DictType*>(get())->d_;
  }

  const Function& GenericType::as_function() const {
    casadi_assert_dev(is_function());
    return static_cast<const FunctionType*>(get())->d_;
  }

  const std::vector<Function>& GenericType::as_function_vector() const {
    casadi_assert_dev(is_function_vector());
    return static_cast<const FunctionVectorType*>(get())->d_;
  }

  void* const & GenericType::as_void_pointer() const {
    casadi_assert_dev(is_void_pointer());
    return static_cast<const VoidPointerType*>(get())->d_;
  }

  bool GenericType::to_bool() const {
    if (is_bool()) {
      return as_bool();
    } else if (is_int()) {
      return static_cast<bool>(to_int());
    } else {
      casadi_assert(is_bool(), "type mismatch");
      return false;
    }
  }

  casadi_int GenericType::to_int() const {
    if (is_double()) {
      return static_cast<casadi_int>(to_double());
    } else if (is_bool()) {
      return static_cast<casadi_int>(to_bool());
    } else {
      casadi_assert(is_int(), "type mismatch");
      return as_int();
    }
  }

  double GenericType::to_double() const {
    if (is_int()) {
      return static_cast<double>(to_int());
    } else {
      casadi_assert(is_double(), "type mismatch");
      return as_double();
    }
  }

  string GenericType::to_string() const {
    casadi_assert(is_string(), "type mismatch");
    return as_string();
  }

  vector<int> GenericType::to_int_type_vector() const {
    casadi_assert(is_int_vector(), "type mismatch");
    return casadi::to_int(as_int_vector());
  }

  vector<casadi_int> GenericType::to_int_vector() const {
    casadi_assert(is_int_vector(), "type mismatch");
    return as_int_vector();
  }

  GenericType::operator std::vector<int>() const {
    std::vector<int> ret;
    std::vector<casadi_int> source = to_int_vector();
    return casadi::to_int(source);
  }

  vector<bool> GenericType::to_bool_vector() const {
    casadi_assert(is_int_vector(), "type mismatch");
    vector<casadi_int> v = to_int_vector();
    vector<bool> ret(v.size());
    for (casadi_int i=0; i<v.size(); ++i) {
      casadi_assert(v[i]==0 || v[i]==1, "Entries must be zero or one");
      ret[i] = v[i]==1;
    }
    return ret;
  }

  vector<vector<casadi_int> > GenericType::to_int_vector_vector() const {
    casadi_assert(is_int_vector_vector(), "type mismatch");
    return as_int_vector_vector();
  }

  vector<double> GenericType::to_double_vector() const {
    if (is_int_vector()) {
      auto v = as_int_vector();
      return vector<double>(v.begin(), v.end());
    } else {
      casadi_assert(is_double_vector(), "type mismatch");
      return as_double_vector();
    }
  }

  vector< vector<double> > GenericType::to_double_vector_vector() const {
    if (is_int_vector_vector()) {
      auto v = as_int_vector_vector();
      vector< vector<double> > ret(v.size());
      for (casadi_int i=0;i<v.size();++i)
        ret[i].assign(v[i].begin(), v[i].end());
      return ret;
    } else {
      casadi_assert(is_double_vector_vector(), "type mismatch");
      return as_double_vector_vector();
    }
  }

  vector<string> GenericType::to_string_vector() const {
    if (is_string()) {
      std::string s = as_string();
      return vector<string>(1, s);
    } else {
      casadi_assert(is_string_vector(), "type mismatch");
      return as_string_vector();
    }
  }

  Dict GenericType::to_dict() const {
    casadi_assert(is_dict(), "type mismatch");
    return as_dict();
  }

  Function GenericType::to_function() const {
    casadi_assert(is_function(), "type mismatch");
    return as_function();
  }

  std::vector<Function> GenericType::to_function_vector() const {
    casadi_assert(is_function_vector(), "type mismatch");
    return as_function_vector();
  }

  bool GenericType::operator==(const GenericType& op2) const {
    return !(*this != op2);
  }

  bool GenericType::operator!=(const GenericType& op2) const {
    if (is_string() && op2.is_string()) {
      return to_string().compare(op2.to_string()) != 0;
    }

    if (is_int() && op2.is_int()) {
      return to_int() != op2.to_int();
    }

    if (is_double() && op2.is_double()) {
      return to_double() != op2.to_double();
    }

    if (is_double_vector() && op2.is_double_vector()) {
      const vector<double> &v1 = to_double_vector();
      const vector<double> &v2 = op2.to_double_vector();
      if (v1.size() != v2.size()) return true;
      for (casadi_int i=0; i<v1.size(); ++i)
        if (v1[i] != v2[i]) return true;
      return false;
    }

    if (is_int_vector() && op2.is_int_vector()) {
      const vector<casadi_int> &v1 = to_int_vector();
      const vector<casadi_int> &v2 = op2.to_int_vector();
      if (v1.size() != v2.size()) return true;
      for (casadi_int i=0; i<v1.size(); ++i)
        if (v1[i] != v2[i]) return true;
      return false;
    }

    if (is_int_vector_vector() && op2.is_int_vector_vector()) {
      const vector< vector<casadi_int> > &v1 = to_int_vector_vector();
      const vector< vector<casadi_int> > &v2 = op2.to_int_vector_vector();
      if (v1.size() != v2.size()) return true;
      for (casadi_int i=0; i<v1.size(); ++i) {
        if (v1[i].size() != v2[i].size()) return true;
        for (casadi_int j=0; j<v1[i].size(); ++j) {
          if (v1[i][j] != v2[i][j]) return true;
        }
      }
      return false;
    }

    if (is_double_vector_vector() && op2.is_double_vector_vector()) {
      const vector< vector<double> > &v1 = to_double_vector_vector();
      const vector< vector<double> > &v2 = op2.to_double_vector_vector();
      if (v1.size() != v2.size()) return true;
      for (casadi_int i=0; i<v1.size(); ++i) {
        if (v1[i].size() != v2[i].size()) return true;
        for (casadi_int j=0; j<v1[i].size(); ++j) {
          if (v1[i][j] != v2[i][j]) return true;
        }
      }
      return false;
    }

    // Different types
    return true;
  }

  GenericType::GenericType(const Dict& dict) {
    own(new DictType(dict));
  }

  void* GenericType::to_void_pointer() const {
    casadi_assert(getType()==OT_VOIDPTR, "type mismatch");
    return as_void_pointer();
  }

  GenericType::GenericType(void* ptr) {
    own(new VoidPointerType(ptr));
  }

  TypeID GenericType::getType() const {
    if (is_null()) {
      return OT_NULL;
    } else {
      return static_cast<const GenericTypeBase*>(get())->getType();
    }
  }

} // namespace casadi
