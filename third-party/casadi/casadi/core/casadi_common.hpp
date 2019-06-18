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


#ifndef CASADI_COMMON_HPP
#define CASADI_COMMON_HPP

#include <cmath>
#include <climits>
#include <limits>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#ifdef SWIG
#define SWIG_IF_ELSE(is_swig, not_swig) is_swig
#define SWIG_OUTPUT(arg) OUTPUT
#define SWIG_INOUT(arg) INOUT
#define SWIG_CONSTREF(arg) const arg
#ifdef SWIGMATLAB
#define SWIG_IND1 true
#else // SWIGMATLAB
#define SWIG_IND1 false
#endif // SWIGMATLAB
#else // SWIG
#define SWIG_IF_ELSE(is_swig, not_swig) not_swig
#define SWIG_OUTPUT(arg) arg
#define SWIG_INOUT(arg) arg
#define SWIG_CONSTREF(arg) const arg &
#define SWIG_IND1 false
#endif // SWIG

#include "casadi_types.hpp"

namespace casadi {

  /// Forward declarations
  class SXElem;
  class MX;
  template<class T> class Matrix;
  class Function;
  class Sparsity;
  class CodeGenerator;
  class NlpBuilder;
  struct Variable;
  class DaeBuilder;
  class XmlFile;
  class Importer;

#ifndef SWIG
// Get GCC version if GCC is used
#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#ifdef __GNUC_PATCHLEVEL__
#define GCC_VERSION (__GNUC__ * 10000 +__GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif // __GNUC_PATCHLEVEL__
#endif // __GNUC_MINOR__
#endif // __GNUC__

// Disable some Visual studio warnings
#ifdef _MSC_VER

// warning C4018: '<' : signed/unsigned mismatch
#pragma warning(disable:4018)

// warning C4244: Potential loss of data converting double to int
#pragma warning(disable:4244)

// warinng C4251: Need a dll interface?
#pragma warning(disable:4251)

// warning C4715: Not all control paths return a value
#pragma warning(disable:4715)

// warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning(disable:4800)

// warning C4910: __declspec(dllexport) and extern incompatible on an explicit instantiation
#pragma warning(disable:4910)

// ?
#pragma warning(disable:4996)

#endif // _MSC_VER

  // Macro "minor" is sometimes defined, cf.
  // https://stackoverflow.com/questions/22240973/major-and-minor-macros-defined-in-sys-sysmacros-h-pulled-in-by-iterator
#undef minor

  // Type with a size corresponding to that of double (or smaller) that can be used to hold a set
  // of booleans. If the compiler supports C99 or has defined __SIZEOF_LONG_LONG__,
  // we shall use the long long datatype, which is 64 bits, otherwise long
  #if (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L || defined(__SIZEOF_LONG_LONG__))
  typedef unsigned long long bvec_t;
  #else
  typedef unsigned long bvec_t;
  #endif

  // Number of directions we can deal with at a time
  // the size of bvec_t in bits (CHAR_BIT is the number of bits per byte, usually 8)
  const int bvec_size = CHAR_BIT*sizeof(bvec_t);

  // Make sure that the integer datatype is indeed smaller or equal to the double
  //assert(sizeof(bvec_t) <= sizeof(double)); // doesn't work - very strange

  ///@{
  /** \brief  Function pointer types for the C API */
  typedef void (*signal_t)(void);
  typedef casadi_int (*getint_t)(void);
  typedef const char* (*name_t)(casadi_int i);
  typedef const casadi_int* (*sparsity_t)(casadi_int i);
  typedef void* (*alloc_mem_t)(void);
  typedef int (*init_mem_t)(void* mem);
  typedef void (*free_mem_t)(void* mem);
  typedef int (*work_t)(casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w);
  typedef int (*eval_t)(const double** arg, double** res,
                        casadi_int* iw, double* w, void* mem);
  ///@}

  /// String representation, any type
  template<typename T>
  std::string str(const T& v);

  /// String representation, CasADi type
  template<typename T>
  std::string str(const T& v, bool more);

  /// String representation of vector
  template<typename T>
  std::string str(const std::vector<T>& v, bool more=false);

  /// String representation of set
  template<typename T>
  std::string str(const std::set<T>& v, bool more=false);

  /// String representation of pair
  template<typename T1, typename T2>
  std::string str(const std::pair<T1, T2>& p, bool more=false);

  /// String representation of a map
  template<typename T1, typename T2>
  std::string str(const std::map<T1, T2> &p, bool more=false);

  /// String representation of a dictionary
  template<typename T2>
  std::string str(const std::map<std::string, T2> &p, bool more=false);


  //! \brief Create a list of strings from __VA_ARGS__, no argument
  inline std::vector<std::string> strvec() {
    return {};
  }

  //! \brief Create a list of strings from __VA_ARGS__, one argument
  template<typename T1>
  inline std::vector<std::string> strvec(const T1& t1) {
    return {str(t1)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, two arguments
  template<typename T1, typename T2>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2) {
    return {str(t1), str(t2)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, three arguments
  template<typename T1, typename T2, typename T3>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3) {
    return {str(t1), str(t2), str(t3)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, four arguments
  template<typename T1, typename T2, typename T3, typename T4>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4) {
    return {str(t1), str(t2), str(t3), str(t4)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, five arguments
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4, const T5& t5) {
    return {str(t1), str(t2), str(t3), str(t4), str(t5)};
  }

  //! \brief Create a list of strings from __VA_ARGS__, six arguments
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  inline std::vector<std::string> strvec(const T1& t1, const T2& t2, const T3& t3,
                                         const T4& t4, const T5& t5, const T6& t6) {
    return {str(t1), str(t2), str(t3), str(t4), str(t5), str(t6)};
  }

  //! \brief Create a string from a formated string
  inline std::string fmtstr(const std::string& fmt, const std::vector<std::string>& args) {
    std::string s = fmt;
    for (auto&& e : args) {
      std::string::size_type n = s.find("%s");
      if (n==std::string::npos) return "** Ill-formated string ** " + fmt;
      s.replace(n, 2, e);
    }
    return s;
  }

  // Implementations
  template<typename T>
  std::string str(const T& v) {
    std::stringstream ss;
    ss << v;
    return ss.str();
  }

  template<typename T>
  std::string str(const T& v, bool more) {
    return v.get_str(more);
  }

  template<typename T>
  std::string str(const std::vector<T>& v, bool more) {
    std::stringstream ss;
    ss << "[";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) ss << ", ";
      ss << v[i];
    }
    ss << "]";
    return ss.str();
  }

  template<typename T>
  std::string str(const std::set<T>& v, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int cnt = 0;
    for (const auto& e : v) {
      if (cnt++!=0) ss << ", ";
      ss << e;
    }
    ss << "}";
    return ss.str();
  }

  template<typename T1, typename T2>
  std::string str(const std::pair<T1, T2>& p, bool more) {
    std::stringstream ss;
    ss << "[" << p.first << "," << p.second << "]";
    return ss.str();
  }

  template<typename T1, typename T2>
  std::string str(const std::map<T1, T2>& p, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int count = 0;
    for (auto& e : p) {
      ss << e.first << ": " << e.second;
      if (++count < p.size()) ss << ", ";
    }
    ss << "}";
    return ss.str();
  }

  template<typename T2>
  std::string str(const std::map<std::string, T2>& p, bool more) {
    std::stringstream ss;
    ss << "{";
    casadi_int count = 0;
    for (auto& e : p) {
      ss << "\"" << e.first << "\": " << e.second;
      if (++count < p.size()) ss << ", ";
    }
    ss << "}";
    return ss.str();
  }
#endif // SWIG

} // namespace casadi

#include "casadi_logger.hpp"

#endif // CASADI_COMMON_HPP
