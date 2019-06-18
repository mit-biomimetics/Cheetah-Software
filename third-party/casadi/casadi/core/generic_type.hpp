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


#ifndef CASADI_GENERIC_TYPE_HPP
#define CASADI_GENERIC_TYPE_HPP

#include "shared_object.hpp"
#include "printable.hpp"
#include "casadi_common.hpp"
#include <string>
#include <vector>

namespace casadi {

#if !(defined(SWIG) && !defined(SWIGXML))

  /** \brief  Types of options */
  enum TypeID {
    OT_NULL,
    OT_BOOL,
    OT_INT,
    OT_DOUBLE,
    OT_STRING,
    OT_INTVECTOR,
    OT_INTVECTORVECTOR,
    OT_BOOLVECTOR,
    OT_DOUBLEVECTOR,
    OT_DOUBLEVECTORVECTOR,
    OT_STRINGVECTOR,
    OT_DICT,
    OT_FUNCTION,
    OT_FUNCTIONVECTOR,
    OT_VOIDPTR,
    OT_UNKNOWN};
#endif // SWIG

  /** \brief Generic data type, can hold different types such as bool, casadi_int, string etc.
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT GenericType
    : public SWIG_IF_ELSE(PrintableCommon, Printable<GenericType>)
#if !(defined(SWIG) && !defined(SWIGXML))
    , public SharedObject
#endif // SWIG
  {
  public:
    /// C++ equivalent of Python's dict or MATLAB's struct
    typedef std::map<std::string, GenericType> Dict;

#if !(defined(SWIG) && !defined(SWIGXML))

    /// Default constructor
    GenericType();

    /// Constructors (implicit type conversion)
    GenericType(bool b);
    GenericType(casadi_int i);
    GenericType(int i) : GenericType(static_cast<casadi_int>(i)) {}
    GenericType(double d);
    GenericType(const std::string& s);
    GenericType(const std::vector<bool>& iv);
    GenericType(const std::vector<casadi_int>& iv);
    GenericType(const std::vector<int>& iv);
    GenericType(const std::vector< std::vector<casadi_int> >& ivv);
    GenericType(const std::vector<double>& dv);
    GenericType(const std::vector< std::vector<double> >& dv);
    GenericType(const std::vector<std::string>& sv);
    GenericType(const char s[]);
    GenericType(const Function& f);
    GenericType(const std::vector<Function>& f);
    GenericType(const Dict& dict);
    GenericType(void* ptr);

    /// Public class name
    static std::string type_name() {return "GenericType";}

#ifndef SWIG
    /// Get a description of a type
    static std::string get_type_description(TypeID type);
#endif

    /// Get a description of the object's type
    std::string get_description() const { return get_type_description(getType()); }

    /// Construct a GenericType given an TypeID
#ifndef SWIG
    static GenericType from_type(TypeID type);

    ///@{
    /// Implicit typecasting
    operator bool() const { return to_bool();}
    operator casadi_int() const { return to_int();}
    operator int() const { return to_int();}
    operator double() const { return to_double();}
    operator std::string() const { return to_string();}
    operator std::vector<bool>() const { return to_bool_vector();}
    operator std::vector<casadi_int>() const { return to_int_vector();}
    operator std::vector<int>() const;
    operator std::vector<std::vector<casadi_int> >() const { return to_int_vector_vector();}
    operator std::vector<std::vector<int> >() const;
    operator std::vector<double>() const { return to_double_vector();}
    operator std::vector< std::vector<double> >() const {
      return to_double_vector_vector();
    }
    operator std::vector<std::string>() const { return to_string_vector();}
    operator const Function&() const { return as_function();}
    operator const std::vector<Function>&() const { return as_function_vector();}
    operator const Dict&() const { return as_dict();}
    ///@}

    bool can_cast_to(TypeID other) const;
    bool can_cast_to(const GenericType& other) const { return can_cast_to(other.getType()) ;}
#endif

#if !(defined(SWIG) && !defined(SWIGXML))
    // Get type of object
    TypeID getType() const;
#endif

    ///@{
    /** \brief Check if a particular type */
    bool is_bool() const;
    bool is_int() const;
    bool is_double() const;
    bool is_string() const;
    bool is_empty_vector() const;
    bool is_int_vector() const;
    bool is_int_vector_vector() const;
    bool is_double_vector() const;
    bool is_double_vector_vector() const;
    bool is_bool_vector() const;
    bool is_string_vector() const;
    bool is_dict() const;
    bool is_function() const;
    bool is_function_vector() const;
    bool is_void_pointer() const;
    ///@}

    ///@{
    /** \brief Cast to the internal type */
    const bool& as_bool() const;
    const casadi_int& as_int() const;
    const double& as_double() const;
    const std::string& as_string() const;
    const std::vector<casadi_int>& as_int_vector() const;
    const std::vector<casadi_int>& as_bool_vector() const;
    const std::vector<std::vector<casadi_int> >& as_int_vector_vector() const;
    const std::vector<double>& as_double_vector() const;
    const std::vector< std::vector<double> >& as_double_vector_vector() const;
    const std::vector<std::string>& as_string_vector() const;
    const Dict& as_dict() const;
    const Function& as_function() const;
    const std::vector<Function>& as_function_vector() const;
    void* const & as_void_pointer() const;
    ///@}

    ///@{
    //! \brief Convert to a type
    bool to_bool() const;
    casadi_int to_int() const;
    double to_double() const;
    std::string to_string() const;
    std::vector<casadi_int> to_int_vector() const;
    std::vector<bool> to_bool_vector() const;
    std::vector< std::vector<casadi_int> > to_int_vector_vector() const;
    std::vector<double> to_double_vector() const;
    std::vector< std::vector<double> > to_double_vector_vector() const;
    std::vector<std::string> to_string_vector() const;
    Dict to_dict() const;
    Function to_function() const;
    std::vector<Function> to_function_vector() const;
    void* to_void_pointer() const;
    std::vector<int> to_int_type_vector() const;
    ///@}

    //! \brief Equality
    bool operator==(const GenericType& op2) const;
    bool operator!=(const GenericType& op2) const;
#endif // SWIG
  };

  /// C++ equivalent of Python's dict or MATLAB's struct
  typedef GenericType::Dict Dict;

} // namespace casadi


#endif // CASADI_GENERIC_TYPE_HPP
