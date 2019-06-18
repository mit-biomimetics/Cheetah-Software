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


#ifndef CASADI_CONSTANT_SX_HPP
#define CASADI_CONSTANT_SX_HPP

#include "sx_node.hpp"
#include <cassert>

/// \cond INTERNAL

// Cashing of constants requires a map
#include <unordered_map>
#define CACHING_MAP std::unordered_map

namespace casadi {

/** \brief Represents a constant SX
  \author Joel Andersson
  \date 2010
*/
class ConstantSX : public SXNode {
public:

// Destructor
~ConstantSX() override {}

// Class name
std::string class_name() const override {return "ConstantSX";}

/** \brief  Get the value must be defined */
double to_double() const override = 0;

/** \brief  Properties */
bool is_constant() const override { return true; }

/** \brief  Get the operation */
casadi_int op() const override { return OP_CONST;}

/** \brief Check if two nodes are equivalent up to a given depth */
bool is_equal(const SXNode* node, casadi_int depth) const override {
  const ConstantSX* n = dynamic_cast<const ConstantSX*>(node);
  return n && n->to_double()==to_double();
}

protected:

/** \brief  Print expression */
std::string print(const std::string& arg1, const std::string& arg2) const override {
   std::stringstream ss;
   ss << to_double();
   return ss.str();
 }

};

/** \brief  DERIVED CLASSES */

/** \brief  Represents a constant real SX
  \author Joel Andersson
  \date 2010
*/
class RealtypeSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit RealtypeSX(double value) : value(value) {}

  public:

    /// Destructor
    ~RealtypeSX() override {
      size_t num_erased = cached_constants_.erase(value);
      assert(num_erased==1);
      (void)num_erased;
    }

    /// Static creator function (use instead of constructor)
    inline static RealtypeSX* create(double value) {
      // Try to find the constant
      CACHING_MAP<double, RealtypeSX*>::iterator it = cached_constants_.find(value);

      // If not found, add it,
      if (it==cached_constants_.end()) {
        // Allocate a new object
        RealtypeSX* n = new RealtypeSX(value);

        // Add to hash_table
        cached_constants_.insert(it, std::make_pair(value, n));

        // Return it to caller
        return n;
      } else { // Else, returned the object
        return it->second;
      }
    }

    ///@{
    /** \brief  Get the value */
    double to_double() const override { return value;}
    casadi_int to_int() const override { return static_cast<casadi_int>(value);}
    ///@}

    bool is_almost_zero(double tol) const override { return fabs(value)<=tol; }

  protected:
    /** \brief Hash map of all constants currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static CACHING_MAP<double, RealtypeSX*> cached_constants_;

    /** \brief  Data members */
    double value;
};


/** \brief  Represents a constant integer SX
  \author Joel Andersson
  \date 2010
*/
class IntegerSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit IntegerSX(casadi_int value) : value(static_cast<int>(value)) {
      casadi_assert(value<=std::numeric_limits<int>::max() &&
                    value>=std::numeric_limits<int>::min(), "Integer overflow");
    }

  public:

    /// Destructor
    ~IntegerSX() override {
      size_t num_erased = cached_constants_.erase(value);
      assert(num_erased==1);
      (void)num_erased;
    }

    /// Static creator function (use instead of constructor)
    inline static IntegerSX* create(casadi_int value) {
      // Try to find the constant
      CACHING_MAP<casadi_int, IntegerSX*>::iterator it = cached_constants_.find(value);

      // If not found, add it,
      if (it==cached_constants_.end()) {
        // Allocate a new object
        IntegerSX* n = new IntegerSX(value);

        // Add to hash_table
        cached_constants_.insert(it, std::make_pair(value, n));

        // Return it to caller
        return n;
      } else { // Else, returned the object
        return it->second;
      }
    }

    ///@{
    /** \brief  evaluate function */
    double to_double() const override {  return static_cast<double>(value); }
    casadi_int to_int() const override {  return static_cast<casadi_int>(value); }
    ///@}

    /** \brief  Properties */
    bool is_integer() const override { return true; }

  protected:

    /** \brief Hash map of all constants currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static CACHING_MAP<casadi_int, IntegerSX*> cached_constants_;

    /** \brief  Data members */
    int value;
};

/** \brief  Represents a zero SX
  \author Joel Andersson
  \date 2010
*/
class ZeroSX : public ConstantSX {
public:

  ~ZeroSX() override {}
  explicit ZeroSX() {}

  ///@{
  /** \brief  Get the value */
  double to_double() const override { return 0;}
  casadi_int to_int() const override { return 0;}
  ///@}

  ///@{
  /** \brief  Properties */
  bool is_integer() const override { return true; }
  bool is_zero() const override { return true; }
  bool is_almost_zero(double tol) const override { return true; }
  ///@}
};


/** \brief  Represents a one SX
  \author Joel Andersson
  \date 2010
*/
class OneSX : public ConstantSX {
public:

  explicit OneSX() {}
  ~OneSX() override {}

  /** \brief  Get the value */
  double to_double() const override { return 1;}
  casadi_int to_int() const override { return 1;}

  /** \brief  Properties */
  bool is_integer() const override { return true; }
  bool is_one() const override { return true; }

};


/** \brief  Represents a minus one SX
  \author Joel Andersson
  \date 2010
*/
class MinusOneSX : public ConstantSX {
public:

  explicit MinusOneSX() {}
  ~MinusOneSX() override {}

  ///@{
  /** \brief  Get the value */
  double to_double() const override { return -1;}
  casadi_int to_int() const override { return -1;}
  ///@}

  ///@{
  /** \brief  Properties */
  bool is_integer() const override { return true; }
  bool is_minus_one() const override { return true; }
  ///@}

};


/** \brief  Represents an infinity SX
  \author Joel Andersson
  \date 2010
*/
class InfSX : public ConstantSX {
public:

  explicit InfSX() {}
  ~InfSX() override {}

  /** \brief  Get the value */
  double to_double() const override { return std::numeric_limits<double>::infinity();}

  /** \brief  Properties */
  bool is_inf() const override { return true; }

};


/** \brief  Represents a minus infinity SX
  \author Joel Andersson
  \date 2010
*/
class MinusInfSX : public ConstantSX {
public:

  explicit MinusInfSX() {}
  ~MinusInfSX() override {}

  /** \brief  Get the value */
  double to_double() const override { return -std::numeric_limits<double>::infinity();}

  /** \brief  Properties */
  bool is_minus_inf() const override { return true; }

};


/** \brief  Represents a not-a-number SX
  \author Joel Andersson
  \date 2010
*/
class NanSX : public ConstantSX {
public:

  explicit NanSX() {this->count++;}
  ~NanSX() override {this->count--;}

  /** \brief  Get the value */
  double to_double() const override { return std::numeric_limits<double>::quiet_NaN();}

  /** \brief  Properties */
  bool is_nan() const override { return true; }

};

} // namespace casadi
/// \endcond

#endif // CASADI_CONSTANT_SX_HPP
