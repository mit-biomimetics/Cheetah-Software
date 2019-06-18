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


#ifndef CASADI_POLYNOMIAL_HPP
#define CASADI_POLYNOMIAL_HPP

# include "printable.hpp"

namespace casadi {

  /** \brief Helper class for differentiating and integrating polynomials
      \author Joel Andersson
      \date 2014
  */
  class CASADI_EXPORT Polynomial : public Printable<Polynomial> {
  public:
    /// Floating point type
    typedef double casadi_real;

    /// Construct a constant polynomial
    Polynomial(casadi_real scalar=1);

    /// Construct a linear polynomial
    Polynomial(casadi_real p0, casadi_real p1);

    /// Construct a quadratic polynomial
    Polynomial(casadi_real p0, casadi_real p1, casadi_real p2);

    /// Construct a cubic polynomial
    Polynomial(casadi_real p0, casadi_real p1, casadi_real p2, casadi_real p3);

    /// Construct from a vector of polynomial coefficients
    template<typename T>
    Polynomial(const std::vector<T>& coeff) : p_(coeff.begin(), coeff.end()) {}

    /// Evaluate numerically
    template<typename T>
    T operator()(const T& x) const {
      auto it = p_.rbegin();
      T ret = *it++;
      while (it!=p_.rend()) {
        ret *= x;
        ret += *it++;
      }
      return ret;
    }

    /// Degree of the polynomial
    casadi_int degree() const;

    /// Get scalar value (error if degree()!=0)
    casadi_real scalar() const;

    /// Create a new polynomial for the derivative
    Polynomial derivative() const;

    /// Create a new polynomial for the anti-derivative (primitive function)
    Polynomial anti_derivative() const;

    /// Remove excess zeros
    void trim();

    /// Readable name of the class
    std::string type_name() const {return "Polynomial";}

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    // Add
    Polynomial operator+(const Polynomial& b) const;

    // Add (in-place)
    Polynomial& operator+=(const Polynomial& b);

    // Subtract
    Polynomial operator-(const Polynomial& b) const;

    // Subtract (in-place)
    Polynomial& operator-=(const Polynomial& b);

    // Multiply
    Polynomial operator*(const Polynomial& b) const;

    // Multiply (in-place)
    Polynomial& operator*=(const Polynomial& b);

    // Divide by constant
    Polynomial operator/(casadi_real b) const;

    // Divide by constant (in-place)
    Polynomial& operator/=(casadi_real b);


  protected:
    std::vector<casadi_real> p_;
  };

} // namespace casadi


#endif // CASADI_POLYNOMIAL_HPP
