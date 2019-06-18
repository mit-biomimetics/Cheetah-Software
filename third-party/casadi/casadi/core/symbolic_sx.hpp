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


#ifndef SYMBOLIC_SXElem_HPP
#define SYMBOLIC_SXElem_HPP

#include "sx_node.hpp"
/// \cond INTERNAL

namespace casadi {

/** \brief Represents a scalar symbolic expression
  \author Joel Andersson
  \date 2010
  A regular user is not supposed to work with this Node class.
  This user can call SX(name) instead.
*/
class SymbolicSX : public SXNode {
public:
  explicit SymbolicSX(const std::string &name) : name_(name) {}
  ~SymbolicSX() override {}

  bool is_symbolic() const override { return true; }

  const std::string& name() const override { return name_; }

  /** \brief  Get the operation */
  casadi_int op() const override { return OP_PARAMETER;}

  bool is_op(casadi_int op) const override { return op==OP_PARAMETER; }

  /** \brief  Name */
  std::string name_;

  // Class name
  std::string class_name() const override {return "SymbolicSX";}

  /** \brief  Print expression */
  std::string print(const std::string& arg1, const std::string& arg2) const override {
    return name_;
  }
};

} // namespace casadi
/// \endcond
#endif // SYMBOLIC_SXElem_HPP
