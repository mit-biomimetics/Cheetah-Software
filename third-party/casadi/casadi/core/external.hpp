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


#ifndef CASADI_EXTERNAL_HPP
#define CASADI_EXTERNAL_HPP

#include "function.hpp"
#include "importer.hpp"

namespace casadi {

  /** \brief  Load an external function
   * File name is assumed to be ./<f_name>.so
   */
  CASADI_EXPORT Function external(const std::string& name, const Dict& opts=Dict());

  /** \brief  Load an external function
   * File name given
   */
  CASADI_EXPORT Function external(const std::string& name, const std::string& bin_name,
                                  const Dict& opts=Dict());

  /** \brief  Load a just-in-time compiled external function
   * File name given
   */
  CASADI_EXPORT Function external(const std::string& name, const Importer& compiler,
                                  const Dict& opts=Dict());

} // namespace casadi

#endif // CASADI_EXTERNAL_HPP
