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


#include "xml_file_internal.hpp"

using namespace std;
namespace casadi {

  XmlFileInternal::XmlFileInternal() {
  }

  XmlFileInternal::~XmlFileInternal() {
  }

  std::map<std::string, XmlFileInternal::Plugin> XmlFileInternal::solvers_;

  const std::string XmlFileInternal::infix_ = "xmlfile";

  void XmlFileInternal::disp(ostream &stream, bool more) const {
  }

  XmlNode XmlFileInternal::parse(const std::string& filename) {
    casadi_error("parse not defined for " + class_name());
    return XmlNode();
  }

} // namespace casadi
