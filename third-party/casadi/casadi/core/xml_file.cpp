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

  XmlFile::XmlFile() {
  }

  XmlFile::XmlFile(const std::string& name) {
    own(XmlFileInternal::getPlugin(name).creator());
  }

  XmlFile::~XmlFile() {
  }

  const XmlFileInternal* XmlFile::operator->() const {
    return static_cast<const XmlFileInternal*>(SharedObject::operator->());
  }

  XmlFileInternal* XmlFile::operator->() {
    return static_cast<XmlFileInternal*>(SharedObject::operator->());
  }

  XmlNode XmlFile::parse(const std::string& filename) {
    return (*this)->parse(filename);
  }

  void XmlFile::load_plugin(const std::string& name) {
    XmlFileInternal::load_plugin(name);
  }

  std::string XmlFile::doc(const std::string& name) {
    return XmlFileInternal::getPlugin(name).doc;
  }

} // namespace casadi
