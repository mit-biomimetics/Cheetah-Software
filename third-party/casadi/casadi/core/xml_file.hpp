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


#ifndef CASADI_XML_FILE_HPP
#define CASADI_XML_FILE_HPP

#include "xml_file.hpp"
#include "shared_object.hpp"
#include "xml_node.hpp"
#include "printable.hpp"

namespace casadi {

  /** Forward declaration of internal class */
  class XmlFileInternal;

  /** \brief XML parser
      Can be used for parsing XML files into CasADi data structures.

      \author Joel Andersson
      \date 2014
   */
  class CASADI_EXPORT XmlFile
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<XmlFile>) {
  public:
    /** \brief Get type name */
    static std::string type_name() {return "XmlFile";}

    // Default constructor
    XmlFile();

    // Constructor
    XmlFile(const std::string& name);

    // Destructor
    ~XmlFile();

    /// Load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

#ifndef SWIG
    /** \brief  Access functions of the node */
    XmlFileInternal* operator->();

    /** \brief  Const access functions of the node */
    const XmlFileInternal* operator->() const;

    // Parse an XML file
    XmlNode parse(const std::string& filename);
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_XML_FILE_HPP
