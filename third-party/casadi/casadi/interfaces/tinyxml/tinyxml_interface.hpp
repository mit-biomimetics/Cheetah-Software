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


#ifndef CASADI_TINYXML_INTERFACE_HPP
#define CASADI_TINYXML_INTERFACE_HPP

/** \defgroup plugin_XmlFile_tinyxml
 * XmlFile using TinyXml
*/

/** \pluginsection{XmlFile,tinyxml} */

/// \cond INTERNAL
#include <tinyxml.h>
#include "casadi/core/xml_file_internal.hpp"
#include <casadi/interfaces/tinyxml/casadi_xmlfile_tinyxml_export.h>

namespace casadi {

  /** \brief \pluginbrief{XmlFile,tinyxml}
   * @copydoc XmlFile_doc
   * @copydoc plugin_XmlFile_tinyxml
   */
  class CASADI_XMLFILE_TINYXML_EXPORT TinyXmlInterface : public XmlFileInternal {
  public:

    // Create an XML file
    TinyXmlInterface();

    /** \brief  Create a new XmlFile */
    static XmlFileInternal* creator()
    { return new TinyXmlInterface();}

    // Get name of the plugin
    const char* plugin_name() const override { return "tinyxml";}

    // Get name of the class
    std::string class_name() const override { return "TinyXmlInterface";}

    // Parse an XML file
    XmlNode parse(const std::string& filename) override;

    // Parse an XML tree
    XmlNode addNode(TiXmlNode* n);

    // Destructor
    ~TinyXmlInterface() override;

    /// A documentation string
    static const std::string meta_doc;

    /// XML file
    TiXmlDocument doc_;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_TINYXML_INTERFACE_HPP
