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


#include "importer.hpp"
#include "importer_internal.hpp"

using namespace std;
namespace casadi {

  Importer::Importer() {
  }

  Importer::Importer(const std::string& name,
                           const std::string& compiler,
                           const Dict& opts) {
    if (compiler=="none") {
      own(new ImporterInternal(name));
    } else if (compiler=="dll") {
      own(new DllLibrary(name));
    } else {
      own(ImporterInternal::getPlugin(compiler).creator(name));
    }
    (*this)->construct(opts);
  }

  ImporterInternal* Importer::operator->() {
    return static_cast<ImporterInternal*>(SharedObject::operator->());
  }

  const ImporterInternal* Importer::operator->() const {
    return static_cast<const ImporterInternal*>(SharedObject::operator->());
  }

  bool Importer::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const ImporterInternal*>(ptr)!=nullptr;
  }

  bool Importer::has_plugin(const std::string& name) {
    return ImporterInternal::has_plugin(name);
  }

  void Importer::load_plugin(const std::string& name) {
    ImporterInternal::load_plugin(name);
  }

  std::string Importer::doc(const std::string& name) {
    return ImporterInternal::getPlugin(name).doc;
  }

  std::string Importer::plugin_name() const {
    return (*this)->plugin_name();
  }

  bool Importer::has_function(const std::string& symname) const {
    return (*this)->has_function(symname);
  }

  signal_t Importer::get_function(const std::string& symname) {
    return (*this)->get_function(symname);
  }

  bool Importer::has_meta(const std::string& cmd, casadi_int ind) const {
    return (*this)->has_meta(cmd, ind);
  }

  std::string Importer::get_meta(const std::string& cmd, casadi_int ind) const {
    return (*this)->get_meta(cmd, ind);
  }

  bool Importer::inlined(const std::string& symname) const {
    return (*this)->inlined(symname);
  }

  std::string Importer::body(const std::string& symname) const {
    return (*this)->body(symname);
  }

} // namespace casadi
