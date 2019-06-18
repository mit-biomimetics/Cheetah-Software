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


#include "casadi_meta.hpp"
#include <casadi/config.h>

namespace casadi {
  const char* CasadiMeta::version() { return CASADI_VERSION_STRING;}

  const char* CasadiMeta::git_revision() { return CASADI_GIT_REVISION;}

  const char* CasadiMeta::git_describe() { return CASADI_GIT_DESCRIBE;}

  const char* CasadiMeta::feature_list() { return CASADI_FEATURE_LIST;}

  const char* CasadiMeta::build_type() { return CASADI_BUILD_TYPE;}

  const char* CasadiMeta::compiler_id() { return CASADI_COMPILER_ID;}

  const char* CasadiMeta::compiler() { return CASADI_COMPILER;}

  const char* CasadiMeta::compiler_flags() { return CASADI_COMPILER_FLAGS;}

  const char* CasadiMeta::modules() { return CASADI_MODULES;}

  const char* CasadiMeta::plugins() { return CASADI_PLUGINS;}

  const char* CasadiMeta::install_prefix() { return CASADI_INSTALL_PREFIX;}

}  // namespace casadi
