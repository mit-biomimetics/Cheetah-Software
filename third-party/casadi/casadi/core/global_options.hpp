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


#ifndef CASADI_GLOBAL_OPTIONS_HPP
#define CASADI_GLOBAL_OPTIONS_HPP

#include <iostream>
#include <fstream>

#include <casadi/core/casadi_export.h>
#include "casadi/core/casadi_common.hpp"

namespace casadi {

  /**
  * \brief Collects global CasADi options
  *
  *
  * Note to developers:  \n
  *  - use sparingly. Global options are - in general - a rather bad idea \n
  *  - this class must never be instantiated. Access its static members directly \n
  *
  *  \author Joris Gillis
  *  \date 2012
  */
  class CASADI_EXPORT GlobalOptions {
    private:
      /// No instances are allowed
      GlobalOptions();
    public:

#ifndef SWIG
      /** \brief Indicates whether simplifications should be made on the fly.
      * e.g.   cos(-x) -> cos(x)
      * Default: true
      */
      static bool simplification_on_the_fly;

      static std::string casadipath;

      static bool hierarchical_sparsity;

      static casadi_int max_num_dir;

      static casadi_int start_index;

#endif //SWIG
      // Setter and getter for simplification_on_the_fly
      static void setSimplificationOnTheFly(bool flag) { simplification_on_the_fly = flag; }
      static bool getSimplificationOnTheFly() { return simplification_on_the_fly; }

      // Setter and getter for hierarchical_sparsity
      static void setHierarchicalSparsity(bool flag) { hierarchical_sparsity = flag; }
      static bool getHierarchicalSparsity() { return hierarchical_sparsity; }

      static void setCasadiPath(const std::string & path) { casadipath = path; }
      static std::string getCasadiPath() { return casadipath; }

      static void setMaxNumDir(casadi_int ndir) { max_num_dir=ndir; }
      static casadi_int getMaxNumDir() { return max_num_dir; }

  };

} // namespace casadi

#endif // CASADI_GLOBAL_OPTIONS_HPP
