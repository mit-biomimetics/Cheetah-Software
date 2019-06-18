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


#ifndef CASADI_INTERRUPT_HPP
#define CASADI_INTERRUPT_HPP

#include <casadi/core/casadi_export.h>
#include "exception.hpp"

#include <iostream>
#include <fstream>

namespace casadi {

  /// \cond INTERNAL

  /**
   * \brief Takes care of user interrupts (Ctrl+C)
   *
   * This is an internal class.
   *
   *  \author Joris Gillis
   *  \date 2015
   */
  class CASADI_EXPORT InterruptHandler {
  private:
    /// No implementation - no instances are allowed of this class
    InterruptHandler();

    /// By default, report everything is okay
    static bool checkInterruptedDefault() {
      return false;
    }

  public:
    /// The routine that is used for checking interrupts
    static bool (*checkInterrupted)();

    /// Raises an error if an interrupt was captured.
    static void check() {
      if (checkInterrupted()) throw KeyboardInterruptException();
    }
  };

   /// \endcond INTERNAL

} // namespace casadi

#endif // CASADI_INTERRUPT_HPP
