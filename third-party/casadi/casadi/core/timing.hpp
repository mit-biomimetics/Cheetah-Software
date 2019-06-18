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


#ifndef CASADI_TIMING_HPP
#define CASADI_TIMING_HPP

#include "generic_type.hpp"

#include <ctime>
#include <chrono>

namespace casadi {
  /// \cond INTERNAL

  /**
  Timer class


  FStats hack;
  hack.tic();
  ....
  hack.toc();

  */
  class CASADI_EXPORT FStats {
    private:
      /// Time point used for wall time computation
      std::chrono::time_point<std::chrono::high_resolution_clock> start_wall;

      /// Time point used for proc time computation
      std::clock_t start_proc;

      /// Time point used for wall time computation
      std::chrono::time_point<std::chrono::high_resolution_clock> stop_wall;

      /// Time point used for proc time computation
      std::clock_t stop_proc;

    public:
      /// Constructor
      FStats();

      /// Reset the statistics
      void reset();

      /// Start timing
      void tic();

      /// Stop timing
      void toc();

      /// Accumulated number of calls since last reset
      casadi_int n_call;

      /// Accumulated wall time [s] since last reset
      double t_wall;

      /// Accumulated proc time [s] since last reset
      double t_proc;
  };
/// \endcond
} // namespace casadi

#endif // CASADI_TIMING_HPP
