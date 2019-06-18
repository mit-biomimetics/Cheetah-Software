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


      #include "slicot_expm.hpp"
      #include <string>

      const std::string casadi::SlicotExpm::meta_doc=
      "\n"
"An efficient solver for Discrete Periodic Lyapunov Equations using\n"
"SLICOT\n"
"\n"
"Uses Periodic Schur Decomposition ('psd') and does not assume positive\n"
"definiteness. Based on Periodic Lyapunov equations: some applications\n"
"and new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997.\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| linear_solver   | OT_STRING       | GenericType()   | User-defined    |\n"
"|                 |                 |                 | linear solver   |\n"
"|                 |                 |                 | class. Needed   |\n"
"|                 |                 |                 | for             |\n"
"|                 |                 |                 | sensitivities.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| linear_solver_o | OT_DICT   | GenericType()   | Options to be   |\n"
"| ptions          |                 |                 | passed to the   |\n"
"|                 |                 |                 | linear solver.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| psd_num_zero    | OT_REAL         | 0.000           | Numerical zero  |\n"
"|                 |                 |                 | used in         |\n"
"|                 |                 |                 | Periodic Schur  |\n"
"|                 |                 |                 | decomposition   |\n"
"|                 |                 |                 | with            |\n"
"|                 |                 |                 | slicot.This     |\n"
"|                 |                 |                 | option is       |\n"
"|                 |                 |                 | needed when     |\n"
"|                 |                 |                 | your systems    |\n"
"|                 |                 |                 | has Floquet     |\n"
"|                 |                 |                 | multiplierszero |\n"
"|                 |                 |                 | or close to     |\n"
"|                 |                 |                 | zero            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
">List of available stats\n"
"\n"
"+----------------+\n"
"|       Id       |\n"
"+================+\n"
"| t_linear_solve |\n"
"+----------------+\n"
"| t_psd          |\n"
"+----------------+\n"
"| t_total        |\n"
"+----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
