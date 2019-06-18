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


      #include "collocation.hpp"
      #include <string>

      const std::string casadi::Collocation::meta_doc=
      "\n"
"Fixed-step implicit Runge-Kutta integrator ODE/DAE integrator based on\n"
"collocation schemes\n"
"\n"
"The method is still under development\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| collocation_sch | OT_STRING       | \"radau\"         | Collocation     |\n"
"| eme             |                 |                 | scheme (radau|l |\n"
"|                 |                 |                 | egendre)        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| implicit_solver | OT_STRING       | GenericType()   | An implicit     |\n"
"|                 |                 |                 | function solver |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| implicit_solver | OT_DICT   | GenericType()   | Options to be   |\n"
"| _options        |                 |                 | passed to the   |\n"
"|                 |                 |                 | NLP Solver      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| interpolation_o | OT_INT      | 3               | Order of the    |\n"
"| rder            |                 |                 | interpolating   |\n"
"|                 |                 |                 | polynomials     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| number_of_finit | OT_INT      | 20              | Number of       |\n"
"| e_elements      |                 |                 | finite elements |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
