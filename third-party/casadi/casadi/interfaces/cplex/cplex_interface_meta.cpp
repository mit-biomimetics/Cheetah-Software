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


      #include "cplex_interface.hpp"
      #include <string>

      const std::string casadi::CplexInterface::meta_doc=
      "\n"
"Interface to Cplex solver for sparse Quadratic Programs\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| barrier_maxiter | OT_INT      | 2.100e+09       | Maximum number  |\n"
"|                 |                 |                 | of barrier      |\n"
"|                 |                 |                 | iterations.     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| convex          | OT_BOOL      | true            | Indicates if    |\n"
"|                 |                 |                 | the QP is       |\n"
"|                 |                 |                 | convex or not   |\n"
"|                 |                 |                 | (affects only   |\n"
"|                 |                 |                 | the barrier     |\n"
"|                 |                 |                 | method).        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| dep_check       | OT_STRING       | \"off\"           | Detect          |\n"
"|                 |                 |                 | redundant       |\n"
"|                 |                 |                 | constraints. (a |\n"
"|                 |                 |                 | utomatic:-1|off |\n"
"|                 |                 |                 | :0|begin:1|end: |\n"
"|                 |                 |                 | 2|both:3)       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| dump_filename   | OT_STRING       | \"qp.dat\"        | The filename to |\n"
"|                 |                 |                 | dump to.        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| dump_to_file    | OT_BOOL      | false           | Dumps QP to     |\n"
"|                 |                 |                 | file in CPLEX   |\n"
"|                 |                 |                 | format.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| qp_method       | OT_STRING       | \"automatic\"     | Determines      |\n"
"|                 |                 |                 | which CPLEX     |\n"
"|                 |                 |                 | algorithm to    |\n"
"|                 |                 |                 | use. (automatic |\n"
"|                 |                 |                 | |primal_simplex |\n"
"|                 |                 |                 | |dual_simplex|n |\n"
"|                 |                 |                 | etwork|barrier| |\n"
"|                 |                 |                 | sifting|concurr |\n"
"|                 |                 |                 | ent|crossover)  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| simplex_maxiter | OT_INT      | 2.100e+09       | Maximum number  |\n"
"|                 |                 |                 | of simplex      |\n"
"|                 |                 |                 | iterations.     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| tol             | OT_DOUBLE         | 0.000           | Tolerance of    |\n"
"|                 |                 |                 | solver          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| warm_start      | OT_BOOL      | false           | Use warm start  |\n"
"|                 |                 |                 | with simplex    |\n"
"|                 |                 |                 | methods         |\n"
"|                 |                 |                 | (affects only   |\n"
"|                 |                 |                 | the simplex     |\n"
"|                 |                 |                 | methods).       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
