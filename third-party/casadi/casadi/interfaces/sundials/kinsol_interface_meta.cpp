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


      #include "kinsol_interface.hpp"
      #include <string>

      const std::string casadi::KinsolInterface::meta_doc=
      "\n"
"KINSOL interface from the Sundials suite\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| abstol          | OT_DOUBLE         | 0.000           | Stopping        |\n"
"|                 |                 |                 | criterion       |\n"
"|                 |                 |                 | tolerance       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| disable_interna | OT_BOOL      | false           | Disable KINSOL  |\n"
"| l_warnings      |                 |                 | internal        |\n"
"|                 |                 |                 | warning         |\n"
"|                 |                 |                 | messages        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| exact_jacobian  | OT_BOOL      | true            |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| f_scale         | OT_DOUBLEVECTOR   |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| iterative_solve | OT_STRING       | \"gmres\"         | gmres|bcgstab|t |\n"
"| r               |                 |                 | fqmr            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| linear_solver_t | OT_STRING       | \"dense\"         | dense|banded|it |\n"
"| ype             |                 |                 | erative|user_de |\n"
"|                 |                 |                 | fined           |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| lower_bandwidth | OT_INT      |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| max_iter        | OT_INT      | 0               | Maximum number  |\n"
"|                 |                 |                 | of Newton       |\n"
"|                 |                 |                 | iterations.     |\n"
"|                 |                 |                 | Putting 0 sets  |\n"
"|                 |                 |                 | the default     |\n"
"|                 |                 |                 | value of        |\n"
"|                 |                 |                 | KinSol.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| max_krylov      | OT_INT      | 0               |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| pretype         | OT_STRING       | \"none\"          | (none|left|righ |\n"
"|                 |                 |                 | t|both)         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| strategy        | OT_STRING       | \"none\"          | Globalization   |\n"
"|                 |                 |                 | strategy (none| |\n"
"|                 |                 |                 | linesearch)     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| u_scale         | OT_DOUBLEVECTOR   |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| upper_bandwidth | OT_INT      |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| use_preconditio | OT_BOOL      | false           | precondition an |\n"
"| ner             |                 |                 | iterative       |\n"
"|                 |                 |                 | solver          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
">List of available monitors\n"
"\n"
"+-----------+\n"
"|    Id     |\n"
"+===========+\n"
"| eval_djac |\n"
"+-----------+\n"
"| eval_f    |\n"
"+-----------+\n"
"\n"
"\n"
"\n"
"\n"
;
