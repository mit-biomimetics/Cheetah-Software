#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
# -*- coding: utf-8 -*-
from casadi import *

# Tests the finite difference implementation in CasADi
# Joel Andersson, UW Madison 2017

# Try with different scaling factors for x
a_test = [1e-10, 1e-4, 1., 1e4, 1e10]

for a in a_test:
    print()
    print("*** Testing a=%g ***" % a)

    # The following function is defined for the interval [0,1] and gives NaN outside it
    x = SX.sym('x');
    r = if_else(x<=1, if_else(x>=0, 2*x*a+cos(x*a), log(x)), log(1-x))/a

    # Test points in the domain, on the boundary and outside it
    x_test = [-1., 0., 0.4, 1., 1.4]
    print("%30s" % "x0", end='')
    for x0 in x_test: print("%15g " % x0, end='')
    print()

    # Try all supported schemes
    fd_methods = []
    fd_methods.append(('analytical', 0))
    fd_methods.append(('forward', 0))
    fd_methods.append(('backward', 0))
    fd_methods.append(('central', 0))
    fd_methods.append(('central', 1))
    fd_methods.append(('central', 2))
    fd_methods.append(('central', 5))
    fd_methods.append(('central', 20))
    fd_methods.append(('smoothing', 0))
    fd_methods.append(('smoothing', 1))
    fd_methods.append(('smoothing', 2))
    fd_methods.append(('smoothing', 5))
    fd_methods.append(('smoothing', 20))
    for (fd_method, h_iter) in fd_methods:
        descr = fd_method
        if h_iter!=0: descr += ' with ' + str(h_iter) + ' h updates'
        print("%30s" % descr, end='')
        # Construct function and differentiate
        opts = dict()
        if fd_method!='analytical':
            opts["enable_fd"]=True # enable finite differencing
            opts["enable_forward"]=False # disable forward mode AD
            opts["enable_reverse"]=False # disable reverse mode AD
            opts["enable_jacobian"]=False # disable AD by calculating full Jacobian
            opts["fd_method"]=fd_method # specify FD scheme
            if fd_method in ['central','smoothing']: opts["fd_options"] = dict(h_iter=h_iter)
        f = Function('f', [x], [r], ['x'], ['r'], opts)
        fwd_f = f.forward(1)
        # Evaluate
        for x0 in x_test:
            r0 = f(x0)
            d = fwd_f(x0, r0, 1.)
            print("%15g " % float(d), end='')
        print()
