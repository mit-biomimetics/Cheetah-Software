blockSQP -- Sequential quadratic programming for problems with
            block-diagonal Hessian matrix.
Copyright (c) 2012-2015 Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>


Introduction
============
blockSQP is a sequential quadratic programming method for finding local solutions
of nonlinear, nonconvex optimization problems. It is particularly suited for
---but not limited to---problems whose Hessian matrix has block-diagonal
structure such as problems arising from direct multiple shooting
parameterizations of optimal control or optimum experimental design problems.

blockSQP has been developed around the quadratic programming solver
qpOASES to solve the quadratic subproblems. Gradients of the objective
and the constraint functions must be supplied by the user. Second derivatives
are approximated by a combination of SR1 and BFGS updates. Global convergence
is promoted by the filter line search of Waechter and Biegler that can also
handle indefinite Hessian approximations.


Installation
============
* Download and install qpOASES from

  https://projects.coin-or.org/qpOASES .

  It is recommended to use at least release 3.2.0.
  Alternatively, check out revision 155 from the qpOASES subversion
  repository that is located at

  https://projects.coin-or.org/svn/qpOASES/trunk/ .

  For best performance it is strongly recommended to install the sparse
  solver MA57 from HSL as described in the qpOASES manual, Sec. 2.2.

* In the blockSQP main directory, open `makefile` and set `QPOASESDIR`
  to the correct location of the qpOASES installation.

* Compile blockSQP by calling `make`. This should produce a shared
  library `libblockSQP.so` in `lib/`, as well as executable example problems
  in the `examples/` folder.


Documentation
=============
A user's manual for blockSQP is available under `doc/manual.pdf`.
There is also Doxygen source code documentation available that can
be created by calling

  doxygen doxyfileBLOCKSQP

within the `doc/` directory.


Licensing
=========
blockSQP is published under the very permissive zlib free software
license which should allow you to use the software wherever you need.

This is the full license text (zlib license):

    blockSQP -- Sequential quadratic programming for problems with
                block-diagonal Hessian matrix.
    Copyright (c) 2012-2015 Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>

    This software is provided 'as-is', without any express or implied
    warranty. In no event will the authors be held liable for any
    damages arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute
    it freely, subject to the following restrictions:

        1. The origin of this software must not be misrepresented;
        you must not claim that you wrote the original software.
        If you use this software in a product, an acknowledgment in the
        product documentation would be appreciated but is not required.

        2. Altered source versions must be plainly marked as such,
        and must not be misrepresented as being the original software.

        3. This notice may not be removed or altered from any source
        distribution.
