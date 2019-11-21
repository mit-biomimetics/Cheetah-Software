# The Operator Splitting QP Solver

[![Build status of the master branch on Linux/OSX](https://img.shields.io/travis/oxfordcontrol/osqp/master.svg?label=Linux%20%2F%20OSX%20build)](https://travis-ci.org/oxfordcontrol/osqp)
[![Build status of the master branch on Windows](https://img.shields.io/appveyor/ci/bstellato/osqp/master.svg?label=Windows%20build)](https://ci.appveyor.com/project/bstellato/osqp/branch/master)
[![Code coverage](https://coveralls.io/repos/github/oxfordcontrol/osqp/badge.svg?branch=master)](https://coveralls.io/github/oxfordcontrol/osqp?branch=master)
![License](https://img.shields.io/badge/License-Apache%202.0-brightgreen.svg)


![PyPI - downloads](https://img.shields.io/pypi/dm/osqp.svg?label=Pypi%20downloads)
![Conda - downloads](https://img.shields.io/conda/dn/conda-forge/osqp.svg?label=Conda%20downloads)

[**Join our forum**](https://groups.google.com/forum/#!forum/osqp) for any questions related to the solver!

**The documentation** is available at [**osqp.org**](https://osqp.org/)

The OSQP (Operator Splitting Quadratic Program) solver is a numerical optimization package for solving problems in the form
```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```

where `x in R^n` is the optimization variable. The objective function is defined by a positive semidefinite matrix `P in S^n_+` and vector `q in R^n`. The linear constraints are defined by matrix `A in R^{m x n}` and vectors `l in R^m U {-inf}^m`, `u in R^m U {+inf}^m`.


The latest version is `0.5.0`.

## Citing OSQP

If you are using OSQP for your work, we encourage you to

* [Cite the related papers](https://osqp.org/citing/),
* Put a star on this repository.

**We are looking forward to hearing your success stories with OSQP!** Please [share them with us](mailto:bartolomeo.stellato@gmail.com).


## Bug reports and support

Please report any issues via the [Github issue tracker](https://github.com/oxfordcontrol/osqp/issues). All types of issues are welcome including bug reports, documentation typos, feature requests and so on.


## Numerical benchmarks
Numerical benchmarks against other solvers are available [here](https://github.com/oxfordcontrol/osqp_benchmarks).

