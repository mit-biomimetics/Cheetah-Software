# QDLDL
A free LDL factorisation routine for quasi-definite linear systems: `Ax=b`

[![Build Status](https://travis-ci.org/oxfordcontrol/qdldl.svg?branch=master)](https://travis-ci.org/oxfordcontrol/qdldl)
[![Build status](https://ci.appveyor.com/api/projects/status/ns4br7v6y3i5stai/branch/master?svg=true)](https://ci.appveyor.com/project/bstellato/qdldl-8q1mv/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/oxfordcontrol/qdldl/badge.svg)](https://coveralls.io/github/oxfordcontrol/qdldl)


## Getting started
To start using QDLDL, first clone the repository

```bash
git clone https://github.com/oxfordcontrol/qdldl.git
```

### Build

To build QDLDL, you need to install [cmake](https://cmake.org/) and run

```bash
mkdir build/
cd build/
cmake ..
make
```

This will generate an `out/` folder with contents:

- `qdldl_example`: a **code example** from [`examples/c/example.c`](./examples/c/example.c)
- `libqdldl`: a static and a dynamic versions of the library.

You can include an addition option `-DUNITTESTS=ON` when calling `cmake`, which will result in an additional executable `qdldl_tester` being built in the `out/` folder to test QDLDL on a variety of problems, including those with rank deficient or otherwise ill-formatted inputs.

**N.B.** All files will have file extensions appropriate to your operating system.


### Install/Uninstall

To install (uninstall) the libraries and headers you can simply run `make install` (`make uninstall`) after running the `cmake` command above.


## Calling QDLDL

### Main API

The QDLDL API consists of 5 functions documented in [`include/qdldl.h`](./include/qdldl.h).
For more details and a working example see [`examples/c/example.c`](./examples/c/example.c).

**N.B.** There is **no memory allocation** performed in these routines. The user is assumed to have the working vectors already allocated.

Here is a brief summary.

* `QDLDL_etree`: compute the elimination tree for the quasidefinite matrix factorization `A = LDL'`
* `QDLDL_factor`: return the factors `L`, `D` and `Dinv = 1./D`
* `QDLDL_solve`: solves the linear system `LDL'x = b`
* `QDLDL_Lsolve`: solves `Lx = b`
* `QDLDL_Ltsolve`: solves `L'x = b`

In the above function calls the matrices `A` and `L` are stored in compressed sparse column (CSC) format.   The matrix `A` is assumed to be symmetric and only the upper triangular portion of A should be passed to the API.   The factor `L` is lower triangular with implicit ones on the diagonal (i.e. the diagonal of L is not stored as part of the CSC formatted data.)

The matrices `D` and `Dinv` are both diagonal matrices, with the diagonal values stored in an array.

The matrix input `A` should be quasidefinite.   The API provides some (non-comprehensive) error checking to protect against non-quasidefinite or non-upper triangular inputs.

### Custom types for integer, floats and booleans
QDLDL uses its own internal types for integers, floats and booleans (`QDLDL_int, QDLDL_float, QDLDL_bool`. They can be specified using the cmake options:

- `DFLOAT` (default false): uses float numbers instead of doubles
- `DLONG` (default true): uses long integers for indexing (for large matrices)

The `QDLDL_bool` is internally defined as `unsigned char`.


## Linking QDLDL

### Basic Example
A basic example appears in [`examples/c/example.c`](./examples/c/example.c) and is compiled using cmake and the `CMakeLists.txt` file in the root folder.


### Including in a cmake project

You can include QDLDL in a cmake project `foo` by adding the subdirectory as
```
# Add project
add_subdirectory(qdldl)
```

QDLDL can be linked using a static or dynamic linker
```
# Link static library
target_link_libraries (foo qdldlstatic)

# Link shared library
target_link_libraries (foo qdldl)
```
for dynamic linking the shared library should be available in your path.

There is also the option to include QDLDL as an object library in your project.
The current `CMakeLists.txt` file creates an object library called `qdldlobject`.
This can be added to your project by adding it after your sources.
For example, when creating a library `foo` you can add

```
add_library(foo foo.c foo.h $<TARGET_OBJECTS:qdldlobject>)
```
for more details see the [cmake documentation](https://cmake.org/cmake/help/latest/command/add_library.html#object-libraries).


## The algorithm

The algorithm is an independent implementation of the elimination tree and factorisation procedures outlined in

> T. A Davis. [Algorithm 849: a concise sparse cholesky factorization package](https://dl.acm.org/citation.cfm?id=1114277). ACM Trans. Math. Softw., 31(4):587–591, 2005.


## Credits

- [Paul Goulart](http://users.ox.ac.uk/~engs1373/): main development
- [Bartolomeo Stellato](https://stellato.io/): code refactoring and testing
- [Goran Banjac](http://control.ee.ethz.ch/~gbanjac/)
