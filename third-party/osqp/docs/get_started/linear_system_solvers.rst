.. _linear_system_solvers_installation :

Linear System Solvers
======================

The linear system solver is a core part of the OSQP algorithm.
Depending on the problem instance, different linear system solvers can greatly speedup or reduce the computation time of OSQP.
To set the preferred linear system solver, see :ref:`linear_system_solvers_setting`.

Dynamic shared library loading
------------------------------
OSQP dynamically loads the shared libraries related to the selected external solver. Thus, there is no need to link it at compile time.
The only requirement is that the shared library related to the solver is in the library path of the operating system

+------------------+---------------------------+----------------+
| Operating system | Path variable             | Extension      |
+==================+===========================+================+
| Linux            | :code:`LD_LIBRARY_PATH`   | :code:`.so`    |
+------------------+---------------------------+----------------+
| Mac              | :code:`DYLD_LIBRARY_PATH` | :code:`.dylib` |
+------------------+---------------------------+----------------+
| Windows          | :code:`PATH`              | :code:`.dll`   |
+------------------+---------------------------+----------------+





QDLDL
---------------
OSQP comes with `QDLDL <https://github.com/oxfordcontrol/qdldl>`_ internally installed.
It does not require any external shared library.
QDLDL is a sparse direct solver that works well for most small to medium sized problems.
However, its becomes not really efficient for large scale problems since it is not multi-threaded.


MKL Pardiso
-----------
`MKL Pardiso <https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface>`_ is an efficient multi-threaded linear system solver that works well for large scale problems part of the Intel Math Kernel Library.
Intel offers `free lincenses <https://software.intel.com/en-us/articles/free-mkl>`_ for MKL for most non-commercial applications.

Install with MKL
^^^^^^^^^^^^^^^^
We can install MKL Pardiso by using the standard `MKL installer <https://software.intel.com/en-us/mkl>`_.
The main library to be loaded is called :code:`libmkl_rt`.
To add it, together with its dependencies, to your path, just execute the automatic MKL script.

+------------------+------------------------------------------------+
| Operating system | Script                                         |
+==================+================================================+
| Linux            | :code:`source $MKLROOT/bin/mklvars.sh intel64` |
+------------------+------------------------------------------------+
| Mac              | :code:`source $MKLROOT/bin/mklvars.sh intel64` |
+------------------+------------------------------------------------+
| Windows          | :code:`%MKLROOT%\mklvars.bat intel64`          |
+------------------+------------------------------------------------+

where :code:`MKLROOT` is the MKL installation directory.

Install with Anaconda
^^^^^^^^^^^^^^^^^^^^^
`Anaconda Python distribution <https://www.anaconda.com/download/>`_ comes with the intel MKL libraries preinstalled including MKL Pardiso.
To use this version, the Anaconda libraries folders have to be added to the system path.
In particular, given the Anaconda installation directory :code:`ANACONDA_ROOT`, we need to add the following folders:

* :code:`libmkl_rt`: This library is located in the folder :code:`ANACONDA_ROOT/pkgs/mkl-VERSION/lib` where :code:`VERSION` is the version of the :code:`mkl` package.

* :code:`libiomp`: This library is located in the folder :code:`ANACONDA_ROOT/pkgs/intel-openmp-VERSION/lib` where :code:`VERSION` is the version of the :code:`intel-openmp` package.









