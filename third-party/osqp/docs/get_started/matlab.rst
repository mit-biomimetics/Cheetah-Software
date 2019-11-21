Matlab
======
OSQP Matlab interface requires Matlab 2015b or newer.


Binaries
--------

Precompiled platform-dependent Matlab binaries are available on `Bintray <https://bintray.com/bstellato/generic/OSQP>`_.

To install the interface, just run the following commands:

.. code:: matlab

    websave('install_osqp.m', 'https://dl.bintray.com/bstellato/generic/OSQP/0.5.0/install_osqp.m');
    install_osqp


Sources
-------

You need to install the following (see :ref:`build_from_sources` for more details):

- A supported 64bit `C/C++ compiler <https://www.mathworks.com/support/compilers.html>`_
- `CMake <https://cmake.org/>`_



After you install both, check that your compiler is selected by executing

.. code:: matlab

   mex -setup

.. note::

   **Windows**: If Matlab does not find TDM-GCC, you need to set the environment variable :code:`MW_MINGW64_LOC` as follows

   .. code:: matlab

      setenv('MW_MINGW64_LOC', 'C:\TDM-GCC-64')


   where :code:`C:\TDM-GCC-64` is the installation folder for TDM-GCC.

You can now build the interface by running inside Matlab

.. code:: matlab

   !git clone --recurse-submodules https://github.com/oxfordcontrol/osqp-matlab
   cd osqp-matlab
   make_osqp


Then you can add the interface to the search path by executing from the same directory

.. code:: matlab

   addpath(pwd)
   savepath
