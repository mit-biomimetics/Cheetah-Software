.. _build_from_sources:


Build from sources
==================

Install GCC and CMake
----------------------

The main compilation directives are specified using

- `GCC compiler <https://gcc.gnu.org/>`_ to build the binaries
- `CMake <https://cmake.org/>`__ to create the Makefiles


Linux
^^^^^
Both :code:`gcc` and :code:`cmake` commands are already installed by default.

Mac OS
^^^^^^

Install Xcode and command line tools
""""""""""""""""""""""""""""""""""""

#. Install the latest release of `Xcode <https://developer.apple.com/download/>`_.

#. Install the command line tools by executing from the terminal

    .. code:: bash

        xcode-select --install

Install CMake via Homebrew
"""""""""""""""""""""""""""

#. Install `Homebrew <https://brew.sh/>`_ and update its packages to the latest version.

#. Install cmake by executing

    .. code:: bash

        brew install cmake


Windows
^^^^^^^

#. Install `TDM-GCC <http://tdm-gcc.tdragon.net/download>`_ 32bit or 64bit depending on your platform.

#. Install the latest binaries of `CMake <https://cmake.org/download/#latest>`__.


Build the binaries
------------------

Run the following commands from the terminal

#. You first need to get the sources from one of these two options:

    * Download the compressed file `from bintray.com <https://dl.bintray.com/bstellato/generic/OSQP/0.4.1/osqp-0.4.1.tar.gz>`_.

    * Clone the repository

        .. code:: bash

            git clone --recursive https://github.com/oxfordcontrol/osqp

#. Create :code:`build` directory and change directory

        .. code:: bash

            cd osqp
            mkdir build
            cd build

#. Create Makefiles

    - In Linux and Mac OS run

        .. code:: bash

            cmake -G "Unix Makefiles" ..

    - In Windows run

        .. code:: bash

            cmake -G "MinGW Makefiles" ..


#. Compile OSQP

    .. code:: bash

       cmake --build .


Thanks to CMake, it is possible to create projects for a wide variety of IDEs; see `here <https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html>`_ for more details. For example, to create a project for Visual Studio 14 2015, it is just necessary to run

.. code:: bash

   cmake -G "Visual Studio 14 2015" ..


The compilation will generate the demo :code:`osqp_demo` and the unittests :code:`osqp_tester` executables. In the case of :code:`Unix` or :code:`MinGW` :code:`Makefiles` option they are located in the :code:`build/out/` directory.  Run them to check that the compilation was correct.


Once the sources are built, the generated static :code:`build/out/libosqp.a` and shared :code:`build/out/libosqp.ext` libraries can be used to interface any C/C++ software to OSQP (see :ref:`install_osqp_libs` installation).

.. _install_the_binaries:

Install the binaries
--------------------



To install the generated libraries and headers to a system-wide location compatible with `GNU standards <http://www.gnu.org/prep/standards/html_node/Directory-Variables.html>`_ it is just necessary to run

.. code:: bash

   cmake --build . --target install

This code installs the libraries in :code:`libdir` and the headers into :code:`includedir/osqp`. For mode details see the defaults folders on the `GNU standards <http://www.gnu.org/prep/standards/html_node/Directory-Variables.html>`_ website.
To change the installation prefix, in the "Create Makefiles" step above, you need to specify the destination folder as :code:`cmake -DCMAKE_INSTALL_PREFIX:PATH=myfolder ..`.

.. note:: This step requires write permissions in the destination
	  folders. You might be able to gain access using the
	  :code:`sudo` command.

We provided also an uninstall routine to remove the copied files by running

.. code:: bash

   cmake --build . --target uninstall

Note that this corresponds to running :code:`make install` and :code:`make uninstall` on unix machines.
