Python
======

Python interface supports Python 2.7 and 3.5 or newer.

Anaconda
--------

.. code:: bash

   conda install -c oxfordcontrol osqp


Pip
----

.. code:: bash

   pip install osqp


Sources
---------
You need to install the following (see :ref:`build_from_sources` for more details):

- `GCC compiler <https://gcc.gnu.org/>`_
- `CMake <https://cmake.org/>`_

.. note::

   **Windows**: You need to install **also** the Visual Studio C++ compiler:

   * Python 2: `Visual C++ 9.0 for Python (VC 9.0) <https://www.microsoft.com/en-us/download/details.aspx?id=44266>`_

   * Python 3: `Build Tools for Visual Studio 2017 <https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017>`_


Now you are ready to build OSQP python interface from sources. Run the following in your terminal

.. code:: bash

   git clone --recurse-submodules https://github.com/oxfordcontrol/osqp-python
   cd osqp-python
   python setup.py install
