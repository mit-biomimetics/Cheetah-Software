Code generation
===============

OSQP can generate tailored C code that compiles into a fast and reliable solver
for the given family of QP problems in which the problem data, but not its
dimensions, change between problem instances.


The generated code is:

.. glossary::

      Malloc-free
          It does not perform any dynamic memory allocation.

      Library-free
          It is not linked to any external library.

      Division-free
          There are no division required in the ADMM algorithm


We make a distinction between two cases depending on which of the data are to be
treated as parameters.

.. glossary::

    Vectors as parameters
        Vectors :math:`q`, :math:`l` and :math:`u` can change between problem instances.
        This corresponds to the compiler flag :code:`EMBEDDED=1`.
        :math:`\rho` adaptation is not enabled. 
        The generated code is division-free and has a very small footprint.

    Matrices as parameters
        Both vectors :math:`q`, :math:`l`, :math:`u` and values in matrices
        :math:`P` and :math:`A` can change between problem instances.
        This corresponds to the compiler flag :code:`EMBEDDED=2`.
        :math:`\rho` adaptation is enabled but the frequency does not rely on any timing. 
        We assume that the sparsity patterns of :math:`P` and :math:`A` are fixed.


.. toctree::
   :maxdepth: 1
   :glob:

   *
