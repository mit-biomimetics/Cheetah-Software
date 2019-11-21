.. _cutest_interface:

CUTEst
======

The command to solve a problem in SIF format contained in the file
probname.SIF is

.. code::

    runcutest -p osqp -D probname.SIF

See the man page for runcutest for more details or other options.

The default OSQP settings used in CUTEst appear in the `OSQP.SPC <https://github.com/ralna/CUTEst/blob/master/src/osqp/OSQP.SPC>`_ file. 
Optionally, new parameter values to overwrite the default values can be stored in a file :code:`OSQP.SPC` in the directory where the :code:`runcutest` command is executed.
The format of the file :code:`OSQP.SPC` is the parameter name starting in the first column followed by one or more spaces and then the parameter value. 
The parameter names are case sensitive. 
If the parameter value is :code:`true` or :code:`false`, then use :code:`1` for true and :code:`0` for :code:`false`.

For more details see the `README.osqp <https://github.com/ralna/CUTEst/blob/master/src/osqp/README.osqp>`_ file in the CUTEst repository.

