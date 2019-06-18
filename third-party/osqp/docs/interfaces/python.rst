.. _python_interface:

Python
======

Import
------
The OSQP module can be imported with

.. code:: python

    import osqp


.. _python_setup:

Setup
-----

The solver is initialized by creating an OSQP object

.. code:: python

    m = osqp.OSQP()

The problem is specified in the setup phase by running

.. code:: python

    m.setup(P=P, q=q, A=A, l=l, u=u, **settings)


The arguments :code:`q`, :code:`l` and :code:`u` are numpy arrays. The elements of :code:`l` and :code:`u` can be :math:`\pm \infty` ( using :code:`numpy.inf`).

The arguments :code:`P` and :code:`A` are scipy sparse matrices in CSC format. 
Matrix :code:`P` can be either complete or just the upper triangular
part. OSQP will make use of only the upper triangular part.
If they are sparse matrices are in another format, the interface will attempt to convert them. There is no need to specify all the arguments.


The keyword arguments :code:`**settings` specify the solver settings. The allowed parameters are defined in :ref:`solver_settings`.

Solve
-----

The problem can be solved by

.. code:: python

   results = m.solve()


The :code:`results` object contains the primal solution :code:`x`, the dual solution :code:`y`, certificate of primal infeasibility :code:`prim_inf_cert`, certificate of dual infeasibility :code:`dual_inf_cert` and the :code:`info` object containing the solver statistics defined in the following table


+-----------------------+------------------------------------------------+
| Member                | Description                                    |
+=======================+================================================+
| :code:`iter`          | Number of iterations                           |
+-----------------------+------------------------------------------------+
| :code:`status`        | Solver status                                  |
+-----------------------+------------------------------------------------+
| :code:`status_val`    | Solver status value as in :ref:`status_values` |
+-----------------------+------------------------------------------------+
| :code:`status_polish` | Polishing status                               |
+-----------------------+------------------------------------------------+
| :code:`obj_val`       | Objective value                                |
+-----------------------+------------------------------------------------+
| :code:`pri_res`       | Primal residual                                |
+-----------------------+------------------------------------------------+
| :code:`dua_res`       | Dual residual                                  |
+-----------------------+------------------------------------------------+
| :code:`setup_time`    | Setup time                                     |
+-----------------------+------------------------------------------------+
| :code:`solve_time`    | Solve time                                     |
+-----------------------+------------------------------------------------+
| :code:`update_time`   | Update time                                    |
+-----------------------+------------------------------------------------+
| :code:`polish_time`   | Polish time                                    |
+-----------------------+------------------------------------------------+
| :code:`run_time`      | Total run time: setup/update + solve + polish  |
+-----------------------+------------------------------------------------+
| :code:`rho_estimate`  | Optimal rho estimate                           |
+-----------------------+------------------------------------------------+
| :code:`rho_updates`   | Number of rho updates                          |
+-----------------------+------------------------------------------------+

Note that if multiple solves are executed from single setup, then after the
first one :code:`run_time` includes :code:`update_time` + :code:`solve_time`
+ :code:`polish_time`.


Update
------
Part of problem data and settings can be updated without requiring a new problem setup.

Update problem vectors
^^^^^^^^^^^^^^^^^^^^^^
Vectors :code:`q`, :code:`l` and :code:`u` can be updated with new values :code:`q_new`, :code:`l_new` and :code:`u_new` by just running

.. code:: python

    m.update(q=q_new, l=l_new, u=u_new)


The user does not have to specify all the keyword arguments.


.. _python_update_settings:

Update problem matrices
^^^^^^^^^^^^^^^^^^^^^^^^
Matrices :code:`A` and :code:`P` can be updated by changing the value of their elements but not their sparsity pattern. 
The interface is designed to mimic the :ref:`C/C++ counterpart <c_cpp_update_data>`. 
Note that the new values of :code:`P` represent only the upper triangular part while :code:`A` is always represented as a full matrix.

You can update the values of all the elements of :code:`P` by executing

.. code:: python

    m.update(Px=Px_new)


If you want to update only some elements, you can pass

.. code:: python

    m.update(Px=Px_new, Px_idx=Px_new_idx)

where :code:`Px_new_idx` is the vector of indices of mapping the elements of :code:`Px_new` to the original vector :code:`Px` representing the data of the sparse matrix :code:`P`.

Matrix :code:`A` can be changed in the same way. You can also change both matrices at the same time by running, for example


.. code:: python

    m.update(Px=Px_new, Px_idx=Px_new_idx, Ax=Ax_new, Ax=Ax_new_idx)


Update settings
^^^^^^^^^^^^^^^

Settings can be updated by running

.. code:: python

    m.update_settings(**kwargs)


where :code:`kwargs` are the settings that can be updated which are marked with an * in :ref:`solver_settings`.


Warm start
----------

OSQP automatically warm starts primal and dual variables from the previous QP solution. If you would like to warm start their values manually, you can use

.. code:: python

    m.warm_start(x=x0, y=y0)


where :code:`x0` and :code:`y0` are the new primal and dual variables. 
