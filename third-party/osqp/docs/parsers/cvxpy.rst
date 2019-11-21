CVXPY
=====

`CVXPY <http://www.cvxpy.org/>`_ 1.0 supports the OSQP solver by default. After defining your problem, you can solve it with OSQP by just calling

.. code:: python

   problem.solve()


OSQP is the default QP solver used. To specify it explicitly together with some options you can execute

.. code:: python

   problem.solve(solver=OSQP, max_iter=2000)

where we set the :code:`max_iter` option to :code:`2000`.
