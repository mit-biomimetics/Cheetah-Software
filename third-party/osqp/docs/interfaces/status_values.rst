.. _status_values :

Status values
==============

These are the exit statuses, their respective constants and values returned by the solver as defined in `constants.h <https://github.com/oxfordcontrol/osqp/blob/master/include/constants.h>`_.
The *inaccurate* statuses define when the optimality, primal infeasibility or dual infeasibility conditions are satisfied with tolerances 10 times larger than the ones set.

+------------------------------+-----------------------------------+-------+
| Status                       | Constant                          | Value |
+==============================+===================================+=======+
| solved                       | OSQP_SOLVED                       | 1     |
+------------------------------+-----------------------------------+-------+
| solved inaccurate            | OSQP_SOLVED_INACCURATE            | 2     |
+------------------------------+-----------------------------------+-------+
| maximum iterations reached   | OSQP_MAX_ITER_REACHED             | -2    |
+------------------------------+-----------------------------------+-------+
| primal infeasible            | OSQP_PRIMAL_INFEASIBLE            | -3    |
+------------------------------+-----------------------------------+-------+
| primal infeasible inaccurate | OSQP_PRIMAL_INFEASIBLE_INACCURATE | 3     |
+------------------------------+-----------------------------------+-------+
| dual infeasible              | OSQP_DUAL_INFEASIBLE              | -4    |
+------------------------------+-----------------------------------+-------+
| dual infeasible inaccurate   | OSQP_DUAL_INFEASIBLE_INACCURATE   | 4     |
+------------------------------+-----------------------------------+-------+
| interrupted by user          | OSQP_SIGINT                       | -5    |
+------------------------------+-----------------------------------+-------+
| run time limit reached       | OSQP_TIME_LIMIT_REACHED           | -6    |
+------------------------------+-----------------------------------+-------+
| unsolved                     | OSQP_UNSOLVED                     | -10   |
+------------------------------+-----------------------------------+-------+
| problem non convex           | OSQP_NON_CVX                      | -7    |
+------------------------------+-----------------------------------+-------+

.. note::

   We recommend the user to **check the convexity of their problem before
   passing it to OSQP**! If the user passes a non-convex problem we do not
   assure the solver will be able to detect it.

   OSQP will try to detect **non-convex** problems by checking if the residuals
   diverge or if there are any issues in the initial factorization (if a direct
   method is used). It will detect non-convex problems when one or more of the
   eigenvalues of :code:`P` are "clearly" negative, i.e., when :code:`P + sigma
   * I` is not positive semidefinite. However, it might fail to detect
   non-convexity when :code:`P` has slightly negative eigenvalues, i.e., when
   :code:`P + sigma * I` is positive semidefinite and :code:`P` is not.


