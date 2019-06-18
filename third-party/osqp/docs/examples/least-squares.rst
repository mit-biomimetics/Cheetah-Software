Least-squares
=============

Consider the following constrained least-squares problem

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \|Ax - b\|_2^2 \\
    \mbox{subject to} & 0 \leq x \leq 1
  \end{array}

The problem has the following equivalent form

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} y^T y \\
    \mbox{subject to} & y = A x - b \\
                      & 0 \le x \le 1
  \end{array}



Python
------

.. code:: python

    import osqp
    import numpy as np
    import scipy as sp
    import scipy.sparse as sparse

    # Generate problem data
    sp.random.seed(1)
    m = 30
    n = 20
    Ad = sparse.random(m, n, density=0.7, format='csc')
    b = np.random.randn(m)

    # OSQP data
    P = sparse.block_diag((sparse.csc_matrix((n, n)), sparse.eye(m)), format='csc')
    q = np.zeros(n+m)
    A = sparse.vstack([
            sparse.hstack([Ad, -sparse.eye(m)]),
            sparse.hstack((sparse.eye(n), sparse.csc_matrix((n, m))))
        ]).tocsc()
    l = np.hstack([b, np.zeros(n)])
    u = np.hstack([b, np.ones(n)])

    # Create an OSQP object
    prob = osqp.OSQP()

    # Setup workspace
    prob.setup(P, q, A, l, u)

    # Solve problem
    res = prob.solve()



Matlab
------

.. code:: matlab

   % Generate problem data
   rng(1)
   m = 30;
   n = 20;
   Ad = sprandn(m, n, 0.7);
   b = randn(m, 1);

   % OSQP data
   P = blkdiag(sparse(n, n), speye(m));
   q = zeros(n+m, 1);
   A = [Ad, -speye(m);
        speye(n), sparse(n, m)];
   l = [b; zeros(n, 1)];
   u = [b; ones(n, 1)];

   % Create an OSQP object
   prob = osqp;

   % Setup workspace
   prob.setup(P, q, A, l, u);

   % Solve problem
   res = prob.solve();



CVXPY
-----

.. code:: python

    from cvxpy import *
    import numpy as np
    import scipy as sp
    import scipy.sparse as sparse

    # Generate problem data
    sp.random.seed(1)
    m = 30
    n = 20
    A = sparse.random(m, n, density=0.7, format='csc')
    b = np.random.randn(m)

    # Define problem
    x = Variable(n)
    objective = 0.5*sum_squares(A*x - b)
    constraints = [x >= 0, x <= 1]

    # Solve with OSQP
    Problem(Minimize(objective), constraints).solve(solver=OSQP)



YALMIP
------

.. code:: matlab

   % Generate data
   rng(1)
   m = 30;
   n = 20;
   A = sprandn(m, n, 0.7);
   b = randn(m, 1);

   % Define problem
   x = sdpvar(n, 1);
   objective = 0.5*norm(A*x - b)^2;
   constraints = [ 0 <= x <= 1];

   % Solve with OSQP
   options = sdpsettings('solver','osqp');
   optimize(constraints, objective, options);

