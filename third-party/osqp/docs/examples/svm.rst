Support vector machine (SVM)
============================

*Support vector machine* seeks an affine function that approximately classifies the two sets of points.
The problem can be stated as

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T x + \gamma \sum_{i=1}^{m} \max(0, b_i a_i^T x + 1),
  \end{array}

where :math:`b_i \in \{ -1, +1 \}` is a set label, and :math:`a_i` is a vector of features for the :math:`i`-th point.
The problem has the following equivalent form

.. math::
  \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} x^T x + \gamma \boldsymbol{1}^T t \\
    \mbox{subject to} & t \ge {\rm diag}(b) Ax + 1 \\
                      & t \ge 0,
  \end{array}

where :math:`{\rm diag}(b)` denotes the diagonal matrix with elements of :math:`b` on its diagonal.



Python
------

.. code:: python

    import osqp
    import numpy as np
    import scipy as sp
    import scipy.sparse as sparse

    # Generate problem data
    sp.random.seed(1)
    n = 10
    m = 1000
    N = int(m / 2)
    gamma = 1.0
    b = np.hstack([np.ones(N), -np.ones(N)])
    A_upp = sparse.random(N, n, density=0.5)
    A_low = sparse.random(N, n, density=0.5)
    Ad = sparse.vstack([
            A_upp / np.sqrt(n) + (A_upp != 0.).astype(float) / n,
            A_low / np.sqrt(n) - (A_low != 0.).astype(float) / n
         ]).tocsc()

    # OSQP data
    Im = sparse.eye(m)
    P = sparse.block_diag((sparse.eye(n), sparse.csc_matrix((m, m))), format='csc')
    q = np.hstack([np.zeros(n), gamma*np.ones(m)])
    A = sparse.vstack([
            sparse.hstack([sparse.diags(b).dot(Ad), -Im]),
            sparse.hstack([sparse.csc_matrix((m, n)), Im])
        ]).tocsc()
    l = np.hstack([-np.inf*np.ones(m), np.zeros(m)])
    u = np.hstack([-np.ones(m), np.inf*np.ones(m)])

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
    n = 10;
    m = 1000;
    N = ceil(m/2);
    gamma = 1;
    A_upp = sprandn(N, n, 0.5);
    A_low = sprandn(N, n, 0.5);
    Ad = [A_upp / sqrt(n) + (A_upp ~= 0) / n;
          A_low / sqrt(n) - (A_low ~= 0) / n];
    b = [ones(N, 1); -ones(N,1)];

    % OSQP data
    P = blkdiag(speye(n), sparse(m, m));
    q = [zeros(n,1); gamma*ones(m,1)];
    A = [diag(b)*Ad, -speye(m);
         sparse(m, n), speye(m)];
    l = [-inf*ones(m, 1); zeros(m, 1)];
    u = [-ones(m, 1); inf*ones(m, 1)];

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
    n = 10
    m = 1000
    N = int(m / 2)
    gamma = 1.0
    b = np.hstack([np.ones(N), -np.ones(N)])
    A_upp = sparse.random(N, n, density=0.5)
    A_low = sparse.random(N, n, density=0.5)
    A = sparse.vstack([
            A_upp / np.sqrt(n) + (A_upp != 0.).astype(float) / n,
            A_low / np.sqrt(n) - (A_low != 0.).astype(float) / n
        ]).tocsc()

    # Define problem
    x = Variable(n)
    objective = 0.5*sum_squares(x) + gamma*sum(pos(diag(b)*A*x + 1))

    # Solve with OSQP
    Problem(Minimize(objective)).solve(solver=OSQP)
    



YALMIP
------

.. code:: matlab

    % Generate problem data
    rng(1)
    n = 10;
    m = 1000;
    N = ceil(m/2);
    gamma = 1;
    A_upp = sprandn(N, n, 0.5);
    A_low = sprandn(N, n, 0.5);
    A = [A_upp / sqrt(n) + (A_upp ~= 0) / n;
         A_low / sqrt(n) - (A_low ~= 0) / n];
    b = [ones(N, 1); -ones(N,1)];

    % Define problem
    x = sdpvar(n, 1);
    objective = 0.5*norm(x)^2 + gamma*sum(max(diag(b)*A*x + 1, 0));

    % Solve with OSQP
    options = sdpsettings('solver','osqp');
    optimize([],objective,options);

