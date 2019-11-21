Portfolio optimization
======================


Portfolio optimization seeks to allocate assets in a way that maximizes the risk adjusted return,


.. math::
  \begin{array}{ll}
    \mbox{maximize} & \mu^T x - \gamma \left( x^T \Sigma x \right) \\
    \mbox{subject to} & \boldsymbol{1}^T x = 1 \\
                      & x \ge 0
  \end{array}


where :math:`x \in \mathbf{R}^{n}` represents the portfolio, :math:`\mu \in \mathbf{R}^{n}` the vector of expected returns, :math:`\gamma > 0` the risk aversion parameter, and :math:`\Sigma \in \mathbf{S}^{n}_{+}` the risk model covariance matrix.
The risk model is usually assumed to be the sum of a diagonal and a rank :math:`k < n` matrix,


.. math::
  \Sigma = F F^T + D,


where :math:`F \in \mathbf{R}^{n \times k}` is the factor loading matrix and :math:`D \in \mathbf{S}^{n}_{+}` is a diagonal matrix describing the asset-specific risk.
The resulting problem has the following equivalent form,

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T D x + \frac{1}{2} y^T y - \frac{1}{2\gamma}\mu^T x \\
    \mbox{subject to} & y = F^T x \\
                      & \boldsymbol{1}^T x = 1 \\
                      & x \ge 0
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
    n = 100
    k = 10
    F = sparse.random(n, k, density=0.7, format='csc')
    D = sparse.diags(np.random.rand(n) * np.sqrt(k), format='csc')
    mu = np.random.randn(n)
    gamma = 1

    # OSQP data
    P = sparse.block_diag((D, sparse.eye(k)), format='csc')
    q = np.hstack([-mu / (2*gamma), np.zeros(k)])
    A = sparse.vstack([
            sparse.hstack([F.T, -sparse.eye(k)]),
            sparse.hstack([sparse.csc_matrix(np.ones((1, n))), sparse.csc_matrix((1, k))]),
            sparse.hstack((sparse.eye(n), sparse.csc_matrix((n, k))))
        ]).tocsc()
    l = np.hstack([np.zeros(k), 1., np.zeros(n)])
    u = np.hstack([np.zeros(k), 1., np.ones(n)])

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
    n = 100;
    k = 10;
    F = sprandn(n, k, 0.7);
    D = sparse(diag( sqrt(k)*rand(n,1) ));
    mu = randn(n, 1);
    gamma = 1;

    % OSQP data
    P = blkdiag(D, speye(k));
    q = [-mu/(2*gamma); zeros(k, 1)];
    A = [F', -speye(k);
         ones(1, n), zeros(1, k);
         speye(n), sparse(n, k)];
    l = [zeros(k, 1); 1; zeros(n, 1)];
    u = [zeros(k, 1); 1; ones(n, 1)];

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
    n = 100
    k = 10
    F = sparse.random(n, k, density=0.7, format='csc')
    D = sparse.diags(np.random.rand(n) * np.sqrt(k), format='csc')
    mu = np.random.randn(n)
    gamma = 1
    Sigma = F*F.T + D

    # Define problem
    x = Variable(n)
    objective = mu.T*x - gamma*quad_form(x, Sigma)
    constraints = [sum(x) == 1, x >= 0]

    # Solve with OSQP
    Problem(Maximize(objective), constraints).solve(solver=OSQP)



YALMIP
------

.. code:: matlab

    % Generate problem data
    rng(1)
    n = 100;
    k = 10;
    F = sprandn(n, k, 0.7);
    D = sparse(diag( sqrt(k)*rand(n,1) ));
    mu = randn(n, 1);
    gamma = 1;
    Sigma = F*F' + D;

    % Define problem
    x = sdpvar(n, 1);
    objective = gamma * (x'*Sigma*x) - mu'*x;
    constraints = [sum(x) == 1, x >= 0];

    % Solve with OSQP
    options = sdpsettings('solver', 'osqp');
    optimize(constraints, objective, options);

