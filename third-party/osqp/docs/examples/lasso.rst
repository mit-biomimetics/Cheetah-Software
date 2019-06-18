Lasso
=====


Lasso is a well known technique for sparse linear regression.
It is obtained by adding an :math:`\ell_1` regularization term in the objective,

.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} \| Ax - b \|_2^2 + \gamma \| x \|_1
  \end{array}


where :math:`x \in \mathbf{R}^{n}` is the vector of parameters, :math:`A \in \mathbf{R}^{m \times n}` is the data matrix, and :math:`\gamma > 0` is the weighting parameter.
The problem has the following equivalent form,

.. math::
  \begin{array}{ll}
    \mbox{minimize}   & \frac{1}{2} y^T y + \gamma \boldsymbol{1}^T t \\
    \mbox{subject to} & y = Ax - b \\
                      & -t \le x \le t
  \end{array}


In order to get a good trade-off between sparsity of the solution and quality of the linear fit, we solve the problem for varying weighting parameter :math:`\gamma`.
Since :math:`\gamma` enters only in the linear part of the objective function, we can reuse the matrix factorization and enable warm starting to reduce the computation time.



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
    Ad = sparse.random(m, n, density=0.5)
    x_true = np.multiply((np.random.rand(n) > 0.8).astype(float),
                         np.random.randn(n)) / np.sqrt(n)
    b = Ad.dot(x_true) + 0.5*np.random.randn(m)
    gammas = np.linspace(1, 10, 11)

    # Auxiliary data
    In = sparse.eye(n)
    Im = sparse.eye(m)
    On = sparse.csc_matrix((n, n))
    Onm = sparse.csc_matrix((n, m))

    # OSQP data
    P = sparse.block_diag((On, sparse.eye(m), On), format='csc')
    q = np.zeros(2*n + m)
    A = sparse.vstack([sparse.hstack([Ad, -Im, Onm.T]),
                       sparse.hstack([In, Onm, -In]),
                       sparse.hstack([In, Onm, In])]).tocsc()
    l = np.hstack([b, -np.inf * np.ones(n), np.zeros(n)])
    u = np.hstack([b, np.zeros(n), np.inf * np.ones(n)])

    # Create an OSQP object
    prob = osqp.OSQP()

    # Setup workspace
    prob.setup(P, q, A, l, u, warm_start=True)

    # Solve problem for different values of gamma parameter
    for gamma in gammas:
        # Update linear cost
        q_new = np.hstack([np.zeros(n+m), gamma*np.ones(n)])
        prob.update(q=q_new)

        # Solve
        res = prob.solve()


Matlab
------

.. code:: matlab

    % Generate problem data
    rng(1)
    n = 10;
    m = 1000;
    Ad = sprandn(m, n, 0.5);
    x_true = (randn(n, 1) > 0.8) .* randn(n, 1) / sqrt(n);
    b = Ad * x_true + 0.5 * randn(m, 1);
    gammas = linspace(1, 10, 11);

    % OSQP data
    P = blkdiag(sparse(n, n), speye(m), sparse(n, n));
    q = zeros(2*n+m, 1);
    A = [Ad, -speye(m), sparse(m,n);
        speye(n), sparse(n, m), -speye(n);
        speye(n), sparse(n, m), speye(n);];
    l = [b; -inf*ones(n, 1); zeros(n, 1)];
    u = [b; zeros(n, 1); inf*ones(n, 1)];

    % Create an OSQP object
    prob = osqp;

    % Setup workspace
    prob.setup(P, q, A, l, u, 'warm_start', true);

    % Solve problem for different values of gamma parameter
    for i = 1 : length(gammas)
        % Update linear cost
        gamma = gammas(i);
        q_new = [zeros(n+m,1); gamma*ones(n,1)];
        prob.update('q', q_new);

        % Solve
        res = prob.solve();
    end



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
    A = sparse.random(m, n, density=0.5)
    x_true = np.multiply((np.random.rand(n) > 0.8).astype(float),
                         np.random.randn(n)) / np.sqrt(n)
    b = A.dot(x_true) + 0.5*np.random.randn(m)
    gammas = np.linspace(1, 10, 11)

    # Define problem
    x = Variable(n)
    gamma = Parameter(nonneg=True)
    objective = 0.5*sum_squares(A*x - b) + gamma*norm1(x)
    prob = Problem(Minimize(objective))

    # Solve problem for different values of gamma parameter
    for gamma_val in gammas:
        gamma.value = gamma_val
        prob.solve(solver=OSQP, warm_start=True)


YALMIP
------

.. code:: matlab

    % Generate problem data
    rng(1)
    n = 10;
    m = 1000;
    A = sprandn(m, n, 0.5);
    x_true = (randn(n, 1) > 0.8) .* randn(n, 1) / sqrt(n);
    b = A * x_true + 0.5 * randn(m, 1);
    gammas = linspace(1, 10, 11);

    % Define problem
    x = sdpvar(n, 1);
    gamma = sdpvar;
    objective = 0.5*norm(A*x - b)^2 + gamma*norm(x,1);

    % Solve with OSQP
    options = sdpsettings('solver', 'osqp');
    x_opt = optimizer([], objective, options, gamma, x);

    % Solve problem for different values of gamma parameter
    for i = 1 : length(gammas)
        x_opt(gammas(i));
    end
