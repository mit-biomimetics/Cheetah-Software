Model predictive control (MPC)
==============================

We consider the problem of controlling a linear time-invariant dinamical system to some reference state :math:`x_r \in \mathbf{R}^{n_x}`.
To achieve this we use *constrained linear-quadratic MPC*, which solves at each time step the following finite-horizon optimal control problem

.. math::
  \begin{array}{ll}
    \mbox{minimize}   & (x_N-x_r)^T Q_N (x_N-x_r) + \sum_{k=0}^{N-1} (x_k-x_r)^T Q (x_k-x_r) + u_k^T R u_k \\
    \mbox{subject to} & x_{k+1} = A x_k + B u_k \\
                      & x_{\rm min} \le x_k  \le x_{\rm max} \\
                      & u_{\rm min} \le u_k  \le u_{\rm max} \\
                      & x_0 = \bar{x}
  \end{array}

The states :math:`x_k \in \mathbf{R}^{n_x}` and the inputs :math:`u_k \in \mathbf{R}^{n_u}` are constrained to be between some lower and upper bounds.
The problem is solved repeatedly for varying initial state :math:`\bar{x} \in \mathbf{R}^{n_x}`.


Python
------

.. code:: python

    import osqp
    import numpy as np
    import scipy as sp
    import scipy.sparse as sparse

    # Discrete time model of a quadcopter
    Ad = sparse.csc_matrix([
      [1.,      0.,     0., 0., 0., 0., 0.1,     0.,     0.,  0.,     0.,     0.    ],
      [0.,      1.,     0., 0., 0., 0., 0.,      0.1,    0.,  0.,     0.,     0.    ],
      [0.,      0.,     1., 0., 0., 0., 0.,      0.,     0.1, 0.,     0.,     0.    ],
      [0.0488,  0.,     0., 1., 0., 0., 0.0016,  0.,     0.,  0.0992, 0.,     0.    ],
      [0.,     -0.0488, 0., 0., 1., 0., 0.,     -0.0016, 0.,  0.,     0.0992, 0.    ],
      [0.,      0.,     0., 0., 0., 1., 0.,      0.,     0.,  0.,     0.,     0.0992],
      [0.,      0.,     0., 0., 0., 0., 1.,      0.,     0.,  0.,     0.,     0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      1.,     0.,  0.,     0.,     0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      0.,     1.,  0.,     0.,     0.    ],
      [0.9734,  0.,     0., 0., 0., 0., 0.0488,  0.,     0.,  0.9846, 0.,     0.    ],
      [0.,     -0.9734, 0., 0., 0., 0., 0.,     -0.0488, 0.,  0.,     0.9846, 0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      0.,     0.,  0.,     0.,     0.9846]
    ])
    Bd = sparse.csc_matrix([
      [0.,      -0.0726,  0.,     0.0726],
      [-0.0726,  0.,      0.0726, 0.    ],
      [-0.0152,  0.0152, -0.0152, 0.0152],
      [-0.,     -0.0006, -0.,     0.0006],
      [0.0006,   0.,     -0.0006, 0.0000],
      [0.0106,   0.0106,  0.0106, 0.0106],
      [0,       -1.4512,  0.,     1.4512],
      [-1.4512,  0.,      1.4512, 0.    ],
      [-0.3049,  0.3049, -0.3049, 0.3049],
      [-0.,     -0.0236,  0.,     0.0236],
      [0.0236,   0.,     -0.0236, 0.    ],
      [0.2107,   0.2107,  0.2107, 0.2107]])
    [nx, nu] = Bd.shape

    # Constraints
    u0 = 10.5916
    umin = np.array([9.6, 9.6, 9.6, 9.6]) - u0
    umax = np.array([13., 13., 13., 13.]) - u0
    xmin = np.array([-np.pi/6,-np.pi/6,-np.inf,-np.inf,-np.inf,-1.,
                     -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
    xmax = np.array([ np.pi/6, np.pi/6, np.inf, np.inf, np.inf, np.inf,
                      np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    # Objective function
    Q = sparse.diags([0., 0., 10., 10., 10., 10., 0., 0., 0., 5., 5., 5.])
    QN = Q
    R = 0.1*sparse.eye(4)

    # Initial and reference states
    x0 = np.zeros(12)
    xr = np.array([0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.])

    # Prediction horizon
    N = 10

    # Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
    # - quadratic objective
    P = sparse.block_diag([sparse.kron(sparse.eye(N), Q), QN,
                           sparse.kron(sparse.eye(N), R)]).tocsc()
    # - linear objective
    q = np.hstack([np.kron(np.ones(N), -Q.dot(xr)), -QN.dot(xr),
                   np.zeros(N*nu)])
    # - linear dynamics
    Ax = sparse.kron(sparse.eye(N+1),-sparse.eye(nx)) + sparse.kron(sparse.eye(N+1, k=-1), Ad)
    Bu = sparse.kron(sparse.vstack([sparse.csc_matrix((1, N)), sparse.eye(N)]), Bd)
    Aeq = sparse.hstack([Ax, Bu])
    leq = np.hstack([-x0, np.zeros(N*nx)])
    ueq = leq
    # - input and state constraints
    Aineq = sparse.eye((N+1)*nx + N*nu)
    lineq = np.hstack([np.kron(np.ones(N+1), xmin), np.kron(np.ones(N), umin)])
    uineq = np.hstack([np.kron(np.ones(N+1), xmax), np.kron(np.ones(N), umax)])
    # - OSQP constraints
    A = sparse.vstack([Aeq, Aineq]).tocsc()
    l = np.hstack([leq, lineq])
    u = np.hstack([ueq, uineq])

    # Create an OSQP object
    prob = osqp.OSQP()

    # Setup workspace
    prob.setup(P, q, A, l, u, warm_start=True)

    # Simulate in closed loop
    nsim = 15
    for i in range(nsim):
        # Solve
        res = prob.solve()

        # Check solver status
        if res.info.status != 'solved':
            raise ValueError('OSQP did not solve the problem!')

        # Apply first control input to the plant
        ctrl = res.x[-N*nu:-(N-1)*nu]
        x0 = Ad.dot(x0) + Bd.dot(ctrl)

        # Update initial state
        l[:nx] = -x0
        u[:nx] = -x0
        prob.update(l=l, u=u)



Matlab
------

.. code:: matlab

    % Discrete time model of a quadcopter
    Ad = [1       0       0   0   0   0   0.1     0       0    0       0       0;
          0       1       0   0   0   0   0       0.1     0    0       0       0;
          0       0       1   0   0   0   0       0       0.1  0       0       0;
          0.0488  0       0   1   0   0   0.0016  0       0    0.0992  0       0;
          0      -0.0488  0   0   1   0   0      -0.0016  0    0       0.0992  0;
          0       0       0   0   0   1   0       0       0    0       0       0.0992;
          0       0       0   0   0   0   1       0       0    0       0       0;
          0       0       0   0   0   0   0       1       0    0       0       0;
          0       0       0   0   0   0   0       0       1    0       0       0;
          0.9734  0       0   0   0   0   0.0488  0       0    0.9846  0       0;
          0      -0.9734  0   0   0   0   0      -0.0488  0    0       0.9846  0;
          0       0       0   0   0   0   0       0       0    0       0       0.9846];
    Bd = [0      -0.0726  0       0.0726;
         -0.0726  0       0.0726  0;
         -0.0152  0.0152 -0.0152  0.0152;
          0      -0.0006 -0.0000  0.0006;
          0.0006  0      -0.0006  0;
          0.0106  0.0106  0.0106  0.0106;
          0      -1.4512  0       1.4512;
         -1.4512  0       1.4512  0;
         -0.3049  0.3049 -0.3049  0.3049;
          0      -0.0236  0       0.0236;
          0.0236  0      -0.0236  0;
          0.2107  0.2107  0.2107  0.2107];
    [nx, nu] = size(Bd);

    % Constraints
    u0 = 10.5916;
    umin = [9.6; 9.6; 9.6; 9.6] - u0;
    umax = [13; 13; 13; 13] - u0;
    xmin = [-pi/6; -pi/6; -Inf; -Inf; -Inf; -1; -Inf(6,1)];
    xmax = [ pi/6;  pi/6;  Inf;  Inf;  Inf; Inf; Inf(6,1)];

    % Objective function
    Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);
    QN = Q;
    R = 0.1*eye(4);

    % Initial and reference states
    x0 = zeros(12,1);
    xr = [0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0];

    % Prediction horizon
    N = 10;

    % Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
    % - quadratic objective
    P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
    % - linear objective
    q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
    % - linear dynamics
    Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
    Bu = kron([sparse(1, N); speye(N)], Bd);
    Aeq = [Ax, Bu];
    leq = [-x0; zeros(N*nx, 1)];
    ueq = leq;
    % - input and state constraints
    Aineq = speye((N+1)*nx + N*nu);
    lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
    uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
    % - OSQP constraints
    A = [Aeq; Aineq];
    l = [leq; lineq];
    u = [ueq; uineq];

    % Create an OSQP object
    prob = osqp;

    % Setup workspace
    prob.setup(P, q, A, l, u, 'warm_start', true);

    % Simulate in closed loop
    nsim = 15;
    for i = 1 : nsim
        % Solve
        res = prob.solve();

        % Check solver status
        if ~strcmp(res.info.status, 'solved')
            error('OSQP did not solve the problem!')
        end

        % Apply first control input to the plant
        ctrl = res.x((N+1)*nx+1:(N+1)*nx+nu);
        x0 = Ad*x0 + Bd*ctrl;

        % Update initial state
        l(1:nx) = -x0;
        u(1:nx) = -x0;
        prob.update('l', l, 'u', u);
    end



CVXPY
-----

.. code:: python

    from cvxpy import *
    import numpy as np
    import scipy as sp
    import scipy.sparse as sparse

    # Discrete time model of a quadcopter
    Ad = sparse.csc_matrix([
      [1.,      0.,     0., 0., 0., 0., 0.1,     0.,     0.,  0.,     0.,     0.    ],
      [0.,      1.,     0., 0., 0., 0., 0.,      0.1,    0.,  0.,     0.,     0.    ],
      [0.,      0.,     1., 0., 0., 0., 0.,      0.,     0.1, 0.,     0.,     0.    ],
      [0.0488,  0.,     0., 1., 0., 0., 0.0016,  0.,     0.,  0.0992, 0.,     0.    ],
      [0.,     -0.0488, 0., 0., 1., 0., 0.,     -0.0016, 0.,  0.,     0.0992, 0.    ],
      [0.,      0.,     0., 0., 0., 1., 0.,      0.,     0.,  0.,     0.,     0.0992],
      [0.,      0.,     0., 0., 0., 0., 1.,      0.,     0.,  0.,     0.,     0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      1.,     0.,  0.,     0.,     0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      0.,     1.,  0.,     0.,     0.    ],
      [0.9734,  0.,     0., 0., 0., 0., 0.0488,  0.,     0.,  0.9846, 0.,     0.    ],
      [0.,     -0.9734, 0., 0., 0., 0., 0.,     -0.0488, 0.,  0.,     0.9846, 0.    ],
      [0.,      0.,     0., 0., 0., 0., 0.,      0.,     0.,  0.,     0.,     0.9846]
    ])
    Bd = sparse.csc_matrix([
      [0.,      -0.0726,  0.,     0.0726],
      [-0.0726,  0.,      0.0726, 0.    ],
      [-0.0152,  0.0152, -0.0152, 0.0152],
      [-0.,     -0.0006, -0.,     0.0006],
      [0.0006,   0.,     -0.0006, 0.0000],
      [0.0106,   0.0106,  0.0106, 0.0106],
      [0,       -1.4512,  0.,     1.4512],
      [-1.4512,  0.,      1.4512, 0.    ],
      [-0.3049,  0.3049, -0.3049, 0.3049],
      [-0.,     -0.0236,  0.,     0.0236],
      [0.0236,   0.,     -0.0236, 0.    ],
      [0.2107,   0.2107,  0.2107, 0.2107]])
    [nx, nu] = Bd.shape

    # Constraints
    u0 = 10.5916
    umin = np.array([9.6, 9.6, 9.6, 9.6]) - u0
    umax = np.array([13., 13., 13., 13.]) - u0
    xmin = np.array([-np.pi/6,-np.pi/6,-np.inf,-np.inf,-np.inf,-1.,
                     -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
    xmax = np.array([ np.pi/6, np.pi/6, np.inf, np.inf, np.inf, np.inf,
                      np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    # Objective function
    Q = sparse.diags([0., 0., 10., 10., 10., 10., 0., 0., 0., 5., 5., 5.])
    QN = Q
    R = 0.1*sparse.eye(4)

    # Initial and reference states
    x0 = np.zeros(12)
    xr = np.array([0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.])

    # Prediction horizon
    N = 10

    # Define problem
    u = Variable((nu, N))
    x = Variable((nx, N+1))
    x_init = Parameter(nx)
    objective = 0
    constraints = [x[:,0] == x_init]
    for k in range(N):
        objective += quad_form(x[:,k] - xr, Q) + quad_form(u[:,k], R)
        constraints += [x[:,k+1] == Ad*x[:,k] + Bd*u[:,k]]
        constraints += [xmin <= x[:,k], x[:,k] <= xmax]
        constraints += [umin <= u[:,k], u[:,k] <= umax]
    objective += quad_form(x[:,N] - xr, QN)
    prob = Problem(Minimize(objective), constraints)

    # Simulate in closed loop
    nsim = 15
    for i in range(nsim):
        x_init.value = x0
        prob.solve(solver=OSQP, warm_start=True)
        x0 = Ad.dot(x0) + Bd.dot(u[:,0].value)



YALMIP
------

.. code:: matlab

    % Discrete time model of a quadcopter
    Ad = [1       0       0   0   0   0   0.1     0       0    0       0       0;
          0       1       0   0   0   0   0       0.1     0    0       0       0;
          0       0       1   0   0   0   0       0       0.1  0       0       0;
          0.0488  0       0   1   0   0   0.0016  0       0    0.0992  0       0;
          0      -0.0488  0   0   1   0   0      -0.0016  0    0       0.0992  0;
          0       0       0   0   0   1   0       0       0    0       0       0.0992;
          0       0       0   0   0   0   1       0       0    0       0       0;
          0       0       0   0   0   0   0       1       0    0       0       0;
          0       0       0   0   0   0   0       0       1    0       0       0;
          0.9734  0       0   0   0   0   0.0488  0       0    0.9846  0       0;
          0      -0.9734  0   0   0   0   0      -0.0488  0    0       0.9846  0;
          0       0       0   0   0   0   0       0       0    0       0       0.9846];
    Bd = [0      -0.0726  0       0.0726;
         -0.0726  0       0.0726  0;
         -0.0152  0.0152 -0.0152  0.0152;
          0      -0.0006 -0.0000  0.0006;
          0.0006  0      -0.0006  0;
          0.0106  0.0106  0.0106  0.0106;
          0      -1.4512  0       1.4512;
         -1.4512  0       1.4512  0;
         -0.3049  0.3049 -0.3049  0.3049;
          0      -0.0236  0       0.0236;
          0.0236  0      -0.0236  0;
          0.2107  0.2107  0.2107  0.2107];
    [nx, nu] = size(Bd);

    % Constraints
    u0 = 10.5916;
    umin = [9.6; 9.6; 9.6; 9.6] - u0;
    umax = [13; 13; 13; 13] - u0;
    xmin = [-pi/6; -pi/6; -Inf; -Inf; -Inf; -1; -Inf(6,1)];
    xmax = [ pi/6;  pi/6;  Inf;  Inf;  Inf; Inf; Inf(6,1)];

    % Objective function
    Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);
    QN = Q;
    R = 0.1*eye(4);

    % Initial and reference states
    x0 = zeros(12,1);
    xr = [0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0];

    % Prediction horizon
    N = 10;

    % Define problem
    u = sdpvar(repmat(nu,1,N), repmat(1,1,N));
    x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));
    constraints = [xmin <= x{1} <= xmax];
    objective = 0;
    for k = 1 : N
        objective = objective + (x{k}-xr)'*Q*(x{k}-xr) + u{k}'*R*u{k};
        constraints = [constraints, x{k+1} == Ad*x{k} + Bd*u{k}];
        constraints = [constraints, umin <= u{k}<= umax, xmin <= x{k+1} <= xmax];
    end
    objective = objective + (x{N+1}-xr)'*QN*(x{N+1}-xr);
    options = sdpsettings('solver', 'osqp');
    controller = optimizer(constraints, objective, options, x{1}, [u{:}]);

    % Simulate in closed loop
    nsim = 15;
    for i = 1 : nsim
        U = controller{x0};
        x0 = Ad*x0 + Bd*U(:,1);
    end
