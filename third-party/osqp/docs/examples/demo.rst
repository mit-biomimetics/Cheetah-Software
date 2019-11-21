Demo
====


Consider the following QP


.. math::
  \begin{array}{ll}
    \mbox{minimize} & \frac{1}{2} x^T \begin{bmatrix}4 & 1\\ 1 & 2 \end{bmatrix} x + \begin{bmatrix}1 \\ 1\end{bmatrix}^T x \\
    \mbox{subject to} & \begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix} \leq \begin{bmatrix} 1 & 1\\ 1 & 0\\ 0 & 1\end{bmatrix} x \leq  \begin{bmatrix}1 \\ 0.7 \\ 0.7\end{bmatrix}
  \end{array}



We show below how to solve the problem in Python, Matlab, Julia and C.



Python
------

.. code:: python

    import osqp
    import scipy.sparse as sparse
    import numpy as np

    # Define problem data
    P = sparse.csc_matrix([[4, 1], [1, 2]])
    q = np.array([1, 1])
    A = sparse.csc_matrix([[1, 1], [1, 0], [0, 1]])
    l = np.array([1, 0, 0])
    u = np.array([1, 0.7, 0.7])

    # Create an OSQP object
    prob = osqp.OSQP()

    # Setup workspace and change alpha parameter
    prob.setup(P, q, A, l, u, alpha=1.0)

    # Solve problem
    res = prob.solve()



Matlab
------

.. code:: matlab

    % Define problem data
    P = sparse([4, 1; 1, 2]);
    q = [1; 1];
    A = sparse([1, 1; 1, 0; 0, 1]);
    l = [1; 0; 0];
    u = [1; 0.7; 0.7];

    % Create an OSQP object
    prob = osqp;

    % Setup workspace and change alpha parameter
    prob.setup(P, q, A, l, u, 'alpha', 1);

    % Solve problem
    res = prob.solve();



Julia
------

.. code:: julia

    import OSQP

    # Define problem data
    P = sparse([4. 1.; 1. 2.])
    q = [1.; 1.]
    A = sparse([1. 1.; 1. 0.; 0. 1.])
    u = [1.; 0.7; 0.7]
    l = [1.; 0.; 0.]

    # Crate OSQP object
    prob = OSQP.Model()

    # Setup workspace and change alpha parameter
    OSQP.setup!(prob; P=P, q=q, A=A, l=l, u=u, alpha=1)

    # Solve problem
    results = OSQP.solve!(prob)



C
-

.. code:: c

    #include "osqp.h"

    int main(int argc, char **argv) {
        // Load problem data
        c_float P_x[4] = {4.00, 1.00, 1.00, 2.00, };
        c_int P_nnz = 4;
        c_int P_i[4] = {0, 1, 0, 1, };
        c_int P_p[3] = {0, 2, 4, };
        c_float q[2] = {1.00, 1.00, };
        c_float A_x[4] = {1.00, 1.00, 1.00, 1.00, };
        c_int A_nnz = 4;
        c_int A_i[4] = {0, 1, 0, 2, };
        c_int A_p[3] = {0, 2, 4, };
        c_float l[3] = {1.00, 0.00, 0.00, };
        c_float u[3] = {1.00, 0.70, 0.70, };
        c_int n = 2;
        c_int m = 3;

        // Problem settings
        OSQPSettings * settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

        // Structures
        OSQPWorkspace * work;  // Workspace
        OSQPData * data;  // OSQPData

        // Populate data
        data = (OSQPData *)c_malloc(sizeof(OSQPData));
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;


        // Define Solver settings as default
        osqp_set_default_settings(settings);
        settings->alpha = 1.0; // Change alpha parameter

        // Setup workspace
        work = osqp_setup(data, settings);

        // Solve Problem
        osqp_solve(work);

        // Cleanup
        osqp_cleanup(work);
        c_free(data->A);
        c_free(data->P);
        c_free(data);
        c_free(settings);

        return 0;
    };
