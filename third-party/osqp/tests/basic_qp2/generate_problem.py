import numpy as np
import scipy.sparse as spa
import utils.codegen_utils as cu

P = spa.csc_matrix(np.array([[11., 0.], [0., 0.]]))
q = np.array([3., 4.])

A = spa.csc_matrix(np.array([[-1.0, 0.], [0., -1.], [-1., 3.],
                             [2., 5.], [3., 4]]))
l = -np.inf * np.ones(A.shape[0])
u = np.array([0., 0., -15., 100., 80.])

n = P.shape[0]
m = A.shape[0]

# New data
q_new = np.array([1., 1.])
u_new = np.array([-2., 0., -20., 100., 80.])

# Generate problem solutions
sols_data = {'x_test': np.array([15., -0.]),
             'y_test': np.array([0., 508., 168., 0., 0.]),
             'obj_value_test': 1282.5,
             'status_test': 'optimal',
             'q_new': q_new,
             'u_new': u_new,
             'x_test_new': np.array([20., -0.]),
             'y_test_new': np.array([0., 664., 221., 0., 0.]),
             'obj_value_test_new': 2220.0,
             'status_test_new': 'optimal'}


# Generate problem data
cu.generate_problem_data(P, q, A, l, u, 'basic_qp2', sols_data)
