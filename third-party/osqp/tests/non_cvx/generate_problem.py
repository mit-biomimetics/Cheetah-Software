import numpy as np
import scipy.sparse as spa
import utils.codegen_utils as cu

P = spa.csc_matrix(np.array([[2., 5.], [5., 1.]]))
q = np.array([3., 4.])

A = spa.csc_matrix(np.array([[-1.0, 0.], [0., -1.], [-1., 3.],
                             [2., 5.], [3., 4]]))
l = -np.inf * np.ones(A.shape[0])
u = np.array([0., 0., -15., 100., 80.])

sols_data = {'sigma_new': 5}

# Generate problem data
cu.generate_problem_data(P, q, A, l, u, 'non_cvx', sols_data)
