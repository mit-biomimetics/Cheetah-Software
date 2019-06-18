import numpy as np
import scipy.sparse as spa
import utils.codegen_utils as cu

P = spa.diags([1., 0.]).tocsc()
q = np.array([1., -1.])

A12 = spa.csc_matrix([[1., 1.], [1., 0.], [0., 1.]])
A34 = spa.csc_matrix([[1., 0.], [1., 0.], [0., 1.]])
l = np.array([0., 1., 1.])
u1 = np.array([5., 3., 3.])
u2 = np.array([0., 3., 3.])
u3 = np.array([2., 3., np.inf])
u4 = np.array([0., 3., np.inf])

# Generate problem solutions
data = {'P': P,
        'q': q,
        'A12': A12,
        'A34': A34,
        'l': l,
        'u1': u1,
        'u2': u2,
        'u3': u3,
        'u4': u4,
        'x1': np.array([1., 3.]),
        'y1': np.array([0., -2., 1.]),
        'obj_value1': -1.5,
        'status1': 'optimal',
        'status2': 'primal_infeasible',
        'status3': 'dual_infeasible',
        'status4': 'primal_infeasible'
        }

# Generate problem data
cu.generate_data('primal_dual_infeasibility', data)
