import numpy as np
import scipy.sparse as spa
import utils.codegen_utils as cu

P = spa.diags([0.617022, 0.92032449, 0.20011437, 0.50233257, 0.34675589]).tocsc()
q = np.array([-1.10593508, -1.65451545, -2.3634686, 1.13534535, -1.01701414])
A = spa.csc_matrix((0,5))
l = np.array([])
u = np.array([])

# Generate problem solutions
sols_data = {'x_test': np.array([1.79237542, 1.79775228, 11.81058885, -2.26014678, 2.93293975]),
             'obj_value_test': -19.209752026813277,
             'status_test': 'optimal'}

# Generate problem data
cu.generate_problem_data(P, q, A, l, u, 'unconstrained', sols_data)