import numpy as np
import scipy.sparse as spa
import utils.codegen_utils as cu

# Set numpy seed for reproducibility
np.random.seed(2)


# Define tests
n = 5
m = 8
test_form_KKT_n = n
test_form_KKT_m = m
p = 0.7

test_form_KKT_A = spa.random(test_form_KKT_m, test_form_KKT_n, density=p).tocsc()
test_form_KKT_P = spa.random(n, n, density=p)
test_form_KKT_P = test_form_KKT_P.dot(test_form_KKT_P.T).tocsc() + spa.eye(n).tocsc()
test_form_KKT_Pu = spa.triu(test_form_KKT_P).tocsc()
test_form_KKT_rho = 1.6
test_form_KKT_sigma = 0.1
test_form_KKT_KKT = spa.vstack([
                        spa.hstack([test_form_KKT_P + test_form_KKT_sigma *
                        spa.eye(test_form_KKT_n), test_form_KKT_A.T]),
                        spa.hstack([test_form_KKT_A,
                        -1./test_form_KKT_rho * spa.eye(test_form_KKT_m)])]).tocsc()
test_form_KKT_KKTu = spa.triu(test_form_KKT_KKT).tocsc()


# Create new P, A and KKT
test_form_KKT_A_new = test_form_KKT_A.copy()
test_form_KKT_A_new.data += np.random.randn(test_form_KKT_A_new.nnz)
test_form_KKT_Pu_new = test_form_KKT_Pu.copy()
test_form_KKT_Pu_new.data += 0.1 * np.random.randn(test_form_KKT_Pu_new.nnz)
test_form_KKT_P_new = test_form_KKT_Pu_new + test_form_KKT_Pu_new.T - spa.diags(test_form_KKT_Pu_new.diagonal())

test_form_KKT_KKT_new = spa.vstack([
                        spa.hstack([test_form_KKT_P_new + test_form_KKT_sigma *
                        spa.eye(test_form_KKT_n), test_form_KKT_A_new.T]),
                        spa.hstack([test_form_KKT_A_new,
                        -1./test_form_KKT_rho * spa.eye(test_form_KKT_m)])]).tocsc()
test_form_KKT_KKTu_new = spa.triu(test_form_KKT_KKT_new).tocsc()


# Test solve problem with initial P and A
test_solve_P = test_form_KKT_P.copy()
test_solve_Pu = test_form_KKT_Pu.copy()
test_solve_q = np.random.randn(n)
test_solve_A = test_form_KKT_A.copy()
test_solve_l = -30 + np.random.randn(m)
test_solve_u = 30 + np.random.randn(m)


# Define new P
test_solve_P_new = test_form_KKT_P_new.copy()
test_solve_Pu_new = test_form_KKT_Pu_new.copy()


# Define new A
test_solve_A_new = test_form_KKT_A_new.copy()


# Generate test data and solutions
data = {'test_form_KKT_n': test_form_KKT_n,
        'test_form_KKT_m': test_form_KKT_m,
        'test_form_KKT_A': test_form_KKT_A,
        'test_form_KKT_P': test_form_KKT_P,
        'test_form_KKT_Pu': test_form_KKT_Pu,
        'test_form_KKT_rho': test_form_KKT_rho,
        'test_form_KKT_sigma': test_form_KKT_sigma,
        'test_form_KKT_KKT': test_form_KKT_KKT,
        'test_form_KKT_KKTu': test_form_KKT_KKTu,
        'test_form_KKT_A_new': test_form_KKT_A_new,
        'test_form_KKT_P_new': test_form_KKT_P_new,
        'test_form_KKT_Pu_new': test_form_KKT_Pu_new,
        'test_form_KKT_KKT_new': test_form_KKT_KKT_new,
        'test_form_KKT_KKTu_new': test_form_KKT_KKTu_new,
        'test_solve_P': test_solve_P,
        'test_solve_Pu': test_solve_Pu,
        'test_solve_q': test_solve_q,
        'test_solve_A': test_solve_A,
        'test_solve_l': test_solve_l,
        'test_solve_u': test_solve_u,
        'n': n,
        'm': m,
        'test_solve_x': np.array([-0.34967513, 1.20460722, -0.46259805,
                                  0.59083905, -0.87685541]),
        'test_solve_y': np.zeros(m),
        'test_solve_obj_value': -1.7665127080483103,
        'test_solve_status': 'optimal',
        'test_solve_P_new': test_solve_P_new,
        'test_solve_Pu_new': test_solve_Pu_new,
        'test_solve_P_new_x': np.array([-0.28228879, 1.3527703, -0.69277181,
                                        0.82445911, -1.11688134]),
        'test_solve_P_new_y': np.zeros(m),
        'test_solve_P_new_obj_value': -2.1490899311728526,
        'test_solve_P_new_status': 'optimal',
        'test_solve_A_new': test_solve_A_new,
        'test_solve_A_new_x': np.array([-0.34967513, 1.20460722, -0.46259805,
                                        0.59083905, -0.87685541]),
        'test_solve_A_new_y': np.zeros(m),
        'test_solve_A_new_obj_value': -1.7665127080484808,
        'test_solve_A_new_status': 'optimal',
        'test_solve_P_A_new_x': np.array([-0.28228879, 1.3527703, -0.69277181,
                                          0.82445911, -1.11688134]),
        'test_solve_P_A_new_y': np.zeros(m),
        'test_solve_P_A_new_obj_value': -2.1490899311726253,
        'test_solve_P_A_new_status': 'optimal'
        }


# Generate test data
cu.generate_data('update_matrices', data)
