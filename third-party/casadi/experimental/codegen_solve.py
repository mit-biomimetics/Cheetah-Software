from casadi import *

n = 5
colind = [0, 3, 6, 8, 10, 12]
row = [0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4]
nz = [19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18]
Asp = Sparsity(n, n, colind, row)
A0 = DM(Asp, nz)
b0 = DM.ones(n)
x0 = solve(A0, b0)
print(x0)

# Test sparse QR
[V, R, beta, prinv, pc] = qr_sparse(A0)
x0 = qr_solve(b0, V, R, beta, prinv, pc)
print(x0)

# Test sparse QR, transposed
[V, R, beta, prinv, pc] = qr_sparse(A0.T)
x0 = qr_solve(b0, V, R, beta, prinv, pc, True)
print(x0)

A = MX.sym('A', Asp)
b = MX.sym('b', n)
x = solve(A, b, 'qr')

for jit in [False, True]:
    # Create function
    f = Function('f', [A,b], [x], ['A', 'b'], ['x'], dict(jit=jit))

    # Numerical evaluation
    x0 = f(A0, b0)
    print(x0)
