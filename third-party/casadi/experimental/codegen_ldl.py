from casadi import *

n = 5
colind = [0, 3, 6, 8, 10, 12]
row = [0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4]
nz = [19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18]
Asp = Sparsity(n, n, colind, row)
A0 = DM(Asp, nz)
A0 = mtimes(A0.T, A0)
Asp = A0.sparsity()
b0 = DM.ones(n)
x0 = solve(A0, b0)
print(x0)

# Test sparse QR
[D, LT, p] = ldl(A0)
x0 = ldl_solve(b0, D, LT, p)
print(x0)

A = MX.sym('A', Asp)
b = MX.sym('b', n)
x = solve(A, b, 'ldl')

for jit in [False, True]:
    # Create function
    f = Function('f', [A,b], [x], ['A', 'b'], ['x'], dict(jit=jit))

    # Numerical evaluation
    x0 = f(A0, b0)
    print(x0)
