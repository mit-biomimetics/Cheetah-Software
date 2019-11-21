push!(LOAD_PATH, pwd())
using QDLDL

c_float = Float64
c_int   = Int64
c_bool  = Bool

function qdMatExample()

    i = [1, 7, 2, 3, 6, 10, 2, 3, 10, 4, 8, 5, 2, 6, 1, 7, 9, 4, 8, 7, 9, 2, 3, 10]
    j = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 10]
    s = [1.000000e+00, 3.503210e-01, 4.606408e-01, -1.211894e-01, -2.900576e-02, 1.824522e-01, -1.211894e-01, 4.179283e-01, -1.565056e+00, 1.778279e-01, -8.453948e-02, 1.000000e-01, -2.900576e-02, -1.000000e+00, 3.503210e-01, -4.410924e-01, 1.786632e-01, -8.453948e-02, -3.162278e-01, 1.786632e-01, -2.990768e-01, 1.824522e-01, -1.565056e+00, -1.000000e-01]

    A = sparse(i,j,s)
end

A = qdMatExample()
b = randn(c_float,A.n)
println("Example (backslash)                : ",  norm(qdldl(A)\b - A\b))

F = qdldl(A)
println("Example (with solve())             : ",  norm(solve(F,b) - A\b))

x = copy(b)
solve!(F,x)
println("Example (in place with solve!())   : ",  norm(x - A\b))

F = qdldl(A,perm=collect(A.n:-1:1))
println("Example (User permutation)         : ",  norm(F\b - A\b))

F = qdldl(A,perm=nothing)
println("Example (No permutation)           : ",  norm(F\b - A\b))


#compute a logical factorisation only
Flog = qdldl(A,logical = true)
Fnum = qdldl(A,logical = true)

logicalCheck = all((full(Fnum.L) .!= 0) .== (full(Flog.L) .!= 0))
println("Logical factorisation works?: ", logicalCheck)
