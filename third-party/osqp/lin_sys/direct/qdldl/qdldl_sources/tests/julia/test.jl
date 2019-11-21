#These types must be the same as in the QDLDL header
c_int   = Int64;
c_float = Float64;
c_bool  = Uint8;

function LDLsolve(A::SparseMatrixCSC{c_float,c_int},b::Array{c_float})

    LDLlib = "../../build/out/libqdldl.dylib"

    work   = zeros(c_int,A.n)
    Lnz    = zeros(c_int,A.n)
    etree  = zeros(c_int,A.n)

    #Julia matrices are 1-index.  C code is 0 indexed.
    Ap = A.colptr-1;
    Ai = A.rowval-1;
    Ax = A.nzval;

    #Construct the elimination tree and column counts for A
    sumLnz = ccall((:QDLDL_etree,LDLlib),c_int,(c_int,Ptr{c_int},Ptr{c_int},Ptr{c_int},Ptr{c_int},Ptr{c_int}),A.n,Ap,Ai,work,Lnz,etree)

    if(sumLnz < 0)
        return "Abort on etree";
    end

    #allocate target matrix data for L
    Ln = A.n;
    Lp = zeros(c_int,A.n+1);
    Li = zeros(c_int,sumLnz);
    Lx = zeros(c_float,sumLnz);

    #allocate D and its inverse
    D      = zeros(c_float,A.n)
    Dinv   = zeros(c_float,A.n)

    #scratch workspace values
    fwork            = zeros(c_float,Ln)
    iwork            = zeros(c_int,Ln*3)
    bwork            = zeros(c_bool,Ln)

    #call the C factorization function
    posCountD = ccall((:QDLDL_factor,LDLlib),c_int,
         (     c_int,              #n
               Ptr{c_int},         #Ap
               Ptr{c_int},         #Ai
               Ptr{c_float},       #Ax
               Ptr{c_int},         #Lp
               Ptr{c_int},         #Li
               Ptr{c_float},       #Lx
               Ptr{c_float},       #D
               Ptr{c_float},       #Dinv
               Ptr{c_int},         #Lnz
               Ptr{c_int},         #etree
               Ptr{c_bool},        #bwork
               Ptr{c_int},         #iwork
               Ptr{c_float}),      #fwork
          A.n,Ap,Ai,Ax,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork)

    if(posCountD < 0)
        return "Abort on factor";
    end

    #Julia indexing is 1 based
    L = SparseMatrixCSC(Ln,Ln,Lp+1,Li+1,Lx)

    #call the C solve function
    x = copy(b)
    ccall((:QDLDL_solve,LDLlib),Void,
         (     c_int,              #n
               Ptr{c_int},         #Lp
               Ptr{c_int},         #Li
               Ptr{c_float},       #Lx
               Ptr{c_float},       #Dinv
               Ptr{c_float}),      #x (starts as b)
            Ln,Lp,Li,Lx,Dinv,x)

    return x
end


function qdMatExample()

    i = [1, 7, 2, 3, 6, 10, 2, 3, 10, 4, 8, 5, 2, 6, 1, 7, 9, 4, 8, 7, 9, 2, 3, 10]
    j = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 10]
    s = [1.000000e+00, 3.503210e-01, 4.606408e-01, -1.211894e-01, -2.900576e-02, 1.824522e-01, -1.211894e-01, 4.179283e-01, -1.565056e+00, 1.778279e-01, -8.453948e-02, 1.000000e-01, -2.900576e-02, -1.000000e+00, 3.503210e-01, -4.410924e-01, 1.786632e-01, -8.453948e-02, -3.162278e-01, 1.786632e-01, -2.990768e-01, 1.824522e-01, -1.565056e+00, -1.000000e-01]

    A = sparse(i,j,s)
end


function qdRandom(m,n)

    density = 0.1
    A =  pdRandom(m,density)
    C = -pdRandom(n,density)
    B = sprandn(m,n,density)
    H = [A B;B' C]

end


function pdRandom(n,density)

    X = sprandn(n,n,density)
    X = 0.5*(X+X') + n*speye(n)

end

println("Things that should fail")
println("----------------------------")
A = sparse([0.0 1.0; 1.0 0.0])
b = randn(c_float,A.n);
println("non LDL decomposable   : ",  LDLsolve(triu(A),b))

A = sparse([1.0 1.0; 1.0 1.0])
b = randn(c_float,A.n);
println("Tril structure         : ",  LDLsolve(tril(A),b))

A = sparse([5.0 1.0; 1.0 5.0])
b = randn(c_float,A.n);
println("TRIL and TRIL parts    : ",  LDLsolve(A,b))

A = sparse([0.0 5.0; 5.0 1.0])
b = randn(c_float,A.n);
println("Zero at A[0,0]         : ",  LDLsolve(triu(A),b))

A = sparse([1.0 1.0; 1.0 1.0])
b = randn(c_float,A.n);
println("Rank deficient matrix  : ",  LDLsolve(triu(A),b))

A = qdMatExample();
A[5,5] = 0.0;
b = randn(c_float,A.n);
println("Zero with allocation   : ",  LDLsolve(triu(A),b))


println("\n\n")
println("Things that should work   ")
println("(Errors are relative to A\\b)")
println("-----------------------------")

A = qdMatExample();
b = randn(c_float,A.n);
println("Basic example                : ",  norm(A\b-LDLsolve(triu(A),b)))

A = sparse([4.0 1.0 2.0; 1.0 0.0 1.0;2.0 1.0 -3.0])
b = [6.0,9.0,12.0];
println("Zero mid matrix              : ",  norm(A\b-LDLsolve(triu(A),b)))

A = sparse([1.0 5.0; 5.0 0.0])
b = randn(c_float,A.n);
println("Zero at A[end,end]           : ",  norm(A\b-LDLsolve(triu(A),b)))

A = sparse([1.0 1.0;1.0 -1.0]);
b = cumsum(ones(c_float(A.n)));
println("2x2                          : ", norm(A\b- LDLsolve(triu(A),b)))

A = speye(100);
b = ones(c_float(A.n));
println("Big identity matrix          : ", norm(A\b- LDLsolve(triu(A),b)))

A = qdMatExample();
b = randn(c_float,A.n);
println("example QD                   : ", norm(A\b- LDLsolve(triu(A),b)))

A = speye(1)/5.0;
b = randn(c_float,A.n);
println("singleton                    : ", norm(A\b- LDLsolve(triu(A),b)))

A = sparse([1.0 1.0;1.0 (1.0-1e-8)]);
b = randn(c_float,A.n);
println("ill conditioned              : ", norm(A\b- LDLsolve(triu(A),b)))

A = sparse([1.0 1.0;1.0 -(1.0 -1e-8)]);
b = randn(c_float,A.n);
println("ill conditioned QD           : ", norm(A\b- LDLsolve(triu(A),b)))

A = sparse([1e-8 1.0;1.0 -(1.0 -1e-8)]);
b = randn(c_float,A.n);
println("+Epsilon at A[0,0]           : ", norm(A\b- LDLsolve(triu(A),b)))

A = qdRandom(100,90);
b = randn(c_float,A.n);
println("Big QD                       : ", norm(A\b- LDLsolve(triu(A),b)))

A = qdRandom(0,90);
b = randn(c_float,A.n);
println("Big ND                       : ", norm(A\b- LDLsolve(triu(A),b)))

A = qdRandom(100,0);
b = randn(c_float,A.n);
println("Big PD                       : ", norm(A\b- LDLsolve(triu(A),b)))

A = sparse([1.0 5.0; 5.0 1.0])
b = randn(c_float,A.n);
println("indefinite matrix            : ", norm(A\b- LDLsolve(triu(A),b)))

Ai = [0; 1; 2; 1; 0; 3; 4; 5; 5; 6; 4; 3] + 1
Aj = [0; 1; 2; 2; 2; 3; 4; 5; 6; 6; 6; 6] + 1
Ax = [-0.25; -0.25; 1; 0.513578; 0.529142; -0.25; -0.25; 1.10274; 0.15538; 1.25883; 0.13458; 0.621134]
A = sparse(Ai, Aj, Ax)
b = randn(c_float,A.n);
println("tricky permuted osqp matrix  : ", norm((A + A' - diagm(diag(A)))\b- LDLsolve(A, b)))
