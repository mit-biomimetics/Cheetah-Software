#These types must be the same as in the QDLDL header
c_int   = Int64;
c_float = Float64;
c_bool  = Uint8;


function qdMatExample()

    i = [1, 7, 2, 3, 6, 10, 2, 3, 10, 4, 8, 5, 2, 6, 1, 7, 9, 4, 8, 7, 9, 2, 3, 10]
    j = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 10]
    s = [1.000000e+00, 3.503210e-01, 4.606408e-01, -1.211894e-01, -2.900576e-02, 1.824522e-01, -1.211894e-01, 4.179283e-01, -1.565056e+00, 1.778279e-01, -8.453948e-02, 1.000000e-01, -2.900576e-02, -1.000000e+00, 3.503210e-01, -4.410924e-01, 1.786632e-01, -8.453948e-02, -3.162278e-01, 1.786632e-01, -2.990768e-01, 1.824522e-01, -1.565056e+00, -1.000000e-01]

    A = sparse(i,j,s)
end


ldllibH = Libdl.dlopen("../../build/out/libqdldl." * Base.Libdl.dlext)
QDLDL_etree    = Libdl.dlsym(ldllibH, :QDLDL_etree)
QDLDL_factor   = Libdl.dlsym(ldllibH, :QDLDL_factor)
QDLDL_Lsolve   = Libdl.dlsym(ldllibH, :QDLDL_Lsolve)
QDLDL_Ltsolve  = Libdl.dlsym(ldllibH, :QDLDL_Ltsolve)
QDLDL_solve    = Libdl.dlsym(ldllibH, :QDLDL_solve)


A = qdMatExample()::SparseMatrixCSC{c_float,c_int}
A = sparse([1.0 1.0;1.0 -1.0]);
b = cumsum(ones(c_float(A.n)));

A = A::SparseMatrixCSC{c_float,c_int}

println("\n-------------------")
println("Testing: QDLDL_etree")
println("-------------------")
work   = zeros(c_int,A.n)
Lnz    = zeros(c_int,A.n)
etree  = zeros(c_int,A.n)

#Julia matrices are 1-index.  C code is 0 indexed.
#use only the upper triangle
Atriu = triu(A);
Ap = Atriu.colptr-1;
Ai = Atriu.rowval-1;
Ax = Atriu.nzval;

#Construct the elimination tree and column counts for A
sumLnz = ccall(QDLDL_etree,c_int,(c_int,Ptr{c_int},Ptr{c_int},Ptr{c_int},Ptr{c_int},Ptr{c_int}),A.n,Ap,Ai,work,Lnz,etree)


println("\n-------------------")
println("Testing: QDLDL_factor")
println("-------------------  ")

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
ccall(QDLDL_factor,c_int,
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

#Julia indexing is 1 based
L = SparseMatrixCSC(Ln,Ln,Lp+1,Li+1,Lx)
println("L = ", L)
println("D = ", D)
println("Inf-norm factorisation tolerance : ", norm( (L+I)*diagm(D)*(L'+I) - A, Inf)  );


println("\n-------------------")
println("Testing: QDLDL_solve")
println("-------------------  ")

#call the C solve function
x = copy(b)
ccall(QDLDL_solve,Void,
     (     c_int,              #n
           Ptr{c_int},         #Lp
           Ptr{c_int},         #Li
           Ptr{c_float},       #Lx
           Ptr{c_float},       #Dinv
           Ptr{c_float}),      #x (starts as b)
        Ln,Lp,Li,Lx,Dinv,x)

println("A\b solve tolerance: ", norm(x-A\b, Inf)  );

Libdl.dlclose(ldllibH)
