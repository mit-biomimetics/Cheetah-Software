module QDLDL

export qdldl, \, solve, solve!

using AMD

const QDLDL_UNKNOWN = -1;
const QDLDL_USED   = true;
const QDLDL_UNUSED = false;

struct QDLDLFactorisation{Tf<:AbstractFloat,Ti<:Integer}

    #contains factors L, D^{-1}
    L::SparseMatrixCSC{Tf,Ti}
    Dinv::Array{Tf}
    #permutation vector (nothing if no permutation)
    p
    #is it logical factorisation only?
    logical::Bool
end

#=
#constructor for enforcing no permutation
function qdldl{Tv<:AbstractFloat,Ti<:Integer}(A::SparseMatrixCSC{Tv,Ti},p::Void)
    #no permutation at all
    (L, Dinv) = factor(triu(A));
    return QDLDLFactorisation(L,Dinv,nothing)
end

#constructor for user-specified permutation, including no permutation
function qdldl{Tv<:AbstractFloat,Ti<:Integer}(A::SparseMatrixCSC{Tv,Ti},p::Array{Ti})
    (L, Dinv) = factor(triu(A[p,p]));
    return QDLDLFactorisation(L,Dinv,p)
end

#constructor for the default (AMD) permutation
function qdldl{Tv<:AbstractFloat,Ti<:Integer}(A::SparseMatrixCSC{Tv,Ti})
    #use AMD ordering as a default
    p = amd(A)::Array{Ti};
    return qdldl(A,p);
end
=#

# Usage :
# qdldl(A) uses the default AMD ordering
# qdldl(A,perm=p) uses a caller specified ordering
# qdldl(A,perm = nothing) factors without reordering
#
# qdldl(A,logical=true) produces a logical factorisation only

function qdldl{Tv<:AbstractFloat,Ti<:Integer}(
            A::SparseMatrixCSC{Tv,Ti};
            perm::Union{Array{Ti,1},Void}=amd(A),
            logical::Bool=false
         )
    Atr = perm == nothing ? triu(A) : triu(A[perm,perm])
    (L, Dinv) = factor(Atr,logical)
    return QDLDLFactorisation(L,Dinv,perm,logical)
end


function Base.:\(QDLDL::QDLDLFactorisation,b)
    return solve(QDLDL,b)
end

# Solves Ax = b using LDL factors for A.
# Returns x, preserving b
function solve(QDLDL::QDLDLFactorisation,b)
    x = copy(b)
    solve!(QDLDL,x)
    return x
end

# Solves Ax = b using LDL factors for A.
# Solves in place (x replaces b)
function solve!(QDLDL::QDLDLFactorisation,b)

    #bomb if logical factorisation only
    if QDLDL.logical
        error("Can't solve with logical factorisation only")
    end

    #permute b
    if QDLDL.p != nothing
        permute!(b,QDLDL.p)
    end
    QDLDL_solve!(QDLDL.L.n,
    QDLDL.L.colptr,
    QDLDL.L.rowval,
    QDLDL.L.nzval,
    QDLDL.Dinv,b)
    #inverse permutation
    if QDLDL.p != nothing
        ipermute!(b,QDLDL.p)
    end
    return nothing
end


function factor{Tv<:AbstractFloat,Ti<:Integer}(A::SparseMatrixCSC{Tv,Ti},logical)

    A = triu(A)

    etree  = Array{Ti}(A.n)
    Lnz    = Array{Ti}(A.n)
    iwork  = Array{Ti}(A.n*3)
    bwork  = Array{Bool}(A.n)
    fwork  = Array{Tv}(A.n)

    #compute elimination gree using QDLDL converted code
    sumLnz = QDLDL_etree!(A.n,A.colptr,A.rowval,iwork,Lnz,etree)

    if(sumLnz < 0)
        error("A matrix is not upper triangular or has an empty column")
    end

    #allocate space for the L matrix row indices and data
    Lp = Array{Ti}(A.n + 1)
    Li = Array{Ti}(sumLnz)
    Lx = Array{Tv}(sumLnz)
    #allocate for D and D inverse
    D  = Array{Tv}(A.n)
    Dinv = Array{Tv}(A.n)

    if(logical)
        Lx[:]   = 1
        D[:]    = 1
        Dinv[:] = 1
    end

    #factor using QDLDL converted code
    posDCount = QDLDL_factor!(A.n,A.colptr,A.rowval,A.nzval,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork,logical)

    if(posDCount < 0)
        error("Zero entry in D (matrix is not quasidefinite)")
    end

    L = SparseMatrixCSC(A.n,A.n,Lp,Li,Lx);

    return L, Dinv

end



# Compute the elimination tree for a quasidefinite matrix
# in compressed sparse column form.

function QDLDL_etree!(n,Ap,Ai,work,Lnz,etree)

    @inbounds for i = 1:n
        # zero out Lnz and work.  Set all etree values to unknown
        work[i]  = 0
        Lnz[i]   = 0i
        etree[i] = QDLDL_UNKNOWN

        #Abort if A doesn't have at least one entry
        #one entry in every column
        if(Ap[i] == Ap[i+1])
            return -1
        end
    end

    @inbounds for j = 1:n
        work[j] = j
        @inbounds for p = Ap[j]:(Ap[j+1]-1)
            i = Ai[p]
            if(i > j)
                return -1
            end
            @inbounds while(work[i] != j)
                if(etree[i] == QDLDL_UNKNOWN)
                    etree[i] = j
                end
                Lnz[i] += 1        #nonzeros in this column
                work[i] = j
                i = etree[i]
            end
        end #end for p
    end

    #tally the total nonzeros
    sumLnz = sum(Lnz)

    return sumLnz
end




function QDLDL_factor!(n,Ap,Ai,Ax,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork,logicalFactor)


    positiveValuesInD = 0

    #partition working memory into pieces
    yMarkers        = bwork
    yIdx            = view(iwork,      1:n)
    elimBuffer      = view(iwork,  (n+1):2*n)
    LNextSpaceInCol = view(iwork,(2*n+1):3*n)
    yVals           = fwork;


    Lp[1] = 1 #first column starts at index one / Julia is 1 indexed

    @inbounds for i = 1:n

        #compute L column indices
        Lp[i+1] = Lp[i] + Lnz[i]   #cumsum, total at the end

        # set all Yidx to be 'unused' initially
        #in each column of L, the next available space
        #to start is just the first space in the column
        yMarkers[i]  = QDLDL_UNUSED
        yVals[i]     = 0.0
        D[i]         = 0.0
        LNextSpaceInCol[i] = Lp[i]
    end

    if(!logicalFactor)
        # First element of the diagonal D.
        D[1]     = Ax[1]
        if(D[1] == 0.0) return -1 end
        if(D[1]  > 0.0) positiveValuesInD += 1 end
        Dinv[1] = 1/D[1];
    end

    #Start from 1 here. The upper LH corner is trivially 0
    #in L b/c we are only computing the subdiagonal elements
    @inbounds for k = 2:n

        #NB : For each k, we compute a solution to
        #y = L(0:(k-1),0:k-1))\b, where b is the kth
        #column of A that sits above the diagonal.
        #The solution y is then the kth row of L,
        #with an implied '1' at the diagonal entry.

        #number of nonzeros in this row of L
        nnzY = 0  #number of elements in this row

        #This loop determines where nonzeros
        #will go in the kth row of L, but doesn't
        #compute the actual values
        @inbounds for i = Ap[k]:(Ap[k+1]-1)

            bidx = Ai[i]   # we are working on this element of b

            #Initialize D[k] as the element of this column
            #corresponding to the diagonal place.  Don't use
            #this element as part of the elimination step
            #that computes the k^th row of L
            if(bidx == k)
                D[k] = Ax[i];
                continue
            end

            yVals[bidx] = Ax[i]   # initialise y(bidx) = b(bidx)

            # use the forward elimination tree to figure
            # out which elements must be eliminated after
            # this element of b
            nextIdx = bidx

            if(yMarkers[nextIdx] == QDLDL_UNUSED)  #this y term not already visited

                yMarkers[nextIdx] = QDLDL_USED     #I touched this one
                elimBuffer[1]     = nextIdx  # It goes at the start of the current list
                nnzE              = 1         #length of unvisited elimination path from here

                nextIdx = etree[bidx];

                @inbounds while(nextIdx != QDLDL_UNKNOWN && nextIdx < k)
                    if(yMarkers[nextIdx] == QDLDL_USED) break; end

                    yMarkers[nextIdx] = QDLDL_USED;   #I touched this one
                    #NB: Julia is 1-indexed, so I increment nnzE first here,
                    #no after writing into elimBuffer as in the C version
                    nnzE += 1                   #the list is one longer than before
                    elimBuffer[nnzE] = nextIdx; #It goes in the current list
                    nextIdx = etree[nextIdx];   #one step further along tree

                end #end while

                # now I put the buffered elimination list into
                # my current ordering in reverse order
                @inbounds while(nnzE != 0)
                    #NB: inc/dec reordered relative to C because
                    #the arrays are 1 indexed
                    nnzY += 1;
                    yIdx[nnzY] = elimBuffer[nnzE];
                    nnzE -= 1;
                end #end while
            end #end if

        end #end for i

        #This for loop places nonzeros values in the k^th row
        @inbounds for i = nnzY:-1:1

            #which column are we working on?
            cidx = yIdx[i]

            # loop along the elements in this
            # column of L and subtract to solve to y
            tmpIdx = LNextSpaceInCol[cidx];

            #don't compute Lx for logical factorisation 
            #this is not implemented in the C version
            if(!logicalFactor)
                yVals_cidx = yVals[cidx]
                @inbounds for j = Lp[cidx]:(tmpIdx-1)
                    yVals[Li[j]] -= Lx[j]*yVals_cidx
                end

                #Now I have the cidx^th element of y = L\b.
                #so compute the corresponding element of
                #this row of L and put it into the right place
                Lx[tmpIdx] = yVals_cidx *Dinv[cidx]

                #D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
                D[k] -= yVals_cidx*Lx[tmpIdx]
            end

            #also record which row it went into
            Li[tmpIdx] = k

            LNextSpaceInCol[cidx] += 1

            #reset the yvalues and indices back to zero and QDLDL_UNUSED
            #once I'm done with them
            yVals[cidx]     = 0.0
            yMarkers[cidx]  = QDLDL_UNUSED;

        end #end for i

        #Maintain a count of the positive entries
        #in D.  If we hit a zero, we can't factor
        #this matrix, so abort
        if(D[k] == 0.0) return -1 end
        if(D[k]  > 0.0) positiveValuesInD += 1 end

        #compute the inverse of the diagonal
        Dinv[k]= 1/D[k]

    end #end for k

    return positiveValuesInD

end

# Solves (L+I)x = b, with x replacing b
function QDLDL_Lsolve!(n,Lp,Li,Lx,x)

    @inbounds for i = 1:n
        @inbounds for j = Lp[i]: (Lp[i+1]-1)
            x[Li[j]] -= Lx[j]*x[i];
        end
    end
    return nothing
end


# Solves (L+I)'x = b, with x replacing b
function QDLDL_Ltsolve!(n,Lp,Li,Lx,x)

    @inbounds for i = n:-1:1
        @inbounds for j = Lp[i]:(Lp[i+1]-1)
            x[i] -= Lx[j]*x[Li[j]]
        end
    end
    return nothing
end

# Solves Ax = b where A has given LDL factors,
# with x replacing b
function QDLDL_solve!(n,Lp,Li,Lx,Dinv,b)

    QDLDL_Lsolve!(n,Lp,Li,Lx,b)
    b .*= Dinv;
    QDLDL_Ltsolve!(n,Lp,Li,Lx,b)

end



end #end module
