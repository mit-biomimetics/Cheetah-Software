#include "glob_opts.h"

#include "qdldl.h"
#include "qdldl_interface.h"

#ifndef EMBEDDED
#include "amd.h"
#endif

#include "lin_alg.h"

#if EMBEDDED != 1
#include "kkt.h"
#endif

#ifndef EMBEDDED

// Free LDL Factorization structure
void free_linsys_solver_qdldl(qdldl_solver *s) {
    if (s) {
        if (s->L)         csc_spfree(s->L);
        if (s->P)         c_free(s->P);
        if (s->Dinv)      c_free(s->Dinv);
        if (s->bp)        c_free(s->bp);

        // These are required for matrix updates
        if (s->Pdiag_idx) c_free(s->Pdiag_idx);
        if (s->KKT)       csc_spfree(s->KKT);
        if (s->PtoKKT)    c_free(s->PtoKKT);
        if (s->AtoKKT)    c_free(s->AtoKKT);
        if (s->rhotoKKT)  c_free(s->rhotoKKT);
        if (s->D)         c_free(s->D);
        if (s->etree)     c_free(s->etree);
        if (s->Lnz)       c_free(s->Lnz);
        if (s->iwork)     c_free(s->iwork);
        if (s->bwork)     c_free(s->bwork);
        if (s->fwork)     c_free(s->fwork);
        c_free(s);

    }
}


/**
 * Compute LDL factorization of matrix A
 * @param  A    Matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus (0 is good)
 */
static c_int LDL_factor(csc *A,  qdldl_solver * p, c_int nvar){

    c_int sum_Lnz;
    c_int factor_status;

    // Compute elimination tree
    sum_Lnz = QDLDL_etree(A->n, A->p, A->i, p->iwork, p->Lnz, p->etree);

    if (sum_Lnz < 0){
      // Error
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the elimination tree. A is not perfectly upper triangular");
#endif
      return sum_Lnz;
    }

    // Allocate memory for Li and Lx
    p->L->i = (c_int *)c_malloc(sizeof(c_int)*sum_Lnz);
    p->L->x = (c_float *)c_malloc(sizeof(c_float)*sum_Lnz);

    // Factor matrix
    factor_status = QDLDL_factor(A->n, A->p, A->i, A->x,
                                 p->L->p, p->L->i, p->L->x,
                                 p->D, p->Dinv, p->Lnz,
                                 p->etree, p->bwork, p->iwork, p->fwork);


    if (factor_status < 0){
      // Error
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the nonzero elements. There are zeros in the diagonal matrix");
#endif
      return factor_status;
    } else if (factor_status < nvar) {
      // Error: Number of positive elements of D should be equal to nvar
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the nonzero elements. The problem seems to be non-convex");
#endif
      return -2;
    }

    return 0;

}


static c_int permute_KKT(csc ** KKT, qdldl_solver * p, c_int Pnz, c_int Anz, c_int m, c_int * PtoKKT, c_int * AtoKKT, c_int * rhotoKKT){
    c_float *info;
    c_int amd_status;
    c_int * Pinv;
    csc *KKT_temp;
    c_int * KtoPKPt;
    c_int i; // Indexing

    info = (c_float *)c_malloc(AMD_INFO * sizeof(c_float));

    // Compute permutation metrix P using AMD
    #ifdef DLONG
    amd_status = amd_l_order((*KKT)->n, (*KKT)->p, (*KKT)->i, p->P, (c_float *)OSQP_NULL, info);
    #else
    amd_status = amd_order((*KKT)->n, (*KKT)->p, (*KKT)->i, p->P, (c_float *)OSQP_NULL, info);
    #endif
    if (amd_status < 0) return (amd_status);


    // Inverse of the permutation vector
    Pinv = csc_pinv(p->P, (*KKT)->n);

    // Permute KKT matrix
    if (!PtoKKT && !AtoKKT && !rhotoKKT){  // No vectors to be stored
        // Assign values of mapping
        KKT_temp = csc_symperm((*KKT), Pinv, OSQP_NULL, 1);
    }
    else {
        // Allocate vector of mappings from unpermuted to permuted
        KtoPKPt = c_malloc((*KKT)->p[(*KKT)->n] * sizeof(c_int));
        KKT_temp = csc_symperm((*KKT), Pinv, KtoPKPt, 1);

        // Update vectors PtoKKT, AtoKKT and rhotoKKT
        if (PtoKKT){
            for (i = 0; i < Pnz; i++){
                PtoKKT[i] = KtoPKPt[PtoKKT[i]];
            }
        }
        if (AtoKKT){
            for (i = 0; i < Anz; i++){
                AtoKKT[i] = KtoPKPt[AtoKKT[i]];
            }
        }
        if (rhotoKKT){
            for (i = 0; i < m; i++){
                rhotoKKT[i] = KtoPKPt[rhotoKKT[i]];
            }
        }

        // Cleanup vector of mapping
        c_free(KtoPKPt);
    }

    // Cleanup
    // Free previous KKT matrix and assign pointer to new one
    csc_spfree((*KKT));
    (*KKT) = KKT_temp;
    // Free Pinv
    c_free(Pinv);
    // Free Amd info
    c_free(info);

    return 0;
}


// Initialize LDL Factorization structure
qdldl_solver *init_linsys_solver_qdldl(const csc * P, const csc * A, c_float sigma, c_float * rho_vec, c_int polish){

    c_int i;                     // loop counter

    // Define Variables
    qdldl_solver * p;  // Initialize LDL solver
    c_int n_plus_m;              // Define n_plus_m dimension
    csc * KKT_temp;              // Temporary KKT pointer

    // Allocate private structure to store KKT factorization
    p = c_calloc(1, sizeof(qdldl_solver));

    // Size of KKT
    n_plus_m = P->m + A->m;

    // Sparse matrix L (lower triangular)
    // NB: We don not allocate L completely (CSC elements)
    //      L will be allocated during the factorization depending on the
    //      resulting number of elements.
    p->L = c_malloc(sizeof(csc));
    p->L->m = n_plus_m;
    p->L->n = n_plus_m;
    p->L->nz = -1;

    // Diagonal matrix stored as a vector D
    p->Dinv = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);
    p->D    = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);

    // Permutation vector P
    p->P    = (QDLDL_int *)c_malloc(sizeof(QDLDL_int) * n_plus_m);

    // Working vector
    p->bp   = c_malloc(sizeof(QDLDL_float) * n_plus_m);

    // Elimination tree workspace
    p->etree = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));
    p->Lnz   = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));

    // Preallocate L matrix (Lx and Li are sparsity dependent)
    p->L->p = (c_int *)c_malloc((n_plus_m+1) * sizeof(QDLDL_int));

    //Lx and Li are sparsity dependent, so set them to
    //null initially so we don't try to free them prematurely
    p->L->i = OSQP_NULL;
    p->L->x = OSQP_NULL;


    // Preallocate workspace
    p->iwork = (QDLDL_int *)c_malloc(sizeof(QDLDL_int)*(3*n_plus_m));
    p->bwork = (QDLDL_bool *)c_malloc(sizeof(QDLDL_bool)*n_plus_m);
    p->fwork = (QDLDL_float *)c_malloc(sizeof(QDLDL_float)*n_plus_m);

    // Form and permute KKT matrix
    if (polish){ // Called from polish()
        // Use p->bp for storing param2 = vec(delta)
        for (i = 0; i < A->m; i++){
            p->bp[i] = sigma;
        }

        KKT_temp = form_KKT(P, A, 0, sigma, p->bp, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL);

        // Permute matrix
        permute_KKT(&KKT_temp, p, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL);
    }
    else { // Called from ADMM algorithm

        // Allocate vectors of indices
        p->PtoKKT = c_malloc((P->p[P->n]) * sizeof(c_int));
        p->AtoKKT = c_malloc((A->p[A->n]) * sizeof(c_int));
        p->rhotoKKT = c_malloc((A->m) * sizeof(c_int));

        // Use p->bp for storing param2 = rho_inv_vec
        for (i = 0; i < A->m; i++){
            p->bp[i] = 1. / rho_vec[i];
        }

        KKT_temp = form_KKT(P, A, 0, sigma, p->bp,
                            p->PtoKKT, p->AtoKKT,
                            &(p->Pdiag_idx), &(p->Pdiag_n), p->rhotoKKT);

        // Permute matrix
        permute_KKT(&KKT_temp, p, P->p[P->n], A->p[A->n], A->m, p->PtoKKT, p->AtoKKT, p->rhotoKKT);
    }

    // Check if matrix has been created
    if (!KKT_temp){
        #ifdef PRINTING
            c_eprint("Error forming and permuting KKT matrix");
        #endif
        return OSQP_NULL;
    }

    // Factorize the KKT matrix
    if (LDL_factor(KKT_temp, p, P->n) < 0) {
        csc_spfree(KKT_temp);
        free_linsys_solver_qdldl(p);
        return OSQP_NULL;
    }

    if (polish){ // If KKT passed, assign it to KKT_temp
        // Polish, no need for KKT_temp
        csc_spfree(KKT_temp);
    }
    else { // If not embedded option 1 copy pointer to KKT_temp. Do not free it.
        p->KKT = KKT_temp;
    }

    // Link Functions
    p->solve = &solve_linsys_qdldl;

    #ifndef EMBEDDED
    p->free = &free_linsys_solver_qdldl;
    #endif

    #if EMBEDDED != 1
    p->update_matrices = &update_linsys_solver_matrices_qdldl;
    p->update_rho_vec = &update_linsys_solver_rho_vec_qdldl;
    #endif

    // Assign type
    p->type = QDLDL_SOLVER;
    //
    // Set number of threads to 1 (single threaded)
    p->nthreads = 1;

    return p;
}

#endif  // EMBEDDED


// Permute x = P*b using P
void permute_x( c_int n, c_float * x,	c_float * b, c_int * P) {
    c_int j;
    for (j = 0 ; j < n ; j++) x[j] = b[P[j]];
}

// Permute x = P'*b using P
void permutet_x( c_int n, c_float * x,	c_float * b, c_int * P) {
    c_int j;
    for (j = 0 ; j < n ; j++) x[P[j]] = b[j];
}


static void LDLSolve(c_float *x, c_float *b, csc *L, c_float *Dinv, c_int *P,
              c_float *bp) {
    /* solves PLDL'P' x = b for x */
    permute_x(L->n, bp, b, P);
    QDLDL_solve(L->n, L->p, L->i, L->x, Dinv, bp);
    permutet_x(L->n, x, bp, P);

}


c_int solve_linsys_qdldl(qdldl_solver * s, c_float * b, const OSQPSettings *settings) {
    /* returns solution to linear system */
    /* Ax = b with solution stored in b */
    LDLSolve(b, b, s->L, s->Dinv, s->P, s->bp);

    return 0;
}


#if EMBEDDED != 1
// Update private structure with new P and A
c_int update_linsys_solver_matrices_qdldl(qdldl_solver * s,
		const csc *P, const csc *A, const OSQPSettings *settings){
    c_int kk;

    // Update KKT matrix with new P
    update_KKT_P(s->KKT, P, s->PtoKKT, settings->sigma, s->Pdiag_idx, s->Pdiag_n);

    // Update KKT matrix with new A
    update_KKT_A(s->KKT, A, s->AtoKKT);

    return (QDLDL_factor(s->KKT->n, s->KKT->p, s->KKT->i, s->KKT->x,
        s->L->p, s->L->i, s->L->x, s->D, s->Dinv, s->Lnz,
        s->etree, s->bwork, s->iwork, s->fwork) < 0);

}



c_int update_linsys_solver_rho_vec_qdldl(qdldl_solver * s, const c_float * rho_vec, const c_int m){
    c_int kk, i;

    // Use s->bp for storing param2 = rho_inv_vec
    for (i = 0; i < m; i++){
        s->bp[i] = 1. / rho_vec[i];
    }

    // Update KKT matrix with new rho
    update_KKT_param2(s->KKT, s->bp, s->rhotoKKT, m);

    return (QDLDL_factor(s->KKT->n, s->KKT->p, s->KKT->i, s->KKT->x,
        s->L->p, s->L->i, s->L->x, s->D, s->Dinv, s->Lnz,
        s->etree, s->bwork, s->iwork, s->fwork) < 0);
}


#endif
