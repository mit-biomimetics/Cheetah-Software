#include "pardiso_interface.h"

#define MKL_INT c_int

// Single Dynamic library interface
#define MKL_INTERFACE_LP64  0x0
#define MKL_INTERFACE_ILP64 0x1

// Solver Phases
#define PARDISO_SYMBOLIC  (11)
#define PARDISO_NUMERIC   (22)
#define PARDISO_SOLVE     (33)
#define PARDISO_CLEANUP   (-1)


// Prototypes for Pardiso functions
void pardiso(void**, const c_int*, const c_int*, const c_int*, const c_int*,
             const c_int*, const c_float*, const c_int*, const c_int*,
             const c_int*, const c_int*, c_int*, const c_int*, c_float*,
             c_float*, const c_int*);
c_int mkl_set_interface_layer(c_int);
c_int mkl_get_max_threads();

// Free LDL Factorization structure
void free_linsys_solver_pardiso(pardiso_solver *s) {
    if (s) {

        // Free pardiso solver using internal function
        s->phase = PARDISO_CLEANUP;
        pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
                 &(s->n), &(s->fdum), s->KKT->p, s->KKT->i, &(s->idum), &(s->nrhs),
                 s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));

        // Check each attribute of the structure and free it if it exists
        if (s->KKT)       csc_spfree(s->KKT);
        if (s->KKT_i)     c_free(s->KKT_i);
        if (s->KKT_p)     c_free(s->KKT_p);
        if (s->bp)        c_free(s->bp);
        if (s->Pdiag_idx) c_free(s->Pdiag_idx);
        if (s->PtoKKT)    c_free(s->PtoKKT);
        if (s->AtoKKT)    c_free(s->AtoKKT);
        if (s->rhotoKKT)  c_free(s->rhotoKKT);

        c_free(s);

    }
}


// Initialize factorization structure
pardiso_solver *init_linsys_solver_pardiso(const csc * P, const csc * A, c_float sigma, c_float * rho_vec, c_int polish){
    c_int i;                     // loop counter
    c_int nnzKKT;                // Number of nonzeros in KKT
    // Define Variables
    pardiso_solver * s;          // Pardiso solver structure
    c_int n_plus_m;              // n_plus_m dimension

    // Size of KKT
    n_plus_m = P->m + A->m;

    // Allocate private structure to store KKT factorization
    s = c_calloc(1, sizeof(pardiso_solver));
    s->n = n_plus_m;

    // Working vector
    s->bp = c_malloc(sizeof(c_float) * n_plus_m);

    // Form KKT matrix
    if (polish){ // Called from polish()
        // Use s->bp for storing param2 = vec(delta)
        for (i = 0; i < A->m; i++){
            s->bp[i] = sigma;
        }

        s->KKT = form_KKT(P, A, 1, sigma, s->bp, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL);
    }
    else { // Called from ADMM algorithm

        // Allocate vectors of indices
        s->PtoKKT = c_malloc((P->p[P->n]) * sizeof(c_int));
        s->AtoKKT = c_malloc((A->p[A->n]) * sizeof(c_int));
        s->rhotoKKT = c_malloc((A->m) * sizeof(c_int));

        // Use s->bp for storing param2 = rho_inv_vec
        for (i = 0; i < A->m; i++){
            s->bp[i] = 1. / rho_vec[i];
        }

        s->KKT = form_KKT(P, A, 1, sigma, s->bp,
                          s->PtoKKT, s->AtoKKT,
                          &(s->Pdiag_idx), &(s->Pdiag_n), s->rhotoKKT);
    }

    // Check if matrix has been created
    if (!(s->KKT)) {
    #ifdef PRINTING
	    c_eprint("Error in forming KKT matrix");
    #endif
	    return OSQP_NULL;
    } else {
	    // Adjust indexing for Pardiso
	    nnzKKT = s->KKT->p[s->KKT->m];
	    s->KKT_i = c_malloc((nnzKKT) * sizeof(c_int));
	    s->KKT_p = c_malloc((s->KKT->m + 1) * sizeof(c_int));

	    for(i = 0; i < nnzKKT; i++){
	    	s->KKT_i[i] = s->KKT->i[i] + 1;
	    }
	    for(i = 0; i < n_plus_m+1; i++){
	    	s->KKT_p[i] = s->KKT->p[i] + 1;
	    }

    }

    // Set MKL interface layer (Long integers if activated)
    #ifdef DLONG
    mkl_set_interface_layer(MKL_INTERFACE_ILP64);
    #else
    mkl_set_interface_layer(MKL_INTERFACE_LP64);
    #endif

    // Set Pardiso variables
    s->mtype = -2;        // Real symmetric indefinite matrix
    s->nrhs = 1;          // Number of right hand sides
    s->maxfct = 1;        // Maximum number of numerical factorizations
    s->mnum = 1;          // Which factorization to use
    s->msglvl = 0;        // Do not print statistical information
    s->error = 0;         // Initialize error flag
    for ( i = 0; i < 64; i++ ){
        s->iparm[i] = 0;  // Setup Pardiso control parameters
        s->pt[i] = 0;     // Initialize the internal solver memory pointer
    }
    s->iparm[0] = 1;      // No solver default
    s->iparm[1] = 3;      // Fill-in reordering from OpenMP
    s->iparm[5] = 1;      // Write solution into b
    /* s->iparm[7] = 2;      // Max number of iterative refinement steps */
    s->iparm[7] = 0;      // Number of iterative refinement steps (auto, performs them only if perturbed pivots are obtained)
    s->iparm[9] = 13;     // Perturb the pivot elements with 1E-13
    s->iparm[34] = 0;     // Use Fortran-style indexing for indices
    /* s->iparm[34] = 1;     // Use C-style indexing for indices */

    // Print number of threads
    s->nthreads = mkl_get_max_threads();

    // Reordering and symbolic factorization
    s->phase = PARDISO_SYMBOLIC;
    pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n), s->KKT->x, s->KKT_p, s->KKT_i, &(s->idum), &(s->nrhs),
             s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));
    if ( s->error != 0 ){
        #ifdef PRINTING
            c_eprint("Error during symbolic factorization: %d", (int)s->error);
        #endif
        free_linsys_solver_pardiso(s);
        return OSQP_NULL;
    }

    // Numerical factorization
    s->phase = PARDISO_NUMERIC;
    pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n), s->KKT->x, s->KKT_p, s->KKT_i, &(s->idum), &(s->nrhs),
             s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));
    if ( s->error != 0 ){
        #ifdef PRINTING
            c_eprint("Error during numerical factorization: %d", (int)s->error);
        #endif
        free_linsys_solver_pardiso(s);
        return OSQP_NULL;
    }

    // Link Functions
    s->solve = &solve_linsys_pardiso;
    s->free = &free_linsys_solver_pardiso;
    s->update_matrices = &update_linsys_solver_matrices_pardiso;
    s->update_rho_vec = &update_linsys_solver_rho_vec_pardiso;

    // Assign type
    s->type = MKL_PARDISO_SOLVER;

    return s;
}

// Returns solution to linear system  Ax = b with solution stored in b
c_int solve_linsys_pardiso(pardiso_solver * s, c_float * b, const OSQPSettings *settings) {
    // Back substitution and iterative refinement
    s->phase = PARDISO_SOLVE;
    pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n), s->KKT->x, s->KKT_p, s->KKT_i, &(s->idum), &(s->nrhs),
             s->iparm, &(s->msglvl), b, s->bp, &(s->error));
    if ( s->error != 0 ){
        #ifdef PRINTING
        c_eprint("Error during solution: %d", (int)s->error);
        #endif
        return 1;
    }

    return 0;
}

// Update solver structure with new P and A
c_int update_linsys_solver_matrices_pardiso(pardiso_solver * s,
		const csc *P, const csc *A, const OSQPSettings *settings){

    // Update KKT matrix with new P
    update_KKT_P(s->KKT, P, s->PtoKKT, settings->sigma, s->Pdiag_idx, s->Pdiag_n);

    // Update KKT matrix with new A
    update_KKT_A(s->KKT, A, s->AtoKKT);

    // Perform numerical factorization
    s->phase = PARDISO_NUMERIC;
    pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n), s->KKT->x, s->KKT_p, s->KKT_i, &(s->idum), &(s->nrhs),
             s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));

    // Return exit flag
    return s->error;
}


c_int update_linsys_solver_rho_vec_pardiso(pardiso_solver * s, const c_float * rho_vec, const c_int m){
    c_int i;

    // Use s->bp for storing param2 = rho_inv_vec
    for (i = 0; i < m; i++){
        s->bp[i] = 1. / rho_vec[i];
    }

    // Update KKT matrix with new rho
    update_KKT_param2(s->KKT, s->bp, s->rhotoKKT, m);

    // Perform numerical factorization
    s->phase = PARDISO_NUMERIC;
    pardiso (s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n), s->KKT->x, s->KKT_p, s->KKT_i, &(s->idum), &(s->nrhs),
             s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));

    // Return exit flag
    return s->error;
}
