#include "kkt.h"

#ifndef EMBEDDED


csc* form_KKT(const csc  *P,
              const  csc *A,
              c_int       format,
              c_float     param1,
              c_float    *param2,
              c_int      *PtoKKT,
              c_int      *AtoKKT,
              c_int     **Pdiag_idx,
              c_int      *Pdiag_n,
              c_int      *param2toKKT) {
  c_int  nKKT, nnzKKTmax; // Size, number of nonzeros and max number of nonzeros
                          // in KKT matrix
  csc   *KKT_trip, *KKT;  // KKT matrix in triplet format and CSC format
  c_int  ptr, i, j;       // Counters for elements (i,j) and index pointer
  c_int  zKKT = 0;        // Counter for total number of elements in P and in
                          // KKT
  c_int *KKT_TtoC;        // Pointer to vector mapping from KKT in triplet form
                          // to CSC

  // Get matrix dimensions
  nKKT = P->m + A->m;

  // Get maximum number of nonzero elements (only upper triangular part)
  nnzKKTmax = P->p[P->n] + // Number of elements in P
              P->m +       // Number of elements in param1 * I
              A->p[A->n] + // Number of nonzeros in A
              A->m;        // Number of elements in - diag(param2)

  // Preallocate KKT matrix in triplet format
  KKT_trip = csc_spalloc(nKKT, nKKT, nnzKKTmax, 1, 1);

  if (!KKT_trip) return OSQP_NULL;  // Failed to preallocate matrix

  // Allocate vector of indices on the diagonal. Worst case it has m elements
  if (Pdiag_idx != OSQP_NULL) {
    (*Pdiag_idx) = c_malloc(P->m * sizeof(c_int));
    *Pdiag_n     = 0; // Set 0 diagonal elements to start
  }

  // Allocate Triplet matrices
  // P + param1 I
  for (j = 0; j < P->n; j++) { // cycle over columns
    // No elements in column j => add diagonal element param1
    if (P->p[j] == P->p[j + 1]) {
      KKT_trip->i[zKKT] = j;
      KKT_trip->p[zKKT] = j;
      KKT_trip->x[zKKT] = param1;
      zKKT++;
    }

    for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // cycle over rows
      // Get current row
      i = P->i[ptr];

      // Add element of P
      KKT_trip->i[zKKT] = i;
      KKT_trip->p[zKKT] = j;
      KKT_trip->x[zKKT] = P->x[ptr];

      if (PtoKKT != OSQP_NULL) PtoKKT[ptr] = zKKT;  // Update index from P to
                                                    // KKTtrip

      if (i == j) {                                 // P has a diagonal element,
                                                    // add param1
        KKT_trip->x[zKKT] += param1;

        // If index vector pointer supplied -> Store the index
        if (Pdiag_idx != OSQP_NULL) {
          (*Pdiag_idx)[*Pdiag_n] = ptr;
          (*Pdiag_n)++;
        }
      }
      zKKT++;

      // Add diagonal param1 in case
      if ((i < j) &&                  // Diagonal element not reached
          (ptr + 1 == P->p[j + 1])) { // last element of column j
        // Add diagonal element param1
        KKT_trip->i[zKKT] = j;
        KKT_trip->p[zKKT] = j;
        KKT_trip->x[zKKT] = param1;
        zKKT++;
      }
    }
  }

  if (Pdiag_idx != OSQP_NULL) {
    // Realloc Pdiag_idx so that it contains exactly *Pdiag_n diagonal elements
    (*Pdiag_idx) = c_realloc((*Pdiag_idx), (*Pdiag_n) * sizeof(c_int));
  }


  // A' at top right
  for (j = 0; j < A->n; j++) {                      // Cycle over columns of A
    for (ptr = A->p[j]; ptr < A->p[j + 1]; ptr++) {
      KKT_trip->p[zKKT] = P->m + A->i[ptr];         // Assign column index from
                                                    // row index of A
      KKT_trip->i[zKKT] = j;                        // Assign row index from
                                                    // column index of A
      KKT_trip->x[zKKT] = A->x[ptr];                // Assign A value element

      if (AtoKKT != OSQP_NULL) AtoKKT[ptr] = zKKT;  // Update index from A to
                                                    // KKTtrip
      zKKT++;
    }
  }

  // - diag(param2) at bottom right
  for (j = 0; j < A->m; j++) {
    KKT_trip->i[zKKT] = j + P->n;
    KKT_trip->p[zKKT] = j + P->n;
    KKT_trip->x[zKKT] = -param2[j];

    if (param2toKKT != OSQP_NULL) param2toKKT[j] = zKKT;  // Update index from
                                                          // param2 to KKTtrip
    zKKT++;
  }

  // Allocate number of nonzeros
  KKT_trip->nz = zKKT;

  // Convert triplet matrix to csc format
  if (!PtoKKT && !AtoKKT && !param2toKKT) {
    // If no index vectors passed, do not store KKT mapping from Trip to CSC/CSR
    if (format == 0) KKT = triplet_to_csc(KKT_trip, OSQP_NULL);
    else KKT = triplet_to_csr(KKT_trip, OSQP_NULL);
  }
  else {
    // Allocate vector of indices from triplet to csc
    KKT_TtoC = c_malloc((zKKT) * sizeof(c_int));

    if (!KKT_TtoC) return OSQP_NULL;  // Error in allocating KKT_TtoC vector

    // Store KKT mapping from Trip to CSC/CSR
    if (format == 0) KKT = triplet_to_csc(KKT_trip, KKT_TtoC);
    else KKT = triplet_to_csr(KKT_trip, KKT_TtoC);

    // Update vectors of indices from P, A, param2 to KKT (now in CSC format)
    if (PtoKKT != OSQP_NULL) {
      for (i = 0; i < P->p[P->n]; i++) {
        PtoKKT[i] = KKT_TtoC[PtoKKT[i]];
      }
    }

    if (AtoKKT != OSQP_NULL) {
      for (i = 0; i < A->p[A->n]; i++) {
        AtoKKT[i] = KKT_TtoC[AtoKKT[i]];
      }
    }

    if (param2toKKT != OSQP_NULL) {
      for (i = 0; i < A->m; i++) {
        param2toKKT[i] = KKT_TtoC[param2toKKT[i]];
      }
    }

    // Free mapping
    c_free(KKT_TtoC);
  }

  // Clean matrix in triplet format and return result
  csc_spfree(KKT_trip);

  return KKT;
}

#endif /* ifndef EMBEDDED */


#if EMBEDDED != 1

void update_KKT_P(csc          *KKT,
                  const csc    *P,
                  const c_int  *PtoKKT,
                  const c_float param1,
                  const c_int  *Pdiag_idx,
                  const c_int   Pdiag_n) {
  c_int i, j; // Iterations

  // Update elements of KKT using P
  for (i = 0; i < P->p[P->n]; i++) {
    KKT->x[PtoKKT[i]] = P->x[i];
  }

  // Update diagonal elements of KKT by adding sigma
  for (i = 0; i < Pdiag_n; i++) {
    j                  = Pdiag_idx[i]; // Extract index of the element on the
                                       // diagonal
    KKT->x[PtoKKT[j]] += param1;
  }
}

void update_KKT_A(csc *KKT, const csc *A, const c_int *AtoKKT) {
  c_int i; // Iterations

  // Update elements of KKT using A
  for (i = 0; i < A->p[A->n]; i++) {
    KKT->x[AtoKKT[i]] = A->x[i];
  }
}

void update_KKT_param2(csc *KKT, const c_float *param2,
                       const c_int *param2toKKT, const c_int m) {
  c_int i; // Iterations

  // Update elements of KKT using param2
  for (i = 0; i < m; i++) {
    KKT->x[param2toKKT[i]] = -param2[i];
  }
}

#endif // EMBEDDED != 1
