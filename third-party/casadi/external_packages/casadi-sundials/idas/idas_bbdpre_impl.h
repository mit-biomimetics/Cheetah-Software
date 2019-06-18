/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file (private version) for the IDABBDPRE
 * module, for a band-block-diagonal preconditioner, i.e. a
 * block-diagonal matrix with banded blocks, for use with IDAS
 * and an IDASPILS linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _IDASBBDPRE_IMPL_H
#define _IDASBBDPRE_IMPL_H

#include <idas/idas_bbdpre.h>
#include <sundials/sundials_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Definition of IBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct IBBDPrecDataRec {

  /* passed by user to IDABBDPrecAlloc and used by
     IDABBDPrecSetup/IDABBDPrecSolve functions */

  long int mudq, mldq, mukeep, mlkeep;
  realtype rel_yy;
  IDABBDLocalFn glocal;
  IDABBDCommFn gcomm;

  /* allocated for use by IDABBDPrecSetup */

  N_Vector tempv4;

  /* set by IDABBDPrecon and used by IDABBDPrecSolve */

  DlsMat PP;
  long int *lpivots;

  /* set by IDABBDPrecAlloc and used by IDABBDPrecSetup */

  long int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to ida_mem */

  void *ida_mem;

} *IBBDPrecData;

/*
 * -----------------------------------------------------------------
 * Type: IDABBDPrecDataB
 * -----------------------------------------------------------------
 */

typedef struct IDABBDPrecDataRecB{

  /* BBD user functions (glocB and cfnB) for backward run */
  IDABBDLocalFnB glocalB;
  IDABBDCommFnB  gcommB;
    
  /* BBD prec data */
  /* //!void *bbd_dataB; */

} *IDABBDPrecDataB;


/*
 * -----------------------------------------------------------------
 * IDABBDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBBD_MEM_NULL    "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL    "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. IDABBDPrecInit must be called."
#define MSGBBD_FUNC_FAILED "The Glocal or Gcomm routine failed in an unrecoverable manner."

#define MSGBBD_AMEM_NULL   "idaadj_mem = NULL illegal."
#define MSGBBD_PDATAB_NULL "IDABBDPRE memory is NULL for the backward integration."
#define MSGBBD_BAD_T       "Bad t for interpolation."

#ifdef __cplusplus
}
#endif

#endif
