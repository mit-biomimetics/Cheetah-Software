#include "snopt_cwrap.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void snInit(snProblem* prob, char* name, char* prtfile, int summOn) {
  int leniw, lenrw, len;
  /*
   * snInit - call snInit to initialize workspace for SNOPT
   * On entry:
   *   prob     is the snProblem struct
   *   name     is the name of the problem
   *   prtfile  is the name of the output print file
   *            (empty string for no print file)
   *   summOn   is an integer indicating whether summary output
   *            (to screen) should be turned on (!= 0) or off (== 0)
   *
   * On exit:
   *   Internal workspace for SNOPT is initialized
   */
  init2zero(prob);

  leniw = 500;
  lenrw = 500;

  allocI(prob, leniw);
  allocR(prob, lenrw);

  prob->name   = name;

  len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("   SNOPT  C interface  2.0.0   ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_sninit(prtfile, len, summOn,
	   prob->iw, prob->leniw,
	   prob->rw, prob->lenrw);
  prob->initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void init2zero(snProblem* prob) {
  /*
   * init2zero - initializes internal variables to NULL/0
   */
  prob->name      = NULL;

  prob->memCalled  = 0;
  prob->initCalled = 0;

  prob->snLog     = NULL;
  prob->snLog2    = NULL;
  prob->sqLog     = NULL;
  prob->snSTOP    = NULL;

  prob->leniw     = 0;
  prob->lenrw     = 0;
  prob->iw        = NULL;
  prob->rw        = NULL;

  prob->leniu     = 0;
  prob->lenru     = 0;
  prob->iu        = NULL;
  prob->ru        = NULL;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void allocI(snProblem* prob, int len) {
  /*
   * allocI - allocate integer workspace for SNOPT
   * On entry:
   *   prob  is the snProblem struct
   *   len   is the desired length of the workspace
   */
  prob->leniw = len;
  prob->iw    = malloc(sizeof(int)*len);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void allocR(snProblem* prob, int len) {
  /*
   * allocR - allocate real workspace for SNOPT
   * On entry:
   *   prob  is the snProblem struct
   *   len   is the desired length of the workspace
   */
  prob->lenrw = len;
  prob->rw    = malloc(sizeof(double)*len);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocI (snProblem* prob, int len) {
  /*
   * reallocI - reallocate integer workspace for SNOPT
   * On entry:
   *   prob  is the snProblem struct
   *   len   is the desired length of the workspace
   */
  prob->leniw = len;
  prob->iw    = (int*)realloc(prob->iw, sizeof(int)*prob->leniw);

  setIntParameter(prob, (char*)"Total int workspace", prob->leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocR(snProblem* prob, int len) {
  /*
   * reallocR - reallocate real workspace for SNOPT
   * On entry:
   *   prob  is the snProblem struct
   *   len   is the desired length of the workspace
   */
  prob->lenrw = len;
  prob->rw = (double*)realloc(prob->rw, sizeof(double)*prob->lenrw);

  setIntParameter(prob, (char*)"Total real workspace", prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setPrintfile(snProblem* prob, char *prtname) {
  /*
   * setPrintfile - set name of the output print file
   * On entry:
   *   prob     is the snProblem struct
   *   prtname  is the name of the file
   */
  int len = strlen(prtname);

  assert(prob->initCalled == 1);
  f_snsetprint(prtname, len,
		prob->iw, prob->leniw, prob->rw, prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int setSpecsfile(snProblem* prob, char *spcname) {
  /*
   * setPrintfile - call snSpecs to read in options file
   * On entry:
   *   prob     is the snProblem struct
   *   spcname  is the name of the file
   * On exit:
   *   Returns integer info code
   *     (see SNOPT documentation for snSpecs)
   */
  int inform;
  int len = strlen(spcname);

  assert(prob->initCalled == 1);
  f_snspec(spcname, len, &inform,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int setParameter(snProblem* prob, char stropt[]) {
  /*
   * setParameter - set given SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option and its value
   * On exit:
   *   Returns number of errors encountered
   *     (see SNOPT documentation for snSet)
   */
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_snset(stropt, len, &errors,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int getParameter(snProblem* prob, char stropt[], char strout[]) {
  /*
   * getParameter - returns SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option and its value
   * On exit:
   *   Returns number of errors encountered
   *   strout   contains the option's value
   *     (see SNOPT documentation for snGetC)
   */
  int errors;
  int inlen  = strlen(stropt);
  int outlen = strlen(strout);

  assert(prob->initCalled == 1);
  f_sngetc(stropt, inlen, strout, outlen, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int setIntParameter(snProblem* prob, char stropt[], int opt) {
  /*
   * setIntParameter - set an integer SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option keyword
   *   opt      is the integer containg the option value
   * On exit:
   *   Returns number of errors encountered
   *     (see SNOPT documentation for snSetI)
   */
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_snseti (stropt, len, opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int getIntParameter(snProblem* prob, char stropt[], int opt) {
  /*
   * getIntParameter - get an integer SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option keyword
   * On exit:
   *   Returns number of errors encountered
   *   opt      is the integer containg the option value
   *     (see SNOPT documentation for snGetI)
   */
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sngeti(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int setRealParameter(snProblem* prob, char stropt[], double opt) {
  /*
   * setRealParameter - set a real SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option keyword
   *   opt      is the double containg the option value
   * On exit:
   *   Returns number of errors encountered
   *     (see SNOPT documentation for snSetR)
   */
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_snsetr (stropt, len, opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int getRealParameter(snProblem* prob, char stropt[], double opt) {
  /*
   * getRealParameter - get a real SNOPT option/parameter
   * On entry:
   *   prob     is the snProblem struct
   *   stropt   is the string containing the option keyword
   * On exit:
   *   Returns number of errors encountered
   *   opt      is the double containg the option value
   *     (see SNOPT documentation for snGetR)
   */
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sngetr (stropt, len, &opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setUserI(snProblem* prob, int *iu, int leniu) {
  /*
   * setUserI - sets user integer workspace
   * On entry:
   *   prob       is the snProblem struct
   *   iu, leniu  is the integer workspace and its length
   */
  prob->iu    = iu;
  prob->leniu = leniu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setUserR(snProblem* prob, double *ru, int lenru) {
  /*
   * setUserR - sets user real workspace
   * On entry:
   *   prob       is the snProblem struct
   *   ru, lenru  is the double workspace and its length
   */
  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setUserspace(snProblem* prob, int *iu, int leniu, double *ru, int lenru) {
  /*
   * setUserspace - sets user integer and real workspace
   * On entry:
   *   prob       is the snProblem struct
   *   iu, leniu  is the integer workspace and its length
   *   ru, lenru  is the double workspace and its length
   */
  prob->iu    = iu;
  prob->leniu = leniu;

  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setWorkspace(snProblem* prob, int m, int n, int ne,
		  int negCon, int nnCon, int nnObj, int nnJac) {
  /*
   * setWorkspace - call snMem to determine required amount of
   *                SNOPT integer and real workspace;
   *                workspace is automatically (re-)allocated
   * On entry:
   *   prob       is the snProblem struct
   *   m, n, ne   are integers indicating the number of constraints,
   *              variables and nonzero elements in the Jacobian
   *   negCon     is the number of nonzero elements in the nonlinear
   *              Jacobian gCon
   *   nnCon      is the number of nonlinear constraints
   *   nnObj      is the number of nonlinear objective variables
   *   nnJac      is the number of nonlinear Jacobian variables
   * On exit:
   *   Workspace is allocated to size determined by SNOPT.
   */
  int miniw, minrw, inform, ineG;
  int memGuess = 0;

  assert(prob->initCalled == 1);

  if (negCon < 0) {
    ineG = nnCon*nnJac;
    memGuess = 1;
  } else {
    ineG = negCon;
  }

  f_snmem(&inform, m, n, ne,
	  ineG, nnCon, nnObj, nnJac, &miniw, &minrw,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  if (miniw > prob->leniw) { reallocI(prob, miniw); }
  if (minrw > prob->lenrw) { reallocR(prob, minrw); }

  prob->memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setWorkspaceA(snProblem* prob, int nF, int n, int neA, int neG) {
  /*
   * setWorkspaceA - call snMemA to determine required amount of
   *                 SNOPTA integer and real workspace;
   *                 workspace is automatically (re-)allocated
   * On entry:
   *   prob       is the snProblem struct
   *   nF, n      are integers indicating the number of constraints,
   *              variables
   *   neA        is the number of nonzero elements in the linear
   *              Jacobian matrix A
   *   neG        is the number of nonzero elements in the nonlinear
   *              Jacobian G
   * On exit:
   *   Workspace is allocated to size determined by SNOPT.
   */
  int miniw, minrw, inform, ineG;
  int memGuess = 0;

  assert(prob->initCalled == 1);

  if (neG < 0) {
    ineG     = n*nF;
    memGuess = 1;
  } else {
    ineG = neG;
  }

  f_snmema(&inform, nF, n, neA, ineG,
	   &miniw, &minrw,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  if (miniw > prob->leniw) { reallocI(prob, miniw); }
  if (minrw > prob->lenrw) { reallocR(prob, minrw); }

  prob->memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setLog(snProblem* prob, isnLog snLog, isnLog2 snLog2, isqLog sqLog) {
  /*
   * setLog - sets user-defined log routines for SNOPT
   * On entry:
   *   prob    is the snProblem struct
   *   snLog   is the major iteration log print routine
   *   snLog2  is the minor iteration log print routine
   *   sqLog   is the QP    iteration log print routine
   */
  prob->snLog  = snLog;
  prob->snLog2 = snLog2;
  prob->sqLog  = sqLog;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void setSTOP(snProblem* prob, isnSTOP snSTOP) {
  /*
   * setSTOP - set the snSTOP routine for SNOPT
   * On entry:
   *   prob    is the snProblem struct
   *   snSTOP  is the user-defined routine called at every major iteration
   */
  prob->snSTOP = snSTOP;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int snoptA(snProblem* prob, int start,
	   int nF, int n, double ObjAdd, int ObjRow,
	   snFunA usrfun,
	   int neA, int *iAfun, int *jAvar, double *A,
	   int neG, int *iGfun, int *jGvar,
	   double *xlow, double *xupp, double *Flow, double *Fupp,
	   double *x, int *xstate, double *xmul,
	   double *F, int *Fstate, double *Fmul,
	   int* nS, int* nInf, double* sInf) {
  /*
   * snoptA - call SNOPTA to solve the problem
   * See SNOPT documentation on SNOPTA.
   */

  int i, inform, iObj, miniw, minrw;

  assert(prob->initCalled == 1);

  if (prob->memCalled == 0) { setWorkspaceA(prob, nF, n, neA, neG); }

  for (i = 0; i < neA; i++) {
    iAfun[i]++;
    jAvar[i]++;
  }
  for (i = 0; i < neG; i++) {
    iGfun[i]++;
    jGvar[i]++;
  }


  iObj = ObjRow+1;

  f_snkera(start, prob->name, nF, n, ObjAdd, iObj, usrfun,
	   prob->snLog, prob->snLog2, prob->sqLog, prob->snSTOP,
	   iAfun, jAvar, neA, A,
	   iGfun, jGvar, neG,
	   xlow, xupp, Flow, Fupp,
	   x, xstate, xmul,
	   F, Fstate, Fmul,
	   &inform, nS, nInf, sInf,
	   &miniw, &minrw,
	   prob->iu, prob->leniu, prob->ru, prob->lenru,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < neA; i++) {
    iAfun[i]--;
    jAvar[i]--;
  }
  for (i = 0; i < neG; i++) {
    iGfun[i]--;
    jGvar[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int snJac( snProblem* prob,
	   int nF, int n, snFunA usrfun,
	   double *x, double *xlow, double *xupp,
	   int* neA, int iAfun[], int jAvar[], double A[],
	   int* neG, int iGfun[], int jGvar[] ) {
  /*
   * snJac - compute Jacobian structure for SNOPTA
   * See SNOPT documentation on snJac.
   */

  int i, inform, lenA, lenG, miniw, minrw;

  lenA = n*nF;
  lenG = n*nF;

  if (prob->memCalled == 0) { setWorkspaceA(prob, nF, n, lenA, lenG); }

  f_snjac(&inform, nF, n, usrfun, x, xlow, xupp,
	  iAfun, jAvar, lenA, neA, A,
	  iGfun, jGvar, lenG, neG,
	  &miniw, &minrw,
	  prob->iu, prob->leniu, prob->ru, prob->lenru,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < *neA; i++) {
    iAfun[i]--;
    jAvar[i]--;
  }
  for (i = 0; i < *neG; i++) {
    iGfun[i]--;
    jGvar[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int snoptB(snProblem* prob, int start, int m, int n, int ne,
	   int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	   snConB funcon, snObjB funobj,
	   double *valJ, int *indJ, int *locJ,
	   double *bl, double *bu, int *hs, double *x,
	   double *pi, double *rc, double* objective,
	   int* nS, int* nInf, double* sInf) {
  /*
   * snoptB - call SNOPTB to solve the problem
   * See SNOPT documentation on SNOPTB.
   */

  int i, inform, iiObj, miniw, minrw;

  assert(prob->initCalled == 1);

  if (prob->memCalled == 0) {
    setWorkspace(prob, m, n, ne, -1, nnCon, nnObj, nnJac);
  }


  for (i = 0; i < ne; i++) {
    indJ[i]++;
  }
  for (i = 0; i <= n; i++) {
    locJ[i]++;
  }

  iiObj = iObj+1;

  f_snkerb(start, prob->name, m, n, ne,
	   nnCon, nnObj, nnJac, iiObj,
	   ObjAdd,
	   funcon, funobj,
	   prob->snLog, prob->snLog2, prob->sqLog, prob->snSTOP,
	   valJ, indJ, locJ,
	   bl, bu, hs, x, pi, rc,
	   &inform, nS, nInf, sInf, objective,
	   &miniw, &minrw,
	   prob->iu, prob->leniu, prob->ru, prob->lenru,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < ne; i++) {
    indJ[i]--;
  }
  for (i = 0; i <= n; i++) {
    locJ[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int snoptC(snProblem* prob, int start, int m, int n, int ne,
	   int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	   snFunC usrfun,
	   double *valJ, int *indJ, int *locJ,
	   double *bl, double *bu, int *hs, double *x,
	   double *pi, double *rc, double* objective,
	   int* nS, int* nInf, double* sInf) {
  /*
   * snoptC - call SNOPTC to solve the problem
   * See SNOPT documentation on SNOPTC.
   */

  int i, inform, iiObj, miniw, minrw;

  assert(prob->initCalled == 1);

  if (prob->memCalled == 0) {
    setWorkspace(prob, m, n, ne, -1, nnCon, nnObj, nnJac);
  }

  for (i = 0; i < ne; i++) {
    indJ[i]++;
  }
  for (i = 0; i <= n; i++) {
    locJ[i]++;
  }
  iiObj = iObj+1;

  f_snkerc(start, prob->name, m, n, ne,
	   nnCon, nnObj, nnJac, iiObj, ObjAdd,
	   usrfun,
	   prob->snLog, prob->snLog2, prob->sqLog, prob->snSTOP,
	   valJ, indJ, locJ,
	   bl, bu, hs, x, pi, rc,
	   &inform, nS, nInf, sInf, objective,
	   &miniw, &minrw,
	   prob->iu, prob->leniu, prob->ru, prob->lenru,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < ne; i++) {
    indJ[i]--;
  }
  for (i = 0; i <= n; i++) {
    locJ[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int solveA(snProblem* prob, int start,
	   int nF, int n, double ObjAdd, int ObjRow,
	   snFunA usrfun,
	   int neA, int *iAfun, int *jAvar, double *A,
	   int neG, int *iGfun, int *jGvar,
	   double *xlow, double *xupp, double *Flow, double *Fupp,
	   double *x, int *xstate, double *xmul,
	   double *F, int *Fstate, double *Fmul,
	   int* nS, int* nInf, double* sInf) {
  /*
   * solveA - call SNOPTA to solve the problem
   * See SNOPT documentation on SNOPTA.
   */

  int inform;
  inform = snoptA( prob, start, nF, n, ObjAdd, ObjRow, usrfun,
		   neA, iAfun, jAvar, A,
		   neG, iGfun, jGvar,
		   xlow, xupp, Flow, Fupp,
		   x, xstate, xmul,
		   F, Fstate, Fmul,
		   nS, nInf, sInf );
  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int solveB(snProblem* prob, int start, int m, int n, int ne,
	   int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	   snConB funcon, snObjB funobj,
	   double *valJ, int *indJ, int *locJ,
	   double *bl, double *bu, int *hs, double *x,
	   double *pi, double *rc, double* objective,
	   int* nS, int* nInf, double* sInf) {
  /*
   * solveB - call SNOPTB to solve the problem
   * See SNOPT documentation on SNOPTB.
   */

  int inform;
  inform = snoptB( prob, start, m, n, ne,
		   nnCon, nnObj, nnJac, iObj, ObjAdd,
		   funcon, funobj, valJ, indJ, locJ,
		   bl, bu, hs, x, pi, rc, objective,
		   nS, nInf, sInf);
  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int solveC(snProblem* prob, int start, int m, int n, int ne,
	   int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	   snFunC usrfun,
	   double *valJ, int *indJ, int *locJ,
	   double *bl, double *bu, int *hs, double *x,
	   double *pi, double *rc, double* objective,
	   int* nS, int* nInf, double* sInf) {
  /*
   * solveC - call SNOPTC to solve the problem
   * See SNOPT documentation on SNOPTC.
   */

  int inform;
  inform = snoptC( prob, start, m, n, ne,
		   nnCon, nnObj, nnJac, iObj, ObjAdd,
		   usrfun, valJ, indJ, locJ,
		   bl, bu, hs, x, pi, rc, objective,
		   nS, nInf, sInf);
  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void deleteSNOPT(snProblem* prob) {
  /*
   * deleteSNOPT - frees internal memory associated with SNOPT
   */
  f_snend(prob->iw, prob->leniw, prob->rw, prob->lenrw);

  free(prob->iw);
  free(prob->rw);

  init2zero(prob);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
