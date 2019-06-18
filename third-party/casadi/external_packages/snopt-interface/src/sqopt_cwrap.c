#include "sqopt_cwrap.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqInit(sqProblem* prob, char* name, char* prtfile, int summOn) {
  int leniw, lenrw, len;

  init2zeroQ(prob);

  leniw = 500;
  lenrw = 500;

  allocIQ(prob, leniw);
  allocRQ(prob, lenrw);

  prob->name   = name;

  len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("   SQOPT  C interface  2.0.0   ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_sqinit(prtfile, len, summOn,
	   prob->iw, prob->leniw,
	   prob->rw, prob->lenrw);
  prob->initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void init2zeroQ(sqProblem* prob) {
  prob->name       = NULL;

  prob->memCalled  = 0;
  prob->initCalled = 0;

  prob->sqLog      = NULL;

  prob->leniw      = 0;
  prob->lenrw      = 0;
  prob->iw         = NULL;
  prob->rw         = NULL;

  prob->leniu      = 0;
  prob->lenru      = 0;
  prob->iu         = NULL;
  prob->ru         = NULL;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void allocIQ(sqProblem* prob, int len) {
  prob->leniw = len;
  prob->iw    = malloc(sizeof(int)*len);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void allocRQ(sqProblem* prob, int len) {
  prob->lenrw = len;
  prob->rw    = malloc(sizeof(double)*len);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocIQ(sqProblem* prob, int len) {
  prob->leniw = len;
  prob->iw    = (int*)realloc(prob->iw, sizeof(int)*prob->leniw);

  sqSetIntParameter(prob, (char*)"Total int workspace", prob->leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocRQ(sqProblem* prob, int len) {
  prob->lenrw = len;
  prob->rw = (double*)realloc(prob->rw, sizeof(double)*prob->lenrw);

  sqSetIntParameter(prob, (char*)"Total real workspace", prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqPrintfile(sqProblem* prob, char *prtname) {
  int len = strlen(prtname);

  assert(prob->initCalled == 1);
  f_sqsetprint(prtname, len,
		prob->iw, prob->leniw, prob->rw, prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqSpecsfile(sqProblem* prob, char *spcname) {
  int inform;
  int len = strlen(spcname);

  assert(prob->initCalled == 1);
  f_sqspec(spcname, len, &inform,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqSetParameter(sqProblem* prob, char stropt[]) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sqset(stropt, len, &errors,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqGetParameter(sqProblem* prob, char stropt[], char strout[]) {
  int errors;
  int inlen  = strlen(stropt);
  int outlen = strlen(strout);

  assert(prob->initCalled == 1);
  f_sqgetc(stropt, inlen, strout, outlen, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqSetIntParameter(sqProblem* prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sqseti (stropt, len, opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqGetIntParameter(sqProblem* prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sqgeti(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqSetRealParameter(sqProblem* prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sqsetr (stropt, len, opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqGetRealParameter(sqProblem* prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_sqgetr (stropt, len, &opt, &errors,
	     prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqUserI(sqProblem* prob, int *iu, int leniu) {
  prob->iu    = iu;
  prob->leniu = leniu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqUserR(sqProblem* prob, double *ru, int lenru) {
  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqUserspace(sqProblem* prob, int *iu, int leniu, double *ru, int lenru) {
  prob->iu    = iu;
  prob->leniu = leniu;

  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqWorkspace(sqProblem* prob, int m, int n,
		 int neA, int ncObj, int nnH) {
  int miniw, minrw, inform;

  assert(prob->initCalled == 1);

  f_sqmem(&inform, m, n, neA, ncObj, nnH,
	  &miniw, &minrw,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  if (miniw > prob->leniw) { reallocIQ(prob, miniw); }
  if (minrw > prob->lenrw) { reallocRQ(prob, minrw); }

  prob->memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void sqSetLog(sqProblem* prob, isqLog sqLog) {
  prob->sqLog  = sqLog;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int solveQ(sqProblem* prob, int start, sqFunHx qpHx,
	   int m, int n, int neA, int ncObj, int nnH,
	   int iObj, double ObjAdd,
	   double *valA, int *indA, int *locA,
	   double *bl, double *bu, double *cObj,
	   int *eType, int *hs, double *x, double *pi, double *rc,
	   double* objective,
	   int* nS, int* nInf, double* sInf) {

  int i, inform, iiObj, miniw, minrw;

  assert(prob->initCalled == 1);

  if (prob->memCalled == 0) {
    sqWorkspace(prob, m, n, neA, ncObj, nnH);
  }

  for (i = 0; i < neA; i++) {
    indA[i]++;
  }
  for (i = 0; i <= n; i++) {
    locA[i]++;
  }
  iiObj = iObj+1;

  f_snkerq(start, qpHx, prob->sqLog,
	   m, n, neA, ncObj, nnH,
	   iiObj, ObjAdd, prob->name,
	   valA, indA, locA,
	   bl, bu, cObj, eType, hs, x, pi, rc,
	   &inform, nS, nInf, sInf, objective,
	   &miniw, &minrw,
	   prob->iu, prob->leniu, prob->ru, prob->lenru,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < neA; i++) {
    indA[i]--;
  }
  for (i = 0; i <= n; i++) {
    locA[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int sqopt(sqProblem* prob, int start, sqFunHx qpHx,
	  int m, int n, int neA, int ncObj, int nnH,
	  int iObj, double ObjAdd,
	  double *valA, int *indA, int *locA,
	  double *bl, double *bu, double *cObj,
	  int *eType, int *hs, double *x, double *pi, double *rc,
	  double* objective,
	  int* nS, int* nInf, double* sInf) {

  int i, inform, iiObj, miniw, minrw;

  assert(prob->initCalled == 1);

  if (prob->memCalled == 0) {
    sqWorkspace(prob, m, n, neA, ncObj, nnH);
  }

  for (i = 0; i < neA; i++) {
    indA[i]++;
  }
  for (i = 0; i <= n; i++) {
    locA[i]++;
  }
  iiObj = iObj+1;

  f_snkerq(start, qpHx, prob->sqLog,
	   m, n, neA, ncObj, nnH,
	   iiObj, ObjAdd, prob->name,
	   valA, indA, locA,
	   bl, bu, cObj, eType, hs, x, pi, rc,
	   &inform, nS, nInf, sInf, objective,
	   &miniw, &minrw,
	   prob->iu, prob->leniu, prob->ru, prob->lenru,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  for (i = 0; i < neA; i++) {
    indA[i]--;
  }
  for (i = 0; i <= n; i++) {
    locA[i]--;
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void deleteSQOPT(sqProblem* prob) {
  f_sqend(prob->iw, prob->leniw, prob->rw, prob->lenrw);

  free(prob->iw);
  free(prob->rw);

  init2zeroQ(prob);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
