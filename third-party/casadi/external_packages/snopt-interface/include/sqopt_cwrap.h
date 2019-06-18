#ifndef SQOPT_C_H
#define SQOPT_C_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "snopt.h"

/* File: sqopt_cwrap.h
 *   C interface for SQOPT.
 */
typedef struct {
  char   *name;

  int     memCalled;
  int     initCalled;
  isqLog  sqLog;

  int     lenrw, leniw;
  int    *iw;
  double *rw;

  int     lenru, leniu;
  int    *iu;
  double *ru;

} sqProblem;

void sqInit         (sqProblem* prob, char* name, char* prtfile, int summOn);
void init2zeroQ     (sqProblem* prob);

void allocIQ        (sqProblem* prob, int len);
void allocRQ        (sqProblem* prob, int len);
void reallocIQ      (sqProblem* prob, int len);
void reallocRQ      (sqProblem* prob, int len);

void sqPrintfile   (sqProblem* prob, char* prtname);
int  sqSpecsfile   (sqProblem* prob, char* spcname);

int sqSetParameter    (sqProblem* prob, char stropt[]);
int sqGetParameter    (sqProblem* prob, char stropt[], char strout[]);
int sqSetIntParameter (sqProblem* prob, char stropt[], int opt);
int sqGetIntParameter (sqProblem* prob, char stropt[], int opt);
int sqSetRealParameter(sqProblem* prob, char stropt[], double opt);
int sqGetRealParameter(sqProblem* prob, char stropt[], double opt);

void sqUserI       (sqProblem* prob, int *iu, int leniu);
void sqUserR       (sqProblem* prob, double *ru, int lenru);
void sqUserspace   (sqProblem* prob, int *iu, int leniu, double *ru, int lenru);
void sqSetLog      (sqProblem* prob, isqLog sqLog);
void sqWorkspace   (sqProblem* prob, int m, int n, int neA, int ncObj, int nnH);

int solveQ(sqProblem* prob, int start, sqFunHx qpHx,
	   int m, int n, int neA, int ncObj, int nnH,
	   int iObj, double ObjAdd,
	   double *valA, int *indA, int *locA,
	   double *bl, double *bu, double *cObj,
	   int *eType, int *hs, double *x, double *pi, double *rc,
	   double* objective,
	   int* nS, int* nInf, double* sInf);

int sqopt(sqProblem* prob, int start, sqFunHx qpHx,
	  int m, int n, int neA, int ncObj, int nnH,
	  int iObj, double ObjAdd,
	  double *valA, int *indA, int *locA,
	  double *bl, double *bu, double *cObj,
	  int *eType, int *hs, double *x, double *pi, double *rc,
	  double* objective,
	  int* nS, int* nInf, double* sInf);

void deleteSQOPT(sqProblem* prob);

#endif /* SQOPT_C_H */
