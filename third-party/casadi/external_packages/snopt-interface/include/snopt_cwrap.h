#ifndef SNOPT_C_H
#define SNOPT_C_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "snopt.h"

/* File: snopt_cwrap.h
 *   C interface for SNOPTB and SNOPTC.
 */
typedef struct {
  char   *name;

  int     memCalled;
  int     initCalled;

  isnSTOP snSTOP;
  isnLog  snLog;
  isnLog2 snLog2;
  isqLog  sqLog;

  int     lenrw, leniw;
  int    *iw;
  double *rw;

  int     lenru, leniu;
  int    *iu;
  double *ru;

} snProblem;

void snInit         ( snProblem* prob, char* name, char* prtfile, int summOn );
void init2zero      ( snProblem* prob );

void allocI         ( snProblem* prob, int len );
void allocR         ( snProblem* prob, int len );
void reallocI       ( snProblem* prob, int len );
void reallocR       ( snProblem* prob, int len );

void setPrintfile   ( snProblem* prob, char* prtname );
int  setSpecsfile   ( snProblem* prob, char* spcname );

int setParameter    ( snProblem* prob, char stropt[] );
int getParameter    ( snProblem* prob, char stropt[], char strout[] );
int setIntParameter ( snProblem* prob, char stropt[], int opt );
int getIntParameter ( snProblem* prob, char stropt[], int opt );
int setRealParameter( snProblem* prob, char stropt[], double opt );
int getRealParameter( snProblem* prob, char stropt[], double opt );

void setUserI       ( snProblem* prob, int *iu, int leniu );
void setUserR       ( snProblem* prob, double *ru, int lenru );
void setUserspace   ( snProblem* prob, int *iu, int leniu,
		      double *ru, int lenru );

void setLog         ( snProblem* prob, isnLog snLog, isnLog2 snLog2, isqLog sqLog );
void setSTOP        ( snProblem* prob, isnSTOP snSTOP );

void setWorkspace   ( snProblem* prob, int m, int n, int ne,
		      int negCon, int nnCon, int nnObj, int nnJac);
void setWorkspaceA  ( snProblem* prob, int nF, int n, int neA, int neG);

int snJac( snProblem* prob,
	   int nF, int n, snFunA usrfun,
	   double *x, double *xlow, double *xupp,
	   int *neA, int *iAfun, int *jAvar, double *A,
	   int *neG, int *iGfun, int *jGvar );

int solveA( snProblem* prob, int start,
	    int nF, int n, double ObjAdd, int ObjRow,
	    snFunA usrfun,
	    int neA, int *iAfun, int *jAvar, double *A,
	    int neG, int *iGfun, int *jGvar,
	    double *xlow, double *xupp, double *Flow, double *Fupp,
	    double *x, int *xstate, double *xmul,
	    double *F, int *Fstate, double *Fmul,
	    int* nS, int* nInf, double* sInf );

int snoptA( snProblem* prob, int start,
	    int nF, int n, double ObjAdd, int ObjRow,
	    snFunA usrfun,
	    int neA, int *iAfun, int *jAvar, double *A,
	    int neG, int *iGfun, int *jGvar,
	    double *xlow, double *xupp, double *Flow, double *Fupp,
	    double *x, int *xstate, double *xmul,
	    double *F, int *Fstate, double *Fmul,
	    int* nS, int* nInf, double* sInf );

int solveB( snProblem* prob, int start, int m, int n, int ne,
	    int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	    snConB funcon, snObjB funobj,
	    double *valJ, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs, double *x,
	    double *pi, double *rc, double* objective,
	    int* nS, int* nInf, double* sInf );

int snoptB( snProblem* prob, int start, int m, int n, int ne,
	    int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	    snConB funcon, snObjB funobj,
	    double *valJ, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs, double *x,
	    double *pi, double *rc, double* objective,
	    int* nS, int* nInf, double* sInf );

int solveC( snProblem* prob, int start, int m, int n, int ne,
	    int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	    snFunC usrfun,
	    double *valJ, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs, double *x,
	    double *pi, double *rc, double* objective,
	    int* nS, int* nInf, double* sInf );

int snoptC( snProblem* prob, int start, int m, int n, int ne,
	    int nnCon, int nnObj, int nnJac, int iObj, double ObjAdd,
	    snFunC usrfun,
	    double *valJ, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs, double *x,
	    double *pi, double *rc, double* objective,
	    int* nS, int* nInf, double* sInf );

void deleteSNOPT    ( snProblem* prob );

#endif /* SNOPT_C_H */
