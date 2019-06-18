#ifndef SNOPTPROBLEM_H
#define SNOPTPROBLEM_H

#include "snopt.h"

/* File snoptProblem.hpp
 *   C++ interface for SNOPT and SQOPT
 *
 */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblem {
protected:
  snoptProblem();
  snoptProblem(const char*name);
  ~snoptProblem();

  void init2zero();

  char    Prob[30];

  int     initCalled, memCalled;

  int     leniw, lenrw;
  double *rw;
  int    *iw;

  int     lenru, leniu;
  double *ru;
  int    *iu;

  void allocI    (int leniw);
  void allocR    (int lenrw);
  void reallocI  (int leniw);
  void reallocR  (int lenrw);

public:
  void setProbName    (const char *Prob);
  void setPrintFile   (const char *prtname);

  int getParameter    (const char *stroptin, char *stroptout);
  int getIntParameter (const char *stropt,   int    &opt);
  int getRealParameter(const char *stropt,   double &opt);
  int setParameter    (const char *stroptin);
  int setIntParameter (const char *stropt,   int     opt);
  int setRealParameter(const char *stropt,   double  opt);

  void setUserI       (int    *iu, int leniu);
  void setUserR       (double *ru, int lenru);
  void setUserspace   (int    *iu, int leniu, double *ru, int lenru);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemABC : public snoptProblem  {
protected:
  snoptProblemABC();
  snoptProblemABC(const char*name);
  ~snoptProblemABC();

  void init2zero();

  isnLog  snLog;
  isnLog2 snLog2;
  isqLog  sqLog;
  isnSTOP snSTOP;

public:
  void initialize     (const char *prtfile, int summOn);
  int  setSpecsFile   (const char *specname);

  void setLog         (isnLog snLog, isnLog2 snLog2, isqLog sqLog);
  void setSTOP        (isnSTOP snSTOP);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemA : public snoptProblemABC {
public:
  snoptProblemA();
  snoptProblemA(const char *name);
  ~snoptProblemA();

  void setWorkspace(int neF, int n, int neA, int neG);
  int  computeJac(int neF, int n, snFunA usrfunA,
		  double *x, double *xlow, double*xupp,
		  int *&iAfun, int *&jAvar, double *&A, int &neA,
		  int *&iGfun, int *&jGvar, int &neG);

  int  solve(int starttype, int nF, int n, double ObjAdd,
	     int ObjRow, snFunA usrfunA,
	     double *xlow, double *xupp, double *Flow, double *Fupp,
	     double *x, int *xstate, double *xmul,
	     double *F, int *Fstate, double *Fmul,
	     int &nS, int &nInf, double &sInf);

  int  solve(int starttype, int nF, int n, double ObjAdd,
	     int ObjRow, snFunA usrfunA,
	     int *iAfun, int *jAvar, double *A, int neA,
	     int *iGfun, int *jGvar, int neG,
	     double *xlow, double *xupp, double *Flow, double *Fupp,
	     double *x, int *xstate, double *xmul,
	     double *F, int *Fstate, double *Fmul,
	     int &nS, int &nInf, double &sInf);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemC : public snoptProblemABC {
public:
  snoptProblemC();
  snoptProblemC(const char*name);
  ~snoptProblemC();

  void setWorkspace(int m, int n, int ne, int negCon,
		    int nnCon, int nnObj, int nnJac);

  int solve(int starttype, int m, int n, int ne, int negCon,
	    int nnCon, int nnObj, int nnJac,
	    int iObj, double ObjAdd, snFunC usrfunC,
	    double *Jval, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs,
	    double *x, double *pi, double *rc,
	    int &nS, int &nInf, double &sInf, double &objective);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class snoptProblemB : public snoptProblemC {
public:
  snoptProblemB();
  snoptProblemB(const char*name);
  ~snoptProblemB();

  int solve(int starttype, int m, int n, int ne, int negCon,
	    int nnCon, int nnObj, int nnJac,
	    int iObj, double ObjAdd,
	    snConB funcon, snObjB funobj,
	    double *Jval, int *indJ, int *locJ,
	    double *bl, double *bu, int *hs,
	    double *x, double *pi, double *rc,
	    int &nS, int &nInf, double &sInf, double &objective);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class sqoptProblem : public snoptProblem {
private:
  isqLog  sqLog;
  void init2zero();

public:
  sqoptProblem();
  sqoptProblem(const char*name);
  ~sqoptProblem();

  void initialize  (const char *prtfile, int summOn);
  int  setSpecsFile(const char *specname);
  void setLog      (isqLog sqLog);
  void setWorkspace(int m, int n, int neA, int ncObj, int nnH);

  int solve(int starttype, sqFunHx qpHx,
	    int m, int n, int neA, int ncObj, int nnH,
	    int iObj, double ObjAdd,
	    double *A, int *indA, int *locA,
	    double *bl, double *bu, double *cObj,
	    int *eType, int *hs, double *x, double *pi, double *rc,
	    int &nS, int &nInf, double &sInf, double &objective);
};

#endif /* SNOPTPROBLEM_H */
