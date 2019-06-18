#ifndef SNOPT_H
#define SNOPT_H

/*
 * File:  snopt.h
 *   Header file for the SNOPT and SQOPT functions.
 */

#ifdef __cplusplus
extern "C" {
#endif

  typedef void (*isnLog)
  ( int *iAbort, int *info, int *HQNType, int KTcond[], int *MjrPrt, int *minimz,
    int *n, int *nb, int *nnCon0, int *nS,
    int *itn, int *nMajor, int *nMinor, int *nSwap,
    double *condHz, int *iObj, double *sclObj, double *ObjAdd,
    double *fMrt, double *PenNrm, double *step,
    double *prInf, double *duInf, double *vimax, double *virel, int hs[],
    int *ne, int nlocJ[], int locJ[], int indJ[], double Jcol[],
    double Ascale[], double bl[], double bu[], double fCon[], double yCon[], double x[],
    char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
    char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw );

  typedef void (*isnLog2)
  ( int *Prob, char *ProbTag, int *Elastc, int *gotR, int *jstFea, int *feasbl,
    int *m, int *mBS, int *nnH, int *nS, int *jSq, int *jBr, int *jSr,
    int *linesP, int *linesS, int *itn, int *itQP, int *kPrc, int *lvlInf,
    double *pivot, double *step, int *nInf, double *sInf, double *wtInf,
    double *ObjPrt, double *condHz, double *djqPrt, double *rgNorm,
    int kBS[], double xBS[],
    int iw[], int *leniw );

  typedef void (*isqLog)
  ( int *Prob, char *ProbTag, int *Elastc, int *gotR, int *jstFea, int *feasbl,
    int *m, int *mBS, int *nnH, int *nS, int *jSq, int *jBr, int *jSr,
    int *linesP, int *linesS, int *itn, int *itQP, int *kPrc, int *lvlInf,
    double *pivot, double *step, int *nInf, double *sInf, double *wtInf,
    double *ObjPrt, double *condHz, double *djqPrt, double *rgNorm,
    int kBS[], double xBS[],
    int iw[], int *leniw );

  typedef void (*isnSTOP)
  ( int *iAbort, int KTcond[], int *MjrPrt, int *minimz,
    int *m, int *maxS, int *n, int *nb,
    int *nnCon0, int *nnCon, int *nnObj0, int *nnObj, int *nS,
    int *itn, int *nMajor, int *nMinor, int *nSwap,
    double *condHz, int *iObj, double *sclObj, double *ObjAdd,
    double *fMrt, double *PenNrm, double *step,
    double *prInf, double *duInf, double *vimax, double *virel, int hs[],
    int *ne, int *nlocJ, int locJ[], int indJ[], double Jcol[], int *negCon,
    double Ascale[], double bl[], double bu[], double fCon[], double gCon[], double gObj[],
    double yCon[], double pi[], double rc[], double rg[], double x[],
    char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
    char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw );


  typedef void (*snFunA)
  ( int *Status, int *n,
    double x[],  int *needF, int *neF,  double F[],
    int *needG,  int *neG,  double G[],
    char   cu[], int *lencu,
    int    iu[], int *leniu,
    double ru[], int *lenru );

  typedef void (*snFunC)
  ( int *mode, int *nnObj, int *nnCon,
    int *nnJac, int *nnL, int *negCon, double x[],
    double *fObj,  double gObj[],
    double fCon[], double gCon[], int *Status,
    char   cu[], int *lencu,
    int    iu[], int *leniu,
    double ru[], int *lenru );

  typedef void (*snObjB)
  ( int    *mode, int   *nnObj,  double x[],
    double *fObj, double gObj[], int    *nState,
    char    cu[], int   *lencu,
    int     iu[], int   *leniu,
    double  ru[], int   *lenru );

  typedef void (*snConB)
  ( int     *mode, int    *nnCon,  int *nnJac,
    int   *negCon, double x[],
    double fCon[], double gCon[],  int *nState,
    char     cu[], int   *lencu,
    int      iu[], int   *leniu,
    double   ru[], int   *lenru );


  typedef void (*sqFunHx)
  ( int    *nnH, double x[], double Hx[], int *nState,
    char     cu[], int   *lencu,
    int      iu[], int   *leniu,
    double   ru[], int   *lenru );


  /* SNOPT */
  void f_sninit ( const char *name, int len, int summOn,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_snspec ( const char *specfile, int len, int *inform,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_sngetc ( const char *buffer, int lenb, char *ivalue,
		  int lenc, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sngeti ( const char *buffer, int len, int  *ivalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sngetr ( const char *buffer, int len, double *ivalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snset  ( const char *buffer, int len, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_snseti ( const char *buffer, int len, int iopt, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_snsetr ( const char *buffer, int len, double rvalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snsetprint ( const char *name, int len,
		      int iw[], int leniw, double rw[], int lenrw );
  void f_snend ( int iw[], int leniw, double rw[], int lenrw );


  /* SQOPT */
  void f_sqinit ( const char *name, int len, int summOn,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sqspec ( const char *specfile, int len, int *inform,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_sqgetc ( const char *buffer, int lenb, char *ivalue,
		  int lenc, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sqgeti ( const char *buffer, int len, int  *ivalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sqgetr ( const char *buffer, int len, double *ivalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_sqset  ( const char *buffer, int len, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sqseti ( const char *buffer, int len, int iopt, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );
  void f_sqsetr ( const char *buffer, int len, double rvalue, int *errors,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_sqsetprint ( const char *name, int len,
		      int iw[], int leniw, double rw[], int lenrw );
  void f_sqend ( int iw[], int leniw, double rw[], int lenrw );


  /* SNOPTA */
  void f_snopta ( int start, const char *name,
		  int nf, int n, double objadd, int objrow,
		  snFunA usrfun,
		  int iAfun[], int jAvar[], int neA, double A[],
		  int iGfun[], int jGvar[], int neG,
		  double xlow[], double xupp[],
		  double flow[], double fupp[],
		  double x[], int xstate[], double xmul[],
		  double f[], int fstate[], double fmul[],
		  int *inform, int *ns, int *ninf, double *sinf,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snkera ( int start, const char *name,
		  int nf, int n, double objadd, int objrow,
		  snFunA usrfun, isnLog snLog, isnLog2 snLog2,
		  isqLog sqLog, isnSTOP snSTOP,
		  int iAfun[], int jAvar[], int neA, double A[],
		  int iGfun[], int jGvar[], int neG,
		  double xlow[], double xupp[],
		  double flow[], double fupp[],
		  double x[], int xstate[], double xmul[],
		  double f[], int fstate[], double fmul[],
		  int *inform, int *ns, int *ninf, double *sinf,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snjac  ( int *info, int nf, int n, snFunA usrfun,
		  double x[], double xlow[], double xupp[],
		  int iAfun[], int jAvar[], int lenA, int *neA, double A[],
		  int iGfun[], int jGvar[], int lenG, int *neG,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snmema ( int *info, int nf, int n, int neA, int neG,
		  int *miniw, int *minrw,
		  int iw[], int leniw, double rw[], int lenrw );

  /* SNOPTB/C */
  void f_snoptb ( int start, const char *name,
		  int m, int n, int ne, int nncon,
		  int nnobj, int nnjac, int iobj, double objadd,
		  snConB funcon, snObjB funobj,
		  double jval[], int indj[], int locj[],
		  double bl[], double bu[], int hs[], double x[],
		  double pi[], double rc[], int *inform,
		  int *ns, int *ninf, double *sinf, double *obj,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snkerb ( int start, const char *name,
		  int m, int n, int ne, int nncon,
		  int nnobj, int nnjac, int iobj, double objadd,
		  snConB funcon, snObjB funobj,
		  isnLog snLog, isnLog2 snLog2,
		  isqLog sqLog, isnSTOP snSTOP,
		  double jval[], int indj[], int locj[],
		  double bl[], double bu[], int hs[], double x[],
		  double pi[], double rc[], int *inform,
		  int *ns, int *ninf, double *sinf, double *obj,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snoptc ( int start, const char *name,
		  int m, int n, int ne, int nncon,
		  int nnobj, int nnjac, int iobj, double objadd,
		  snFunC usrfun,
		  double jval[], int indj[], int locj[],
		  double bl[], double bu[], int hs[], double x[],
		  double pi[], double rc[], int *inform,
		  int *ns, int *ninf, double *sinf, double *obj,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snkerc ( int start, const char *name,
		  int m, int n, int ne, int nncon,
		  int nnobj, int nnjac, int iobj, double objadd,
		  snFunC usrfun,
		  isnLog snLog, isnLog2 snLog2,
		  isqLog sqLog, isnSTOP snSTOP,
		  double jval[], int indj[], int locj[],
		  double bl[], double bu[], int hs[], double x[],
		  double pi[], double rc[], int *inform,
		  int *ns, int *ninf, double *sinf, double *obj,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_snmem  ( int *info, int m, int n, int ne, int negcon,
		  int nncon, int nnobj, int nnjac,
		  int *miniw, int *minrw,
		  int iw[], int leniw, double rw[], int lenrw );

  /* SQOPT */
  void f_snkerq ( int start, sqFunHx qpHx, isqLog sqLog,
		  int m, int n, int neA, int ncObj, int nnH,
		  int iobj, double objadd, const char *name,
		  double A[], int indA[], int locA[],
		  double bl[], double bu[], double cObj[],
		  int eType[], int hs[], double x[],
		  double pi[], double rc[], int *inform,
		  int *ns, int *ninf, double *sinf, double *obj,
		  int *miniw, int *minrw,
		  int iu[], int leniu, double ru[], int lenru,
		  int iw[], int leniw, double rw[], int lenrw );

  void f_sqmem  ( int *info, int m, int n, int neA, int ncObj, int nnH,
		  int *miniw, int *minrw,
		  int iw[], int leniw, double rw[], int lenrw );
#ifdef __cplusplus
}
#endif

#endif /* SNOPT_H */
