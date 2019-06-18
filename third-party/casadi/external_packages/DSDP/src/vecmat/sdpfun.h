
/*
 * common functions
 */
#if !defined (min)
#define min(a,b) ((a <= b)? (a) : (b))
#endif
#if !defined (max)
#define max(a,b) ((a >= b)? (a) : (b))
#endif
#if !defined (sign)
#define sign(a) ((a<0)? (-1) : (1))
#endif

/*
 * functions in sdpallo.c
 */
int       iAlloc(int,char*,int**);
void      iFree(int**);
int       dAlloc(int,char*,double**);
void      dFree(double**);
void      cFree(char**);
int       LvalAlloc(chfac*,char*);
int       CfcAlloc(int,char*,chfac**);
void      CfcFree(chfac**);
int       dPtAlloc(int,char*,double ***); 
void      dPtFree(double***);

/* 
 * functions in sdpmain.c 
 */ 
void    GetUhat(chfac*,double*,double*);

/*
 * functions in sdpshut.c
 */
void    ShutDown(void);
int     ExitProc(int,char*);


/*
 * functions in sdplib.a
 */
void    DotProd(double*,chfac*);

void    ForwSubst(chfac*,double*,double*);
int     iSum(int,int*);
void    dCopy(int,double*,double*);
int     ChlFact(chfac*,int*,double*,int);
void    ChlSolve(chfac*,double*,double*);
void    ChlSolveForward(chfac*, double*, double*);
void    ChlSolveBackward(chfac*, double*, double*);
int     SymbProc(int*,int*,int,chfac**);
void    iCopy(int,int*,int*);

int    OdAlloc(int,int,char*,order**);
void   OdFree(order**);
void   OdInit(order*,int*);
void   OdIndex(order*,int,int);
void   OdProc(order*,xlist*,int*,int*,int*,int*,int*,
              int*,int*,int*,int*,int*,int*,int*,int*);
int    GetOrder(order*,int*);

void   DotProd(double*,chfac*);

/* void   CfcInit(chfac*,symat*,double*); */
int    ChlFact(chfac*,int*,double*,int);
void   copyChl(chfac *, chfac *);

/* void   PermSmatx(smatx*,int*,int*); */

int    XtAlloc(int,int,char*,xlist**);
void   XtFree(xlist**);
int    XtSucc(xlist*);
void   XtDel(xlist*,int);
void   XtPut(xlist*,int,int);
int    XtLeast(xlist*);
int    XtGet(xlist*,int*,int*);

int    IptAlloc(int,int,int**,char*);
void   IptFree(int,int**);
int    LocIntPos(int,int,int*);
void   PermTransSym(int,int*,int*,int*,int*,int,int*,int*,int*);

void   iZero(int,int*,int*);
void   iFill(int,int,int*,int*,int*);
void   iSwap(int,int,int*);
void   iCopy(int,int*,int*);
int    iSum(int,int*);
void   dZero(int,double*,int*,int*);
void   dCopy(int,double*,double*);
void   dCat(int,int*,double*,double*);
double dSum(int,double*);
void   PlusByOne(int,int*,int*,int*);

void iSet(int, int, int *, int *);
void plusXs(int, int*, int*);

int MatMult4(chfac *,double *,double*,int);
int Mat4LogDet(chfac *,double *);
int MatZeroEntries4(chfac *);
int MatSolve4(chfac *,double *,double*,int);
void  CfcFree(chfac**);
int MatSetColumn4(chfac *, double *, int);
int MatAddColumn4(chfac *, double,double *, int);
int MchlSetup2(int m, chfac** A);
int MatSetValue4(chfac *, int,int,double, int);
int Mat4GetDiagonal(chfac*, double *,int);
int Mat4SetDiagonal(chfac*, double *,int);
int Mat4AddDiagonal(chfac*, double *,int);
int MatAddDiagonalElement(chfac*,int, double);
int Mat4View(chfac *);
int Mat4DiagonalShift(chfac*, double);
