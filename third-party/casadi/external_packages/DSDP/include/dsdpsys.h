#if !defined(__DSDP_KERNAL_H) 
#define __DSDP_KERNAL_H

/*!
\file dsdpsys.h
\brief Error handling, printing, and profiling
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

/* Define some macros for error checking */
#ifdef __FUNCT__
#undef __FUNCT__
#endif
#define __FUNCT__ "DSDPUnknownFunction"

/*
#ifdef __cplusplus
#define DSDPBEGINCROUTINES  extern "C" {
#define DSDPENDCROUTINES  }
#else
#define DSDPBEGINCROUTINES {
#define DSDPENDCROUTINES  }
#endif
*/

#ifdef __cplusplus
extern "C" {
#endif

extern void DSDPTime(double*);


extern int DSDPLogInfoAllow(int, char*);

extern void DSDPError(const char*, int, const char*);
extern void DSDPLogFInfo(void *vobj, int outlevel, const char message[], ...);
extern int DSDPFError(void *vobj, const char *func, int linen,const char *filef, const char message[], ...);

extern void DSDPMemoryLog(void);
extern int DSDPEventLogBegin(int);
extern int DSDPEventLogEnd(int);
extern int DSDPEventLogRegister(const char*, int*);
extern int DSDPEventLogInitialize(void);
extern int DSDPEventLogSummary(void);
extern int DSDPEventLogInitialize(void);

#ifdef __cplusplus
}
#endif


#ifndef DSDPCHKERR
#define DSDPCHKERR(a)  { if (a){ DSDPError(__FUNCT__,__LINE__,__FILE__); return a; } }
#endif

#ifdef DSDPFunctionReturn
#undef  DSDPFunctionReturn
#endif
#define DSDPFunctionReturn return

#ifdef DSDPFunctionBegin
#undef  DSDPFunctionBegin
#endif
#define DSDPFunctionBegin { }

#ifdef DSDPMATLAB
#include "mex.h"
#define DSDPPrintf   mexPrintf
#define DSDPErrorPrintf  mexPrintf
#endif

#include "dsdpmem.h"

#ifndef DSDPPrintf
#define DSDPPrintf printf
#endif

#ifndef DSDPErrorPrintf
#define DSDPErrorPrintf printf
#endif

#define DSDPLogInfo DSDPLogFInfo
/*#define DSDPLogInfo if(0)DSDPLogFInfo */



#define DSDPSETERR(a,b)         {DSDPFError(0,__FUNCT__,__LINE__,__FILE__,b); return (a); }
#define DSDPSETERR1(a,b,c)      {DSDPFError(0,__FUNCT__,__LINE__,__FILE__,b,c); return (a); }
#define DSDPSETERR2(a,b,c,d)    {DSDPFError(0,__FUNCT__,__LINE__,__FILE__,b,c,d); return (a); }
#define DSDPSETERR3(a,b,c,d,e)  {DSDPFError(0,__FUNCT__,__LINE__,__FILE__,b,c,d,e); return (a); }


/*
*/
#define DSDPMin(a,b) ((a <= b)? (a) : (b))
#define DSDPMax(a,b) ((a >= b)? (a) : (b))


#endif
