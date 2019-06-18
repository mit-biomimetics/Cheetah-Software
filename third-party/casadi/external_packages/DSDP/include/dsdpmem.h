
#if !defined(__DSDP_MEM_H) 
#define __DSDP_MEM_H
/*! \file dsdpmem.h
    \brief Memory allocation in DSDP
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

extern int DSDPMMalloc(const char*, size_t, void**);
extern int DSDPFFree(void**);
/* Define some macros for memory management */

/*
#define DSDPFree(a)   0;printf("FREE: %s\n",__FUNCT__); free((*a))
#define DSDPMalloc(a,b)  0;printf("Malloc: %s\n",__FUNCT__); if (b){*(b)=malloc((a)); }
*/

#ifdef DSDPMATLAB

#define DSDPCALLOC1(VAR,TYPE,MERR) { \
   *(VAR) = (TYPE*)mxMalloc(sizeof(TYPE)); \
   *MERR=0; \
   if ( *(VAR)==0){*(MERR)=1;} \
   else {memset(*(VAR),0,sizeof(TYPE));} }

#define DSDPCALLOC2(VAR,TYPE,SIZE,MERR) { \
   *MERR=0; \
   *VAR=0; \
   if (SIZE>0){ \
     *(VAR) = (TYPE*)mxMalloc((SIZE)*sizeof(TYPE)); \
     if (*(VAR)==0){ *(MERR)=1;} \
     else {memset(*(VAR),0,(SIZE)*sizeof(TYPE));} \
   } \
}

#define DSDPFREE(VAR,MERR)   {if (*(VAR)){mxFree(*(VAR));}*(VAR)=0;*(MERR)=0;}
#endif

/*
#ifndef DSDPCALLOC1
#define DSDPCALLOC1(VAR,TYPE,MERR)  { (*MERR)=DSDPMMalloc(__FUNCT__,sizeof(TYPE),(void**)(VAR));}
#endif

#ifndef DSDPCALLOC2
#define DSDPCALLOC2(VAR,TYPE,SIZE,MERR)  { (*MERR)=DSDPMMalloc(__FUNCT__,(SIZE)*sizeof(TYPE),(void**)(VAR));}
#endif

#ifndef DSDPFREE
#define DSDPFREE(a,b)   {*(b)=DSDPFFree((void**)(a));}
#endif
*/

#ifndef DSDPCALLOC1
#define DSDPCALLOC1(VAR,TYPE,MERR) { \
   *(VAR) = (TYPE*)calloc(1, sizeof(TYPE)); \
   *MERR=0; \
   if ( *(VAR)==0){*(MERR)=1;} \
   else {    memset(*(VAR),0,sizeof(TYPE)); } \
}
#endif

#ifndef DSDPCALLOC2
#define DSDPCALLOC2(VAR,TYPE,SIZE,MERR) { \
   *MERR=0; \
   *VAR=0; \
   if (SIZE>0){ \
     *(VAR) = (TYPE*)calloc(SIZE, sizeof(TYPE)); \
     if (*(VAR)==0){ *(MERR)=1;} \
     else { memset(*(VAR),0,(SIZE)*sizeof(TYPE)); } \
   } \
}
#endif

#ifndef DSDPFREE
#define DSDPFREE(VAR,MERR)   {if (*(VAR)){free(*(VAR));}*(VAR)=0;*(MERR)=0;}
#endif

/*
#ifndef DSDPCALLOC1
#define DSDPCALLOC1(VAR,TYPE,SIZE,MERR) {*(VAR) = new TYPE;*MERR=0;}
#endif
#ifndef DSDPCALLOC2
#define DSDPCALLOC2(VAR,TYPE,SIZE,MERR) {*(VAR) = new TYPE[SIZE];*MERR=0;}
#endif
#ifndef DSDPFREE
#define DSDPFREE(a,b)   {delete(*(a));*(a)=0;*(b)=0;}
#endif
*/


#endif
