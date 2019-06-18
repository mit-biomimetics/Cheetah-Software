/* DSDPTime could return 0 and still work */
/*!
\file dsdptime.c
\brief Timing routines for GNU and Microsoft compilers.
*/
/* 
#define DSDP_TIME 
*/

#include "dsdpsys.h"

#ifdef DSDP_MS_TIME
#include <ctype.h>
#include <time.h>
void DSDPTime(double * ttime) {  /* MICROSOFT COMPILER */
  clock_t t=clock();
  double tscale=0.001;
  (*ttime)=((double)t) * tscale;
}
#else
#ifdef DSDP_TIME
#include <sys/time.h> 
void DSDPTime(double * ttime) {  /* USED IN LINUX */
  static struct timeval _tp;
  *ttime=0;
  gettimeofday(&_tp,(struct timezone *)0);
  (*ttime)=((double)_tp.tv_sec)+(1.0e-6)*(_tp.tv_usec);
}
#else
void DSDPTime(double * ttime) { *ttime=0; return; } /* NO TIME */
#endif
#endif
/* for Microsoft */
/*
*/
