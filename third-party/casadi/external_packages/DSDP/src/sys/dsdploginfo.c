/*$Id: ploginfo.c,v 1.22 2001/03/23 23:20:50 balay Exp $*/
/*
      DSDPLogInfo() is contained in a different file from the other profiling to 
   allow it to be replaced at link time by an alternative routine.
*/
#include <stdarg.h>
#include <sys/types.h>
#include <stdlib.h>
/*#include <malloc.h> */ /* Commented out since it is obsolete and causes compilation problems on OS X.*/
#include "dsdpsys.h"
#include "dsdpbasictypes.h"

#define DSDP_NULL 0
#define DSDP_MAX_PATH_LEN 200

/*!
\file dsdploginfo.c
\brief Print the progress of the DSDP solver.
*/

typedef void* DSDPObject;
/*
extern FILE *dsdp_history;
*/
static FILE *dsdp_history=0;
/*
  The next three variables determine which, if any, DSDPLogInfo() calls are used.
  If DSDPLogPrintInfo is zero, no info messages are printed. 
  If DSDPLogPrintInfoNull is zero, no info messages associated with a null object are printed.

  If DSDPLogInfoFlags[OBJECT_COOKIE - DSDP_COOKIE] is zero, no messages related
  to that object are printed. OBJECT_COOKIE is, for example, MAT_COOKIE.
*/
static int DSDPLogPrintInfo     = 0;
static int DSDPLogPrintInfoNull = 0;
static FILE      *DSDPLogInfoFile      = DSDP_NULL;
static int rank=0;
 
void DSDPSetRank(int rrank){
  rank=rrank;
  return;
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPLogInfoAllow"
/*@C
    DSDPLogInfoAllow - Causes DSDPLogInfo() messages to be printed to standard output.

    Not Collective, each processor may call this seperately, but printing is only
    turned on if the lowest processor number associated with the DSDPObject associated
    with the call to DSDPLogInfo() has called this routine.

    Input Parameter:
+   flag - DSDP_TRUE or DSDP_FALSE
-   filename - optional name of file to write output to (defaults to stdout)

    Options Database Key:
.   -log_info [optional filename] - Activates DSDPLogInfoAllow()

    Level: advanced

   Concepts: debugging^detailed runtime information
   Concepts: dumping detailed runtime information

.seealso: DSDPLogInfo()
@*/
int DSDPLogInfoAllow(int flag, char *filename)
{
  char fname[DSDP_MAX_PATH_LEN], tname[5];
  int  prank=0;
  char*  ierr;

  DSDPFunctionBegin;
  if (flag && filename) {
    sprintf(tname, ".%d", prank);
    ierr = strcat(fname, tname);
  } else if (flag) {
    DSDPLogInfoFile = stdout;
  }
  DSDPLogPrintInfo     = flag;
  DSDPLogPrintInfoNull = flag;
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPLogInfo"
/*@C
    DSDPLogInfo - Logs informative data, which is printed to standard output
    or a file when the option -log_info <file> is specified.

    Collective over DSDPObject argument

    Input Parameter:
+   vobj - object most closely associated with the logging statement
-   message - logging message, using standard "printf" format

    Options Database Key:
$    -log_info : activates printing of DSDPLogInfo() messages 

    Level: intermediate

    Example of Usage:
$
$     double alpha
$     DSDPLogInfo(0,"Matrix uses parameter alpha=%g\n",alpha);
$

   Concepts: runtime information

.seealso: DSDPLogInfoAllow()
@*/
void  DSDPLogFInfo(void *vobj, int outlevel, const char message[], ...)  
{
  va_list     Argp;
  int         urank;
  size_t      len;
  char        string[8*1024];

  DSDPFunctionBegin;
  DSDPLogInfoFile = stdout; 
  if (DSDPLogPrintInfo < outlevel) return;
  if ((DSDPLogPrintInfoNull < outlevel) && !vobj) return;

  urank = 0;
  if (rank>0) return;

  va_start(Argp, message);
  sprintf(string, "[%d][%2d] DSDP: ", rank,outlevel);
  len = strlen(string);
#if defined(DSDP_HAVE_VPRINTF_CHAR)
  vsprintf(string+len, message, (char *) Argp);
#else
  vsprintf(string+len, message, Argp);
#endif
  fprintf(DSDPLogInfoFile, "%s", string);
  fflush(DSDPLogInfoFile);
  if (dsdp_history) {
#if defined(DSDP_HAVE_VPRINTF_CHAR)
    vfprintf(dsdp_history, message, (char *) Argp);
#else
    vfprintf(dsdp_history, message, Argp);
#endif
  }
  va_end(Argp);
  return;
  /*  DSDPFunctionReturn(0); */
}
