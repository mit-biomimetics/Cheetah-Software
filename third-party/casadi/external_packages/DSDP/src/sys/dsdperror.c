#include "dsdpsys.h"
/*!
\file dsdperror.c
\brief Error codes returned for most subroutines
*/
#include <stdarg.h>
#include <sys/types.h>


#define DSDP_NULL 0
#define DSDP_MAX_PATH_LEN 200

static FILE      *DSDPLogInfoFile      = DSDP_NULL;
static FILE *dsdp_history=0;

typedef void* DSDPObject;

int DSDPFError(void *vobj, const char functionname[], int linen, const char filename[], const char message[], ...)
{
  va_list     Argp;
  int         urank;
  size_t      len;
  char        string[8*1024];

  DSDPFunctionBegin;
  DSDPLogInfoFile = stdout; 
  /*
 if (DSDPLogPrintInfo == DSDP_FALSE) DSDPFunctionReturn(0);
  if ((DSDPLogPrintInfoNull == DSDP_FALSE) && !vobj) DSDPFunctionReturn(0);
  */
  urank = 0;

  va_start(Argp, message);
  sprintf(string, "[%d] DSDP: %s(): Line %d in file %s ", urank,functionname,linen,filename);
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
  DSDPFunctionReturn(0);
}


void DSDPError(const char * functionname, int linen, const char filename[]){
  DSDPErrorPrintf("DSDP Error in function: %s , line %d of file %s \n",functionname, linen,filename);
  return;
}
typedef struct{
  void *mem;
  char fname[20];
  size_t size;
  int freed;
} DSDPMemory;

#define DSDPMEMMAX 1
static DSDPMemory DSDPMemoryTable[DSDPMEMMAX];

static long int mmmem=0;
#undef __FUNCT__  
#define __FUNCT__ "DSDPMMalloc"
int DSDPMMalloc(const char* fname, size_t size, void** mmem){
  void*tt=0;
  DSDPFunctionBegin;
  if (size>0){
#ifdef DSDPMATLAB
    tt=(void*)mxMalloc(size);
#else
    tt=(void*)malloc(size);
#endif    
    if (tt==NULL){ 
      *mmem=0;
      DSDPSETERR3(100,"Memory Error in routine '%s'. Cannot allocate %d bytes, %d MB\n",fname,size,(int)(size/1000000)); 
    }
    memset(tt,0,size);
    *mmem=tt;

    if (mmmem<DSDPMEMMAX){
      DSDPMemoryTable[mmmem].size=size;
      DSDPMemoryTable[mmmem].freed=0;
      strncpy(DSDPMemoryTable[mmmem].fname,fname,19);
      DSDPMemoryTable[mmmem].mem=tt;
    }
    mmmem++;

  } else {
    *mmem=0;
  }
  DSDPFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DSDPFFree"
int DSDPFFree(void** mmem){
  int j,gotit=0;
  DSDPFunctionBegin;
  if (mmem && *mmem){
    for (j=0;j<DSDPMEMMAX;j++){
      if (*mmem==DSDPMemoryTable[j].mem){
	DSDPMemoryTable[j].freed=1;
	gotit=1;
     }
    }  
    mmmem--;
#ifdef DSDPMATLAB
    mxFree(*mmem);
#else
    free(*mmem);
#endif    
    *mmem=0;
    if (0==1 && gotit==0){
      printf(" DSDP MEMORY Error: Already Freed? \n");
    }
  }
  DSDPFunctionReturn(0);
}

void DSDPMemoryLog(){
  int j;
  DSDPFunctionBegin;
  for (j=0;j<DSDPMEMMAX;j++){
    if (DSDPMemoryTable[j].size>0 && DSDPMemoryTable[j].freed==0){
      printf("%d, MEMORY Not FREED: %s, %d \n",j,DSDPMemoryTable[j].fname,(int)DSDPMemoryTable[j].size);
    }
  }
  
  DSDPLogInfo(0,2,"  MEMORY MALLOC NOT FREED: %ld\n",mmmem);

  return;
}
