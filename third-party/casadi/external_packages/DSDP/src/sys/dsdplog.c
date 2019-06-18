#include "dsdpsys.h"
#include <stdio.h>
#include <stdlib.h>

#define DSDP_MAX_EVENT_NAME_LENGTH 50
#define DSDP_MAX_EVENTS            30

/*! \file dsdplog.c
\brief Profile the performance of DSDP

*/
FILE *dsdpoutputfile;

typedef struct {
  int counter;
  double begintime;
  double totaltime;
  char ename[DSDP_MAX_EVENT_NAME_LENGTH];
} EventInfo;


typedef struct{
  EventInfo event[DSDP_MAX_EVENTS];
  int nevents;
  int neventsmax;
  double time0;
} EventLog;

static EventLog eventlog;


#undef __FUNCT__
#define __FUNCT__ "DSDPEventLogBegin"
int DSDPEventLogBegin(int eventid){
  double tt;
  DSDPTime(&tt);
  if (eventid<=0){return 0;}
  if (eventlog.event[eventid].begintime!=0 && eventid!=DSDP_MAX_EVENTS-1){
    DSDPPrintf("Timing error: id: %d %s.  Call begin without calling end.%4.4e\n",eventid,eventlog.event[eventid].ename,eventlog.event[eventid].begintime);
  };
  eventlog.event[eventid].begintime=tt;
  eventlog.event[eventid].counter++;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPEventLogEnd"
int DSDPEventLogEnd(int eventid){
  double tt;
  DSDPTime(&tt);
  if (eventid<=0){return 0;}
  tt=tt-eventlog.event[eventid].begintime;
  eventlog.event[eventid].totaltime+=tt;
  eventlog.event[eventid].begintime=0;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPEventLogRegister"
int DSDPEventLogRegister(const char *ename, int *eventid){
  int id;
  id=eventlog.nevents;
  if (id<0 || id>=DSDP_MAX_EVENTS){ *eventid=DSDP_MAX_EVENTS-1;return 0;}
  eventlog.event[id].begintime=0;
  eventlog.event[id].totaltime=0;
  eventlog.event[id].counter=0;
  strncpy(eventlog.event[id].ename,ename,DSDP_MAX_EVENT_NAME_LENGTH-1);
  eventlog.nevents++;
  *eventid=id;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPEventLogInitialize"
int DSDPEventLogInitialize(void){
  int i;
  double t0;
  DSDPTime(&t0);
  eventlog.time0=t0;
  for (i=0;i<DSDP_MAX_EVENTS;i++){
    eventlog.event[i].begintime=0;
    eventlog.event[i].totaltime=0;
    eventlog.event[i].counter=0;
    strncpy(eventlog.event[i].ename,"",DSDP_MAX_EVENT_NAME_LENGTH-1);
  }
  eventlog.nevents=1;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DSDPEventLogSummary"
int DSDPEventLogSummary(void){
  int i;
  double etime,ttime,tfinal;
  DSDPTime(&tfinal);
  if (tfinal==0){
    DSDPPrintf("DSDP Timing is not turned on.  Check installation and recompile. \n\n");
  }
  ttime=tfinal-eventlog.time0;
  /*  DSDPMemoryLog(); */
  DSDPPrintf("PERFORMANCE SUMMARY\n");
  DSDPPrintf("                     Event                      Calls    Time(s)   Time(%%)\n");
  DSDPPrintf("--------------------------------------------------------------------------\n");
  for (i=1;i<eventlog.nevents;i++){
    etime=eventlog.event[i].totaltime;
    if (etime==0 || eventlog.event[i].counter==0) continue;
    DSDPPrintf(" %40s   %9d   %4.4e  %5.2f\n",eventlog.event[i].ename,eventlog.event[i].counter,etime,100*etime/ttime);
  }
  DSDPPrintf("--------------------------------------------------------------------------\n");

  if (dsdpoutputfile){
    fprintf(dsdpoutputfile,"PERFORMANCE SUMMARY\n");
    fprintf(dsdpoutputfile,"                     Event                      Calls    Time(s)   Time(%%)\n");
    fprintf(dsdpoutputfile,"--------------------------------------------------------------------------\n");
    for (i=1;i<eventlog.nevents;i++){
      etime=eventlog.event[i].totaltime;
      if (etime==0 || eventlog.event[i].counter==0) continue;
      fprintf(dsdpoutputfile," %40s   %9d   %4.4e  %5.2f\n",eventlog.event[i].ename,eventlog.event[i].counter,etime,100*etime/ttime);
    }
    fprintf(dsdpoutputfile,"--------------------------------------------------------------------------\n");
  }
  fflush(NULL);
  return 0;
}

