#include "numchol.h"

static int ChkXlist(xlist *xt,
                    int   ind)
{
  return (xt->port[ind]==xt->idep);
} /* ChkXlist */

static void XtClear(xlist *xt)
{
  int i,sze;
  
  sze     =xt->last;
  xt->idep=xt->most+1;
  xt->lowp=xt->idep;
  xt->cure=sze;
  xt->ntot=0;
  
  for (i=0; i<xt->idep; i++)
    xt->head[i]=xt->last;
    
  for (i=0; i<sze; i++) {
    xt->port[i]=xt->idep;
    xt->fwrd[i]=xt->last;
    xt->bwrd[i]=xt->last;
  }
} /* XtClear */

int XtAlloc(int  last,
               int  most,
               char *info,
	       xlist**rr)
{
  xlist *r;
  int ierr=0;

  r=(xlist*) calloc(1,sizeof(xlist));
  if (!r) ExitProc(OutOfSpc,info);

  
  r->loca=TRUE;
  r->last=last;
  r->most=most;
  r->ntot=0;
    
  ierr=iAlloc(most+1,info,&r->head); if(ierr) return 1;
  ierr=iAlloc(last,info,&r->port); if(ierr) return 1;
  ierr=iAlloc(last,info,&r->fwrd); if(ierr) return 1;
  ierr=iAlloc(last,info,&r->bwrd); if(ierr) return 1;
  
  XtClear(r);
  
  *rr=r;
  return (0);
} /* XtAlloc */

void XtFree(xlist **xt)
{
  xlist *r=*xt;
  
  if (r) {
    if (r->loca) {
      iFree(&r->head);
      iFree(&r->port);
      iFree(&r->fwrd);
      iFree(&r->bwrd);
    }
    free(r);
    *xt=NULL;
  }
} /* XtFree */

int XtSucc(xlist *xt)
{
  int t,last=xt->last,most=xt->most,
      *head;

  if (xt->cure==last)
    return (FALSE);
    
  if (xt->fwrd[xt->cure]!=last)
    xt->cure=xt->fwrd[xt->cure];
      
  else {
    head=xt->head;
     
    for(t=xt->port[xt->cure]+1; t<=most && head[t]==last; ++t);
    
    if (t>most)
      xt->cure=last;
    else
      xt->cure=xt->head[t];
  }

  return (TRUE);
} /* XtSucc */

void XtDel(xlist *xt,
           int    e)
{
  int p;

  if (!ChkXlist(xt,e)) {
    
    if (xt->ntot<=0)
      ExitProc(SysError,NULL);
     
    xt->ntot--;
    
    if (xt->cure==e) {
      if (xt->ntot)
        XtSucc(xt);
      else
        xt->cure=xt->last;
    }
    
    p          =xt->port[e];
    xt->port[e]=xt->idep;
    
    if (xt->bwrd[e]!=xt->last)
      xt->fwrd[xt->bwrd[e]]=xt->fwrd[e];
    else
      xt->head[p]=xt->fwrd[e];
    
    if (xt->fwrd[e]!=xt->last)
      xt->bwrd[xt->fwrd[e]]=xt->bwrd[e];

    if (xt->head[p]==xt->last &&
         xt->lowp==p) {
      xt->lowp=xt->idep;
      if (xt->ntot) {
        for(++p; p<=xt->most; ++p){
          if (xt->head[p]!=xt->last){
            xt->lowp=p;
            break;
          }
        }
      }
    }
  }
} /* XtDel */

void XtPut(xlist *xt,
           int   e,
           int   p)
{
  if (0<=e && e<xt->last && 0<=p && p<=xt->most) {
    XtDel(xt,e);
    xt->ntot++;
    xt->port[e] =p;
    xt->fwrd[e] =xt->head[p];
    xt->bwrd[e]=xt->last;
    
    if (xt->head[p]!=xt->last)
      xt->bwrd[xt->head[p]]=e;
    
    xt->head[p]=e;
    xt->lowp    =min(p,xt->lowp);
  }
  
  else
    ExitProc(SysError,NULL);
} /* XtPut */

int XtLeast(xlist *xt)
{
  if (xt->lowp==xt->idep) {
    if (xt->ntot!=0)
      ExitProc(SysError,NULL);
    
    xt->cure=xt->last;
    return FALSE;
  }
  
  else {
    if (xt->ntot<=0)
      ExitProc(SysError,NULL);
    
    xt->cure=xt->head[xt->lowp];
    return TRUE;
  }
} /* XtLeast */

int XtGet(xlist *xt,
          int   *e,
          int   *p)
{
  if (xt->cure>xt->last)
    ExitProc(SysError,NULL);
  
  if (xt->cure==xt->last)
    return FALSE;
    
  *e=xt->cure;
  *p=xt->port[*e];
  
  return TRUE;
} /* XtGet */
