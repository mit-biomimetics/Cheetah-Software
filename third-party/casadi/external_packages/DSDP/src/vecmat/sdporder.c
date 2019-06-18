#include "numchol.h"

int OdAlloc(int  nnod,
               int  nn0,
               char *info,
	       order **rr)
{
  order *r;
  int ierr=0;

  r=(order*)calloc(1,sizeof(order));
  if (!r) ExitProc(OutOfSpc,info);
  
  r->nnod=nnod;
  r->nn0 =nn0;
  
  ierr=iAlloc(nn0,info,&r->adjn); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->rbeg); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->rexs); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->rlen); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->rend); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->pres); if(ierr) return 1;
  ierr=iAlloc(nnod,info,&r->succ); if(ierr) return 1;
  *rr=r;
  return (0);
} /* OdAlloc */

void OdFree(order **od)
{
  order *r;
  
  if (*od) {
    r=*od;
    iFree(&r->adjn);
    iFree(&r->rbeg);
    iFree(&r->rexs);
    iFree(&r->rlen);
    iFree(&r->rend);
    iFree(&r->pres);
    iFree(&r->succ);
    free(*od);
    *od=NULL;  
  }
} /* OdFree */

void OdInit(order *od,
            int   *nnzi)
{
  int i,n=od->nnod;
  
  if (n) {
    od->rexs[0]=nnzi[0];
    od->rlen[0]=nnzi[0];
    od->rbeg[0]=0;
    od->pres[0]=n;
    od->succ[0]=1;
    for(i=1; i<od->nnod; ++i) {
      od->pres[i]=i-1;
      od->succ[i]=i+1;
      od->rexs[i]=nnzi[i];
      od->rlen[i]=nnzi[i];
      od->rbeg[i]=od->rbeg[i-1]+od->rlen[i-1];
    }

    od->succ[n-1]=n;
    od->last     =n-1;

    od->raft=od->rbeg[n-1]+nnzi[n-1];

    if (od->raft>od->nn0)
      ExitProc(OutOfSpc,"InitMmd");
  }  
} /* OdInit */

static void OdSet(order *od,
                  int   allow_eli,
                  xlist *elist,
                  int   *node_status,
                  int   *marker,
                  int   *isize,
                  int   *ilink,
                  int   *oinfo,
                  int   *osize,
                  int   *e,
                  int   *p)
{
  int i,n,deg,*rbeg,*rexs,*rend;

  n   =od->nnod;
  rbeg=od->rbeg;
  rexs=od->rexs;
  rend=od->rend;
  *e  =0;
  
  for (i=0; i<n; i++) {
    isize[i] =0;
    ilink[i] =n;
    osize[i] =0;
	oinfo[i] =n;
    marker[i]=0;
  }
  
  for(i=0; i<n; ++i) {
    rbeg[i] -= rexs[i];
    rend[i]  = 0;
  }
   
  for(i=0; i<n; ++i) {
    deg = rexs[i];
    if (!allow_eli||deg) {
      node_status[i]=1;
      XtPut(elist,i,deg);
    }
    else {
      node_status[i] = 0;
      marker[i]      = TRUE;
      p[*e]          = i;
      (*e)++;
    }
  }
} /* OdSet */

void OdIndex(order *od,
             int   i,
             int   j)
{
  if (i!=j) {
    od->adjn[od->rbeg[i]++]=j;
    od->adjn[od->rbeg[j]++]=i;
  }
} /* OdIndex */

static void OdArriv(order *od,
                    int   *node_status,
                    int   *marker,
                    int   *isize,
                    int   x,
                    int   *xdeg,
                    int   *rsze,
                    int   *esze,
                    int   *rchset)
{
  int *visited,i,n,y,z,l,s,t,f,stopt,
      stops,*adjn,*rbeg,*rexs,*rend;

  n    =od->nnod;
  adjn =od->adjn;
  rbeg =od->rbeg;
  rexs =od->rexs;
  rend =od->rend;
  *rsze=0;
  *esze=0;
  
  if (rexs[x]) {
    l=n;
    
    visited=marker;
    visited[x]=TRUE;
    
    for(t=rbeg[x], stopt=rbeg[x]+rend[x]; t<stopt; ++t) {
      y=adjn[t];
      
      if (node_status[y]!=0) {
        l--;
        rchset[l]=y;
        visited[y]=TRUE;
        
        for(s=rbeg[y], stops=rbeg[y]+rexs[y]; s<stops; ++s) {
          z=adjn[s];
          
          if (node_status[z]!=0) {
            if (!visited[z]) {
              visited[z]=TRUE;
               
              rchset[*rsze]=z;
              (*rsze)++;
            }
          }
        }
      }
    }
    
    f=rbeg[x]+rend[x];
    for(t=f, stopt=rbeg[x]+rexs[x]; t<stopt; ++t) {
      y=adjn[t];
      if (!visited[y]) {
        adjn[f++]=y;
        visited[y]=TRUE;
        
        rchset[*rsze]=y;
        (*rsze)++;
      }
    }
    
    rexs[x]=f-rbeg[x];
        
    *esze=n-l;
    visited[x] = FALSE;
    iZero(*rsze,visited,rchset);
    iZero(n-l,visited,rchset+l);
  }

  if (xdeg) {
    *xdeg = *rsze+isize[x];
    for(i=0; i<*rsze; ++i)
      *xdeg+=isize[rchset[i]];
  }
} /* OdArriv */

static void OdRenew(order *od,
                    int   *ilink,
                    int   x,
                    int   xdeg,
                    int   *e,
                    int   *p)
{
  int c,n;

  n        =od->nnod;
  od->ntot+=xdeg--;
  p[*e]    =x;
  (*e)++;
  for(c=x; ilink[c]!=n; c=ilink[c]) {
    od->ntot+=xdeg--;
    p[*e]    =ilink[c];
    (*e)++;
  }
} /* OdRenew */

static void OdCheck(order *od,
                    int   *node_status)
{
  int f,i,t,stopt,rnew,z,previous,n,*adjn,
      *rbeg,*rexs,*rlen,*rend,*pres,*succ;
    
  n   =od->nnod;
  adjn=od->adjn;
  rbeg=od->rbeg;
  rexs=od->rexs;
  rlen=od->rlen;
  rend=od->rend;
  pres=od->pres;
  succ=od->succ;
  
  f=0;
  previous=n;
  for(i=od->head; i!=n; i=succ[i]) {
    if (node_status[i]!=0) {
      rnew=f;
      
      for(t=rbeg[i], stopt=rbeg[i]+rend[i]; t<stopt; ++t) {
        z=adjn[t];
        if (node_status[z]==3)
          adjn[f++]=z;
      }
      
      rend[i]=f-rnew;
      
      for(stopt=rbeg[i]+rexs[i]; t<stopt; ++t) {
        z=adjn[t];
        if (node_status[z]!=0)
          adjn[f++]=z;
      }
      
      rexs[i]=f-rnew;
      rlen[i]=rexs[i];
      
      rbeg[i]=rnew;
      
      if (previous==n) {
        od->head=i;
        pres[i]=n;
      }
      
      else {
        succ[previous]=i;
        pres[i]=previous;
      }
      previous=i;
    }
  }
  
  if (previous!=n) {
    succ[previous]=n;
    od->raft=rbeg[previous]+rexs[previous];
  }
  
  od->last=previous;
} /* OdCheck */

static void OdAdd(order *od,
                  int   *node_status,
                  int   x,
                  int   newsze)
{
  int n,*adjn,*rbeg,*rexs,*rlen,*pres,*succ;

  n   =od->nnod;
  adjn=od->adjn;
  rbeg=od->rbeg;
  rexs=od->rexs;
  rlen=od->rlen;
  pres=od->pres;
  succ=od->succ;
      
  if (newsze<=rlen[x])
    return;

  if (od->raft+newsze>od->nn0)
    OdCheck(od,node_status);
    
  if (od->raft+newsze>od->nn0)
    ExitProc(OutOfSpc,"OdAdd");
    
  if (pres[x]!=n)
    rlen[pres[x]]+=rlen[x];
    
  iCopy(rexs[x],adjn+rbeg[x],adjn+od->raft);
  rbeg[x]=od->raft;
  rlen[x]=newsze;
  od->raft+=newsze;
  
  if (pres[x]==n) {
    if (succ[x]==n)
      od->head=x;
    else
      od->head=succ[x];
  }
  
  else {
    if (succ[x]==n)
      succ[pres[x]]=x;
    else
      succ[pres[x]]=succ[x];
  }
  
  if (succ[x]!=n)
    pres[succ[x]]=pres[x];
   
  if (od->last!=x) {
    succ[od->last]=x;
    pres[x]=od->last;
  }
  
  succ[x] =n;  
  od->last=x;
} /* OdAdd */

static int OdComb(order *od,
                  int   *node_status,
                  int   *marker,
                  int   *isize,
                  int   *ilink,
                  int   *osize,
                  int   xsize,
                  int   *xset)
{
  int i,n,rnew,rlen,x,icur;
  
  n   =od->nnod;
  rlen=0;

  if (xsize==0)
    rnew=n;
  else if (xsize==1)
    rnew=xset[0];
  else {
    rnew=xset[0];
    for(i=1; i<xsize; ++i)
      rlen+=1+isize[xset[i]];
    
    node_status[rnew]=1;
    osize[rnew]=0;
     
    for(icur=rnew; ilink[icur]!=n; icur=ilink[icur]);
    isize[rnew]+=rlen;
    
    for(i=1; i<xsize; ++i) {
      x=xset[i];
      
      node_status[x]=0;
      marker[x]=TRUE;
      
      ilink[icur]=x;
      
      for(icur=x; ilink[icur]!=n; icur=ilink[icur]);
      
      isize[x]=0;
    }
  }
  
  return (rnew);
} /* OdComb */

static int OdSelect(order *od,
                    xlist *elist,
                    int   *node_status,
                    int   *marker,
                    int   *isize,
                    int   *ilink,
                    int   *oinfo,
                    int   *osize,
                    int   x,
                    int   *rsze,
                    int   *rchset,
                    int   *ibuf1,
                    int   *ibuf2,
                    int   *mask2,
                    int   *e,
                    int   *p)
{
  int absorp,old,i,j,n,esze,y,z,l,f,t,stopt,s,
      o,stops,indsze,xdeg,e0,ssze,*slist,tsze,
      *tlist,sze,*adjn,*rbeg,*rexs,*rlen,*rend;
  
  adjn =od->adjn;
  rbeg =od->rbeg;
  rexs =od->rexs;
  rlen =od->rlen;
  rend =od->rend;
  n    =od->nnod;
  slist=ibuf1;
  
  e0 = *e;
  OdArriv(od,node_status,marker,isize,x,&xdeg,rsze,&esze,rchset);
  
  XtDel(elist,x);
  
  OdRenew(od,ilink,x,xdeg,e,p);
  
  for(i=n-esze; i<n; ++i) {
    node_status[rchset[i]]=0;
    marker[rchset[i]]=TRUE;
  }
  
  marker[x]=TRUE;
  iSet(*rsze,TRUE,marker,rchset);
  
  ssze=0;
  for(i=0; i<*rsze;) {
    y=rchset[i];
    
    if (node_status[y]==0||node_status[y]==3)
      ExitProc(SysError,NULL);
    
    f=rbeg[y];
    for(t=f, stopt=f+rend[y]; t<stopt; ++t) {
      z=adjn[t];
      if (node_status[z]==3) {
        adjn[f++]=z;
        
        if (!mask2[z]) {
          slist[ssze++]=z;
          mask2[z]=TRUE;
        }
      }
    }
    rend[y]=f-rbeg[y];
    
    for(stopt=rbeg[y]+rexs[y]; t<stopt; ++t) {
      z=adjn[t];
      if (!marker[z])
        adjn[f++]=z;
    }
    
    rexs[y]=f-rbeg[y];
     
    if (rexs[y]==0) {
      OdRenew(od,ilink,y,xdeg-(*e-e0),e,p);
      node_status[y] = 0;
      marker[y]      = TRUE;
      
      (*rsze)--;
      iSwap(i,*rsze,rchset);
    }
    
    else {
      if (rexs[y]>=rlen[y]) 
        ExitProc(SysError,NULL);
      
      if (rexs[y]>rend[y])
        adjn[rbeg[y]+rexs[y]]=adjn[rbeg[y]+rend[y]];
        
      rexs[y]++;
       
      adjn[rbeg[y]+rend[y]]=x;
      rend[y]++;
      
      i++;
    }
  }
  
  iSet(ssze,FALSE,mask2,slist);
  
  if (*rsze==0) {
    node_status[x]=0;
    marker[x]=TRUE;
  }
  
  else {
    node_status[x]=3;
     
    rend[x]=0;
    rexs[x]=0;
    if (*rsze>rlen[x])
      OdAdd(od,node_status,x,*rsze);
     
    rexs[x]=*rsze;
    iCopy(*rsze,rchset,adjn+rbeg[x]);
    
    tsze=0;
    tlist=ibuf2;
    for(i=0; i<ssze; ++i){
      y=slist[i];
      old=marker[y];
      marker[y]=TRUE;
      
      absorp=TRUE;
        
      indsze=n;
      l=n;
      
      f=rbeg[y];
      for(t=f, stopt=f+rexs[y]; t<stopt; ++t) {
        z=adjn[t];
        if (node_status[z]!=0) {
          adjn[f++]=z;
          
          if (marker[z]) {
            l--;
            slist[l]=z;
            
            if (!mask2[z]) {
              for(s=rbeg[z],stops=rbeg[z]+rexs[z];
                  s<stops &&marker[adjn[s]]; ++s);
              
              if (s==stops) {
                indsze--;
                iSwap(l,indsze,slist);
              }
              
              mask2[z]=TRUE;
              tlist[tsze++]=z;
            }
          }
          else
            absorp=FALSE;
        }
      }
      
      marker[y]=old;
      rexs[y]=f-rbeg[y];
      
      if (indsze<n) {
        z=OdComb(od,node_status,marker,
                 isize,ilink,osize,
                 n-indsze,slist+indsze);
        
        node_status[z]=1;
        
        sze=0;
        for(j=l; j<indsze; ++j) {
          o=slist[j];
          sze+=1+isize[o];
          node_status[o]=2;
          oinfo[o]=z;
        }
        osize[z]=max(osize[z],sze);
      }
      
      if (absorp) {
        node_status[y]=0;
        marker[y]=TRUE;
      }
    }
    
    iSet(tsze,FALSE,mask2,tlist);
  }
  
  marker[x]=(node_status[x]==0);
  
  for(t=0; t<*rsze; ++t) {
    z=rchset[t];
    marker[z]=(node_status[z]==0);
  }
  
  return (FALSE);
} /* OdSelect */

static int OdOrder(order *od,
                   int   *node_status,
                   int   *marker,
                   int   *isize,
                   int   x,
                   int   *ibuf1)
{
  int rsze,esze,deg;
  
  OdArriv(od,node_status,marker,isize,
          x,&deg,&rsze,&esze,ibuf1);
  
  return deg;
} /* OdOrder */

static void OdModf(order *od,
                         xlist *elist,
                         int   *node_status,
                         int   *marker,
                         int   *isize,
                         int   *oinfo,
                         int   rsze,
                         int   *rchset,
                         int   *ibuf1)
{

  int i,x,deg;
  
  for(i=0; i<rsze; ++i) {
    x=rchset[i];
    if (node_status[x]==2)
      if (node_status[oinfo[x]]==0||node_status[oinfo[x]]==3)
        node_status[x]=1;

    if (node_status[x]==1) {
      deg=OdOrder(od,node_status,marker,isize,x,ibuf1);
      XtPut(elist,x,deg-isize[x]);
    }
    
    else
      XtDel(elist,x);
  }
} /* mindeg_upddeg */

void OdProc(order *od,
            xlist *xt,
            int   *bbuf1,
            int   *bbuf2,
            int   *ibuf1,
            int   *ibuf2,
            int   *ibuf3,
            int   *ibuf4,
            int   *ibuf5,
            int   *ibuf6,
            int   *ibuf7,
            int   *ibuf8,
            int   *ibuf9,
            int   *intbuf1,
            int   *p)
{
  int    *mask2,total_fillin,use_mtmd=TRUE,*marker=bbuf1,
         *node_status,i,n,e,x,y,rsze,deg,mindeg,xsize,sze,
         j,*isize,*ilink,*oinfo,*osize,*rchset,*dependent,
         *slist;
  xlist *elist;

  elist=xt;

  isize    =ibuf1;
  ilink    =ibuf2;
  oinfo    =ibuf3;
  osize    =ibuf4;
  rchset   =ibuf5;  
  dependent=ibuf6;
  slist    =ibuf7;

  node_status=intbuf1;
  mask2=bbuf2;
  
  OdSet(od,TRUE,elist,node_status,marker,
        isize,ilink,oinfo,osize,&e,p);
  
  n=od->nnod;

  iSet(n,0,dependent,NULL);
  
  total_fillin=FALSE;
  for(; e<n && !total_fillin;) {
    
    XtLeast(elist);
    
    if (!XtGet(elist,&y,&mindeg)) {
      printf("\n No new nodes e=%d  n=%d",e,n);
      printf(" Node status: ");
      
      for(i=0; i<n; ++i)
        if (node_status[i]==1)
          printf("A\n");
        else if (node_status[i]==2)
          printf("\n O%d: rlen=%d oinfo=%d\n",
                 i,isize[i],oinfo[i]);
      
      ExitProc(SysError,NULL);
    }
    
    xsize=0;
    for(;use_mtmd;) {
      if (!XtGet(elist,&x,&deg)||deg>mindeg)
        break;
       
      if (node_status[x]!=1)
        XtDel(elist,x);
      
      else {
        XtSucc(elist);
        if (!dependent[x]) {
          
          total_fillin = OdSelect(od,elist,node_status,marker,
                                  isize,ilink,oinfo,osize,x,
                                  &rsze,rchset,ibuf8,ibuf9,mask2,
                                  &e,p);
          
          if (!total_fillin) {
            dependent[x]=2;
            slist[xsize++]=x;
            
            for(i=0; i<rsze; ++i) {
              y=rchset[i];
              if (!dependent[y]) {
                dependent[y]=1;
                slist[xsize++]=y;
              }
            }
          }
          
          if (!use_mtmd) break;
        }
      }
    }
    
    if (!total_fillin) {
      sze=0;
      for(j=0; j<xsize; ++j) {
        y=slist[j];
        if (dependent[y]==1 && node_status[y]!=0)
          slist[sze++]=y;
        dependent[y]=0;
      }
       
      OdModf(od,elist,node_status,marker,
             isize,oinfo,sze,slist,ibuf8);
        
    }
  }
  
  if (e<n) {
    sze=0;
    for(i=0; i<n; ++i)
      if (node_status[i]==2||node_status[i]==1)
        ibuf8[sze++]=i;
    
    x = OdComb(od,node_status,marker,isize,ilink,osize,
               sze,ibuf8);
     
    OdRenew(od,ilink,x,n-e-1,&e,p);
  }
  
} /* OdProc */

int GetOrder(order *od,
             int   *p)
{
  int   ierr=0,bbufs=2,*bbuf[2]={0},
        ibufs=9,*ibuf[9]={0},
        n,*iwmd;
  xlist *xt;

  n=od->nnod;

  ierr=XtAlloc(n,n+1,"xt, GetOrder",&xt); if(ierr) return FALSE;
  ierr=iAlloc(n,"ibuf21, GetOrder",&iwmd); if(ierr) return FALSE;
  
  IptAlloc(ibufs,n,ibuf,"ibuf, GetOrder");
  IptAlloc(bbufs,n,bbuf,"bbuf, GetOrder");
  
  OdProc(od,xt,ibuf[0],ibuf[1],ibuf[2],ibuf[3],ibuf[4],ibuf[5],
         ibuf[6],ibuf[7],ibuf[8],iwmd,bbuf[0],bbuf[1],p);
    

  /*
    XtFree(&xt);
  */

  free(xt->head);
  free(xt->port);
  free(xt->fwrd);
  free(xt->bwrd);
  free(xt);

  iFree(&iwmd);
  IptFree(ibufs,ibuf);
  IptFree(bbufs,bbuf); 
  return TRUE;
} /* GetOrder */
