#include "numchol.h"

int IptAlloc(int  m,
              int  n,
              int  *ipt[],
              char *info)
{
  int i;
  
  if (n) {
    for (i=0; i<m; i++) {
      ipt[i]=(int*)calloc(n,sizeof(int));
      if (!ipt[i]){ 
	ExitProc(OutOfSpc,info);
	return 1;
      }
    }
  }
  return 0;
} /* AllocIptr */

void IptFree(int m,
             int *ipt[])
{
  int i;
  
  for (i=0; i<m; i++)
    iFree(&ipt[i]);
} /* IptFree */

int LocIntPos(int n,
              int i,
              int *v)
{
  int j;
  
  for(j=0; j<n && i!=v[j]; ++j);
  return (j);
} /* LocIntPos */

static void InsertDplRow(int i,
                         int ifir,
                         int *ifirst,
                         int *ilink)
{
  int temp;
  
  temp=ifirst[ifir];
  ifirst[ifir]=i;
  ilink[i]=temp;
} /* InsertDplRow */

void PermTransSym(int nrow,
                  int *fir,
                  int *nnz,
                  int *sub,
                  int *p,
                  int rowwise,
                  int *firt,
                  int *nnzt,
                  int *subt)
{
  int i,j,s,t,stopt;

  iZero(nrow,nnzt,NULL);
  
  if (rowwise) {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j =p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          nnzt[max(i,j)]++;
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          nnzt[max(i,j)]++;
        }
      }
    }
  }
  
  else {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j =p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          nnzt[min(i,j)]++;
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          nnzt[min(i,j)]++;
        }
      }
    }
  }

  firt[0]=0;
  for(i=1; i<nrow; ++i) {
    firt[i]  =firt[i-1] + nnzt[i-1];
    nnzt[i-1]=0;
  }
  nnzt[nrow-1]=0;

  if (rowwise) {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j=p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          if (i>j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          if (i>j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
  }
  
  else {
    if (p) {
      for(s=0; s<nrow; ++s) {
        j=p[s];
        for(t=fir[s], stopt=t+nnz[s]; t<stopt; ++t) {
          i=p[sub[t]];
          if (i<j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
    else {
      for(j=0; j<nrow; ++j) {
        for(t=fir[j], stopt=t+nnz[j]; t<stopt; ++t) {
          i=sub[t];
          if (i<j) {
            subt[firt[i]+nnzt[i]]=j;
            nnzt[i]++;
          }
          else {
            subt[firt[j]+nnzt[j]]=i;
            nnzt[j]++;
          }
        }
      }
    }
  }
} /* PermTransSym */

static void LocDplRow(int nrow,
                      int ncol,
                      int fir[],
                      int sze[],
                      int *sub,
                      int map[],
                      int ifirst[],
                      int ilink[],
                      int ilist[],
                      int *icount,
                      int i1nrow[])
{
  int i,rnew,n,oisze,isze,s,count,
      temp,k,nexti,*cur=i1nrow;
  
  n      =nrow;
  *icount=0;
  
  for (i=0; i<nrow; i++) {
    cur[i]  =0;
    ilink[i]=n;
    ilist[i]=n;
  }
  iSet(ncol,n,ifirst,NULL);
  
  isze =0;
  count=0;
  oisze=isze;
  rnew  =n;
  for(i=0; i<nrow; ++i) {
    if (map)
      for(; cur[i]<sze[i] && !map[sub[fir[i]+cur[i]]]; ++cur[i]);
    
    if ( cur[i]<sze[i] ) {
      s=sub[fir[i]+cur[i]];
      if ( ifirst[s]==n )
        ilist[isze++]=s;
      
      InsertDplRow(i,s,ifirst,ilink);
      
      cur[i]++;
    }
    
    else {
      temp    =rnew;
      rnew     =i;
      ilink[i]=temp;
    }
  }
  
  for(k=oisze; k<isze; ++k) {
    temp            =ifirst[ilist[k]];
    ifirst[ilist[k]]=n;
    ilist[k]        =temp;
  }
  
  if (rnew!=n) {
    count++;
    ilist[nrow-count]=rnew;
  }
  
  while(isze) {
    isze--;
    oisze      =isze;
    
    i          =ilist[isze];
    ilist[isze]=n;
    
    if (i==n)
      exit(0);
    
    rnew=n;
    if (ilink[i]==n)
      rnew=i;
    else {
      for(; i!=n; i=nexti) {
        nexti   =ilink[i];
        ilink[i]=n;
        
        if ( map )
          for(; cur[i]<sze[i] && !map[sub[fir[i]+cur[i]]]; ++cur[i]);
          
        if (cur[i]<sze[i]) {
          s =sub[fir[i]+cur[i]];
          cur[i]++;
           
          if (ifirst[s]==n)
            ilist[isze++]=s;
           
          temp     =ifirst[s];
          ifirst[s]=i;
          ilink[i] =temp;
        }
        
        else {
          temp    =rnew;
          rnew     =i;
          ilink[i]=temp;
        }
      }
    }
    
    for(k=oisze; k<isze; ++k) {
      temp            =ifirst[ilist[k]];
      ifirst[ilist[k]]=n;
      ilist[k]        =temp;
    }
    
    if (rnew!=n) {
      count++;
      ilist[nrow-count]=rnew;
    }
  }
  
  *icount=count;
  for(k=0; k<count; ++k)
    ilist[k]=ilist[nrow-count+k];
} /* LocDplRow */

static int CompIntElem(const void *e1,
                       const void *e2)
{
  int *i1,*i2;

  i1=(int *) e1;
  i2=(int *) e2;

  if (*i1<*i2)
    return (-1);
  else if(*i1>*i2)
    return (1);
  return (0);
} /* CompIntElem */

static void iSort(int n,
                  int *x)
{
  qsort((void *)x,n,sizeof(int),CompIntElem);
} /* iSort */

static void DetectDenseNodes(chfac *sf,
                             int   *i1nrow,
                             int   *i2nrow,
                             int   *i3nrow,
                             int   *i4nrow,
                             int   *i5nrow,
                             int   *i6nrow)
{
  int ierr=0,j,k,l,t,ndens,nil=sf->nrow,
      *subg=sf->subg,*ujbeg=sf->ujbeg,
      *ujsze=sf->ujsze,*usub=sf->usub,
      *fir,*sze,*ilist,*ilink;
  
  if (!sf->nsnds||
      !i1nrow  || !i2nrow || !i3nrow ||
      !i4nrow  || !i5nrow || !i6nrow) {
    sf->sdens=FALSE;
    return;
  }
  
  sf->sdens  =TRUE;
  fir        =i1nrow;
  sze        =i2nrow;
  ilist      =i3nrow;
  ilink      =i4nrow;
  
  sf->nsndn=0;
  
  l=subg[sf->nsnds-1];
  for(k=0; k+1<sf->nsnds; ++k) {
    j=subg[k];
    for(t=0; t<ujsze[j] && usub[ujbeg[j]+t]<l; ++t);
    
    fir[k] =ujbeg[j]+t;
    sze[k] =ujsze[j]-t;
  }
  
  LocDplRow(sf->nsnds-1,sf->nrow,fir,sze,usub,
            NULL,
            i6nrow,ilink,ilist,&ndens,i5nrow);
  
  ierr=iAlloc(ndens+1,  "sf->dhead, DetectDenseNodes",&sf->dhead);if(ierr) return;
  ierr=iAlloc(sf->nsnds,"sf->dsub, DetectDenseNodes",&sf->dsub);if(ierr) return;
  ierr=iAlloc(sf->nsnds,"sf->dbeg, DetectDenseNodes",&sf->dbeg);if(ierr) return;
  
  nil        =sf->nsnds-1;
  sf->ndens   =0;
  sf->nsndn   =0;
  sf->dhead[0]=0;
  
  for(k=0; k<ndens; ++k) {
    j=ilist[k];
    if ( sze[j] ) {
      sf->dhead[sf->ndens+1]=sf->dhead[sf->ndens];
      sf->ndens++;
      for(; j!=nil; j=ilink[j]) {
        sf->nsndn                    +=  sf->subg[j+1]-sf->subg[j];
        sf->dsub[sf->dhead[sf->ndens]]=j;
        sf->dbeg[sf->dhead[sf->ndens]]=fir[j]-ujbeg[subg[j]];
        sf->dhead[sf->ndens]++;
      }
      iSort(sf->dhead[sf->ndens]-sf->dhead[sf->ndens-1],
            sf->dsub+sf->dhead[sf->ndens-1]);
      
      for(t=sf->dhead[sf->ndens-1]; t<sf->dhead[sf->ndens]; ++t)
        sf->dbeg[t]=fir[sf->dsub[t]]-ujbeg[subg[sf->dsub[t]]];
    }
  }
} /* DetectDenseNodes */

static int ChlSymb(chfac *sf,
                   int   ulnnz)
{
  int ierr=0,chksn,i,j,t,stopt,sze,first,cur,k,
      ffree=0,ipos,nrow=sf->nrow,nil=nrow,
      *nnz,*fir,*pssub,*link,*buf,*mask,
      *usub,*tusub,*i1nrow,*i2nrow,*i3nrow,
      *i4nrow,*p=sf->perm,*invp=sf->invp,
      *ujbeg=sf->ujbeg,*ujsze=sf->ujsze,
      *subg=sf->subg;
  
  ierr=iAlloc(sf->snnz,"pssub, ChlSymb",&pssub);if(ierr) return FALSE;
  
  for(i=0; i<nrow; ++i)
    invp[p[i]]=i;
  
  nnz=sf->uhead;
  fir=sf->subg;
  
  PermTransSym(nrow,sf->shead,sf->ssize,sf->ssub,
               invp,TRUE,fir,nnz,pssub);
  
  PermTransSym(nrow,fir,nnz,pssub,NULL,FALSE,
               sf->shead,sf->ssize,sf->ssub);
  
  iFree(&pssub);
  
  k       =ulnnz+nrow;
  ierr =iAlloc(k,"usub, ChlSymb",&usub);if(ierr) return FALSE;
  buf     =usub+ulnnz;
  
  mask=sf->uhead;
  
  link=invp;
  
  iZero(nrow,mask,NULL);
  iSet(nrow,nil,link,NULL);
  
  ffree   =0;
  sf->nsnds=0;
  subg[0] =0;
  for(i=0; i<nrow; ++i) {
    sze  =sf->ssize[i];
    first=nil;
    cur  =link[i];
    chksn=FALSE;
    
    if (cur==nil) {
      
      subg[sf->nsnds+1] =1 + subg[sf->nsnds];
      ujsze[i]          =sze;
      ujbeg[i]         =ffree;
      ffree            += sze;
      
      iCopy(sze,sf->ssub+sf->shead[i],usub+ujbeg[i]);
      if (sze) {
        first=usub[ujbeg[i]];
        for(cur=first; link[cur]!=nil; cur=link[cur]);
        link[cur]     =sf->nsnds;
        link[sf->nsnds]=nil;
      }
      sf->nsnds++;
    }
    
    else {
      mask[i]=1;
      
      iCopy(sze,sf->ssub+sf->shead[i],buf);
      iSet(sze,1,mask,buf);
        
      for(; cur!=nil; cur=link[cur]) {
        chksn |= (1+cur==sf->nsnds);
        k     =subg[cur];
        
        for(t=ujbeg[k], stopt=t+ujsze[k]; t<stopt; ++t) {
          j=usub[t];
          if ( j>i && !mask[j] ) {
            buf[sze]=j;
            mask[j] =1;
            sze++;
          }
        }
      }
      
      if (chksn) {
        k    =subg[sf->nsnds-1];
        chksn=sze==( ujsze[k]-(subg[sf->nsnds]-subg[sf->nsnds-1]) );
      }
      
      first  =nrow;
      mask[i]=0;
      for(t=0; t<sze; ++t) {
        j     =buf[t];
        mask[j]=0;
        first  =min(j,first);
      }
      
      if (chksn) {
        ipos=LocIntPos(ujsze[i-1],i,usub+ujbeg[i-1]);
         
        if (ipos==ujsze[i-1]) 
          ExitProc(SysError,NULL);
        
        iSwap(ujbeg[i-1],ipos+ujbeg[i-1],usub);

        subg[sf->nsnds]++;
        ujbeg[i] =ujbeg[i-1]+1;
        ujsze[i]  =ujsze[i-1]-1;
        
        if (usub[ujbeg[i]-1]!=i)
          ExitProc(SysError,NULL);

        if ( first!=nil ) {
          for(cur=first; link[cur]!=nil; cur=link[cur]);
          link[cur]       =sf->nsnds-1;
          link[sf->nsnds-1]=nil;
        }
      }
      
      else {
        subg[sf->nsnds+1] =1 + subg[sf->nsnds];
        ujbeg[i]         =ffree;
        ujsze[i]          =sze;
        ffree            += sze;
        
        if (ffree>ulnnz)
          ExitProc(SysError,NULL);
        
        iCopy(sze,buf,usub+ujbeg[i]);
        
        if ( first!=nil ) {
          for(cur=first; link[cur]!=nil; cur=link[cur]);
          link[cur]     =sf->nsnds;
          link[sf->nsnds]=nil;
        }
        sf->nsnds++;
      }
    }
     
    if (ujsze[i]+1==nrow-i)
      break;
  }
  
  for(++i; i<nrow; ++i) {
    ujsze[i] =ujsze[i-1]-1;
    ujbeg[i]=ujbeg[i-1]+1;
     
    subg[sf->nsnds]++;
  }
   
  ierr=iAlloc(ffree,"tusub, ChlSymb",&tusub);if(ierr) return FALSE;
    
  fir=buf;
  nnz=sf->uhead;
  
  iZero(nrow,nnz,NULL);
  
  for(k=0; k<sf->nsnds; ++k) {
    j=subg[k];
    plusXs(ujsze[j],nnz,usub+ujbeg[j]);
  }
  
  fir[0]=0;
  for(k=1; k<nrow; ++k)
    fir[k]=fir[k-1] + nnz[k-1];
    
  iZero(nrow,nnz,NULL);
  
  for(k=0; k<sf->nsnds; ++k) {
    j=subg[k];
    for(t=ujbeg[j], stopt=t+ujsze[j]; t<stopt; ++t) {
      i                    =usub[t];
      tusub[fir[i]+nnz[i]] =j;
      nnz[i]++;
    }
    ujsze[j]=0;
  }
  
  for(i=0; i<nrow; ++i) {
    for(t=fir[i], stopt=t+nnz[i]; t<stopt; ++t) {
      j                      =tusub[t];
      usub[ujbeg[j]+ujsze[j]] =i;
      ujsze[j]++;
    }
  }
  
  iFree(&tusub);

  if (ffree<=sf->ujnz) {
    iCopy(ffree,usub,sf->usub);
    iFree(&usub);
  }
  
  else {
    sf->ujnz=0;
    iFree(&sf->usub);
    
    ierr=iAlloc(ffree,"sf->usub, ChlSymb",&sf->usub);if(ierr) return FALSE;
    iCopy(ffree,usub,sf->usub);
    
    sf->ujnz=ffree;
    iFree(&usub);
  }
  
  ierr=iAlloc(4*nrow,"i1nrow, ChlSymb",&i1nrow);if(ierr) return FALSE;
  i2nrow=i1nrow+nrow;
  i3nrow=i2nrow+nrow;
  i4nrow=i3nrow+nrow;
  
  DetectDenseNodes(sf,sf->uhead,sf->invp,
                   i1nrow,i2nrow,i3nrow,i4nrow);
  
  iFree(&i1nrow);
  
  sf->uhead[0]=0;
  for(i=1; i<nrow; ++i)
    sf->uhead[i]=sf->uhead[i-1]+sf->ujsze[i-1];
  
  for(i=0; i<nrow; ++i)
    invp[p[i]]=i;
  
  for(k=0; k<sf->nsnds; ++k)
    if ( subg[k]+1!=subg[k+1] )
      break;
      
  ierr=iAlloc(3*sf->n,NULL,&sf->iw);if(ierr) return FALSE;
  ierr=dAlloc(2*sf->n,NULL,&sf->rw);if(ierr) return FALSE;
  sf->factor=0;

  return TRUE;
} /* ChlSymb */

int SymbProc(int*    isze,
             int*    jsub,
             int  n,
	     chfac** sf)
{
  int    ierr,i,k,t,snnz,lnnz,*nnzi,nrow,resp;
  chfac* cf;
  order  *od;
  
  /*
   * set data for symmetric matrix to be factorized
   */
  ierr=CfcAlloc(n,"sdt->sf, SymbProc",&cf);if(ierr) return FALSE;

  nrow=cf->nrow;
  for (snnz=0,i=0; i<nrow; i++)
    snnz+=isze[i];
  /*
  if (!snnz)
    return TRUE;
  */
  ierr=iAlloc(snnz,"cf, SymbProc",&cf->ssub); if(ierr) return FALSE;
  cf->snnz=snnz;
  
  nnzi=cf->perm;
  iZero(nrow,nnzi,NULL);
  
  k=0;
  for (i=0; i<nrow; i++) {
    snnz        =isze[i];
    cf->shead[i]=k;
    cf->ssize[i]=snnz;
    k+=snnz;
  }
  iCopy(k,jsub,cf->ssub);
  
  nnzi=cf->perm;
  iZero(nrow,nnzi,NULL);
  
  for (i=0; i<nrow; i++) {
    nnzi[i]+=cf->ssize[i];
    plusXs(cf->ssize[i],nnzi,cf->ssub+cf->shead[i]);
  }
  
  ierr=OdAlloc(nrow,2*cf->snnz,"od, PspSymbo",&od); if(ierr) return FALSE;
  nnzi=cf->perm;
  
  OdInit(od,nnzi);
  for (i=0; i<nrow; i++)
    for (t=0; t<cf->ssize[i]; ++t)
      OdIndex(od,i,cf->ssub[cf->shead[i]+t]);
  
  GetOrder(od,cf->perm);
  
  lnnz=od->ntot;
  OdFree(&od);
  
  resp=ChlSymb(cf,lnnz);
  LvalAlloc(cf,"cf, PspSymb");
  /* sdt->sf=cf; */

  *sf=cf;
  return resp;
} /* SymbProc */

int MchlSetup2(int m, chfac** A)
{
  int   ierr,i,j,k,lnnz,mnnz;
  chfac *mf;
  
  ierr =CfcAlloc(m,NULL,&mf); if(ierr) return 1;
  *A=mf;

  mnnz=m*(m-1)/2;
  if (!mnnz && m>1)
    return TRUE;
  
  lnnz=mnnz;
  ierr=iAlloc(mnnz,NULL,&mf->ssub);if(ierr) return 1;

  mf->snnz=mnnz;
  k=0;
  for (i=0; i<m; i++)
  {
    mnnz         =m-(i+1);
    mf->shead[i] =k;
    mf->ssize[i] =mnnz;
    for (j=0; j<mnnz; j++)
      mf->ssub[k+j]=(j+1+i);
    k           +=mnnz;
    mf->perm[i]=i;
  }


  k=ChlSymb(mf,lnnz);

  iFree(&mf->ssub);
  iFree(&mf->shead);
  iFree(&mf->ssize);

  /* This part is different */
  mf->alldense=1;
  iFree(&mf->invp);
  mf->invp=mf->perm;

  iFree(&mf->ujbeg);
  mf->ujbeg=mf->perm;

  iFree(&mf->usub);
  mf->usub=mf->perm+1;

  ierr=LvalAlloc(mf,"cf, PspSymb");if(ierr) return 1;

  return 0;
  
} /* MchlSetup2 */
