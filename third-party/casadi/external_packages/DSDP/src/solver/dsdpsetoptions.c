#include "dsdp5.h"
/*! \file dsdpsetoptions.c
   \brief Set DSDP options from file or command line arguments
*/

/*
 static char[] ZGAPTOL="-gaptol";
 static char[] ZPRINT="-print";
 static char[] ZPENALTY="-penalty";
 static char[] ZBIGM="-bigM";
 static char[] ZMAXIT="-maxit";
 static char[] ZR0="-r0";
 static char[] ZZBAR="-zbar";
 static char[] ZINFDTOL="-infdtol";
 static char[] ZINFPTOL="-infptol";
 static char[] ZRHO="-rho";
 static char[] ZDRHO="-drho";
 static char[] ZBOUNDY="-boundy";
 static char[] ZSTEPTOL="-steptol";
 static char[] ZREUSE="-reuse";
 static char[] ADADD="-dadd";
 static char[] ZDBOUND="-dbound";
 static char[] ZMU0="-mu0";
 static char[] DOBJMIN="-dobjmin";
*/

/*!
\fn int DSDPSetOptions(DSDP dsdp,char *runargs[], int nargs);
\brief Read command line arguments to set options in DSDP.
\param dsdp is the solver
\param runargs is the array of strings representing the options
\param nargs is the number of arguments
\sa DSDPReadOptions()
\sa DSDPPrintOptions()
\sa DSDPSetGapTolerance()
\sa DSDPSetStandardMonitor()
\sa DSDPSetPenaltyParameter()
\sa DSDPSetMaxIts()
\sa DSDPSetR0()
\sa DSDPSetYBounds()
\sa DSDPSetPotentialParameter()
\ingroup DSDPBasic
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPSetOptions"
int DSDPSetOptions(DSDP dsdp,char *runargs[],int nargs){

  int kk, info,reuse;
  int maxit,rpos,drho,iloginfo;
  double penalty,rho,zbar,cc,r0,mu0,gaptol,dbound,dd;
  double ylow,yhigh,maxtrust,steptol,inftol,infptol,pnormtol;

  DSDPFunctionBegin;

  for (kk=0; kk<nargs-1; kk++){
    if (strncmp(runargs[kk],"-gaptol",5)==0){
      gaptol=atof(runargs[kk+1]);
      info=DSDPSetGapTolerance(dsdp,gaptol);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-penalty",7)==0){
      penalty=atof(runargs[kk+1]);
      info=DSDPSetPenaltyParameter(dsdp,penalty); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-bigM",5)==0){
      rpos=atoi(runargs[kk+1]);
      info=DSDPUsePenalty(dsdp,rpos); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-maxit",6)==0){
      maxit=atoi(runargs[kk+1]);
      info=DSDPSetMaxIts(dsdp,maxit); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-r0",3)==0){
      r0=atof(runargs[kk+1]);
      info=DSDPSetR0(dsdp,r0); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-zbar",5)==0){
      zbar=atof(runargs[kk+1]);
      info=DSDPSetZBar(dsdp,zbar);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-infdtol",7)==0){
      inftol=atof(runargs[kk+1]);
      info=DSDPSetRTolerance(dsdp,inftol);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-infptol",7)==0){
      infptol=atof(runargs[kk+1]);
      info=DSDPSetPTolerance(dsdp,infptol);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-rho",4)==0){
      rho=atof(runargs[kk+1]);
      info=DSDPSetPotentialParameter(dsdp,rho); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-drho",5)==0){
      drho=atoi(runargs[kk+1]);
      info=DSDPUseDynamicRho(dsdp,drho);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-mu0",4)==0){
      mu0=atof(runargs[kk+1]);
      info=DSDPSetBarrierParameter(dsdp,mu0);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-maxtrustradius",7)==0){
      maxtrust=atof(runargs[kk+1]);
      info=DSDPSetMaxTrustRadius(dsdp,maxtrust); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-boundy",6)==0){
      yhigh=fabs(atof(runargs[kk+1]));ylow=-yhigh;
      info=DSDPSetYBounds(dsdp,ylow,yhigh);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-steptol",7)==0){
      steptol=fabs(atof(runargs[kk+1]));
      info=DSDPSetStepTolerance(dsdp,steptol); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-pnormtol",7)==0){
      pnormtol=fabs(atof(runargs[kk+1]));
      info=DSDPSetPNormTolerance(dsdp,pnormtol); DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-reuse",6)==0){
      reuse=atoi(runargs[kk+1]);
      info=DSDPReuseMatrix(dsdp,reuse);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-dadd",6)==0){
      cc=atof(runargs[kk+1]);
      info=DSDPAddObjectiveConstant(dsdp,cc);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-dbound",6)==0){
      dbound=atof(runargs[kk+1]);
      info=DSDPSetDualBound(dsdp,dbound);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-fix",4)==0){
      info=DSDPSetFixedVariable(dsdp,1,atof(runargs[kk+1]));DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-dobjmin",7)==0){
      dd=atof(runargs[kk+1]);
      info = DSDPSetDualLowerBound(dsdp,dd);DSDPCHKERR(info);
    } else if (strncmp(runargs[kk],"-dloginfo",8)==0){
      iloginfo=atoi(runargs[kk+1]);
      info=DSDPLogInfoAllow(iloginfo,0);
    }
  }

  for (kk=0; kk<nargs; kk++){
    if (0){
    } else if (strncmp(runargs[kk],"-help",5)==0){
      info=DSDPPrintOptions();
    }
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPReadOptions(DSDP dsdp, char filename[]);
\brief Read DSDP parameters from a file.
\param dsdp is the solver
\param filename is the name of a file
\sa DSDPSetOptions()
\sa DSDPPrintOptions()
\sa DSDPView()
\ingroup DSDPSolver
*/
#define MAXOPTIONS 40
#define STRLENGTH  40
#define BUFFERSIZ 100
#undef __FUNCT__
#define __FUNCT__ "DSDPReadOptions"
int DSDPReadOptions(DSDP dsdp, char filename[]){

  int i,info,line=0;
  char thisline[BUFFERSIZ]="%",doption[STRLENGTH],dvalue[STRLENGTH];
  char fargs[2*MAXOPTIONS][STRLENGTH];
  char *fargs2[2*MAXOPTIONS];
  FILE *fp;

  DSDPFunctionBegin;

  for (i=0;i<2*MAXOPTIONS;i++){fargs2[i]=fargs[i];}
  
  fp=fopen(filename,"r");
  if (fp){
    while(!feof(fp) ){
      if (line>=MAXOPTIONS) break;
      fgets(thisline,BUFFERSIZ,fp);
      if (sscanf(thisline,"%s %s",doption,dvalue)>=2){
	if (doption[0]!='%'){
	  strncpy(fargs[2*line],doption,STRLENGTH-1);
	  strncpy(fargs[2*line+1],dvalue,STRLENGTH-1);
	  line++;
	}
      }
      thisline[0]='%';
    }
    
    info=DSDPSetOptions(dsdp,fargs2,2*line);    
    fclose(fp);
  }
  DSDPFunctionReturn(0);
}

/*!
\fn int DSDPView(DSDP dsdp);
\brief Print many of the parameters currently set in DSDP.
\param dsdp is the solver
\sa DSDPSetOptions()
\sa DSDPGetPenaltyParameter()
\sa DSDPGetSolutionType()
\sa DSDPGetGapTolerance()
\ingroup DSDPBasic
*/
#undef __FUNCT__
#define __FUNCT__ "DSDPView"
int DSDPView(DSDP dsdp){

  int info,reuse,m,maxit;
  double penalty,rho,mu0,gaptol,dnorm[3],derror[6],potential,ymax;
  double ylow,yhigh,maxtrust,steptol,inftol,infptol,pnormtol,dbound,tracex;
  DSDPSolutionType pdfeasible;

  DSDPFunctionBegin;
  info=DSDPGetMaxIts(dsdp,&maxit); DSDPCHKERR(info);
  printf("Terminate DSDP after %d iterations.\n",maxit);
  info=DSDPGetDualBound(dsdp,&dbound); DSDPCHKERR(info);
  printf("Terminate DSDP if dual objective is greater than %8.4e\n",dbound);
  info=DSDPGetGapTolerance(dsdp,&gaptol);DSDPCHKERR(info);
  printf("Terminate DSDP if the relative duality gap is less than %8.4e\n",gaptol);
  info=DSDPGetStepTolerance(dsdp,&steptol); DSDPCHKERR(info);
  printf("Terminate DSDP if step length in D less than %8.4e\n",steptol);
  info=DSDPGetPNormTolerance(dsdp,&pnormtol); DSDPCHKERR(info);
  printf("Terminate DSDP only if Pnorm less than %8.4e\n",pnormtol);
  info=DSDPGetMaxTrustRadius(dsdp,&maxtrust); DSDPCHKERR(info);
  printf("Max Trust Radius is %8.4e\n",maxtrust); 
  info=DSDPGetReuseMatrix(dsdp,&reuse);DSDPCHKERR(info);
  printf("Reapply Hessian of Barrier up to %d times per iteration.\n",reuse);

  info=DSDPGetDataNorms(dsdp,dnorm);DSDPCHKERR(info);
  printf("The norms of C: %8.4e, A: %4.4e, and b: %8.4e\n",dnorm[0],dnorm[1],dnorm[2]);
  info=DSDPGetNumberOfVariables(dsdp,&m);DSDPCHKERR(info);
  printf("There are %d y variables:  ",m);
  info=DSDPGetYMaxNorm(dsdp,&ymax); DSDPCHKERR(info);
  printf("largest is %8.4e, ",ymax);
  info=DSDPGetYBounds(dsdp,&ylow,&yhigh);DSDPCHKERR(info);
  printf("bounded below by %8.4e and above by %8.4e. \n",ylow,yhigh);
  info=DSDPGetTraceX(dsdp,&tracex);DSDPCHKERR(info);
  printf("The X variables have a trace of %8.4e ",tracex);
  info=DSDPGetPenaltyParameter(dsdp,&penalty); DSDPCHKERR(info);
  printf("bounded by penalty parameter: %8.4e\n",penalty);
  info=DSDPGetBarrierParameter(dsdp,&mu0);DSDPCHKERR(info);
  printf("Current Barrier Parameter: %8.4e\n",mu0);
  info=DSDPGetPotentialParameter(dsdp,&rho); DSDPCHKERR(info);
  printf("Potential Parameter: %8.4e ( times dimension) \n",rho);
  info=DSDPGetPotential(dsdp,&potential);DSDPCHKERR(info);
  printf("The value of the potential function is %8.4e\n",potential);
  info=DSDPGetRTolerance(dsdp,&inftol); DSDPCHKERR(info);
  printf("(D) Feasible only if R < %8.4e\n",inftol);
  info=DSDPGetPTolerance(dsdp,&infptol); DSDPCHKERR(info);
  printf("(P) Feasible only if Pinfeas < %8.4e\n",infptol);
  info=DSDPGetSolutionType(dsdp,&pdfeasible);DSDPCHKERR(info);
  if (pdfeasible==DSDP_PDFEASIBLE){
    printf(" DSDP Solutions are both feasible and bounded\n");
  } else if (pdfeasible==DSDP_UNBOUNDED){
    printf(" (D) is unbounded and (P) is infeasible\n");
  } else if (pdfeasible==DSDP_INFEASIBLE){
    printf(" (D) is infeasible and (D) is unbounded\n");
  } else if (pdfeasible==DSDP_PDUNKNOWN){
    printf(" Hmm.  Not clear whether either solution is feasible.\n");
  }
  info=DSDPGetFinalErrors(dsdp,derror);DSDPCHKERR(info);
  printf("The errors: %8.4e, %4.4e, %8.4e, ",derror[0],derror[1],derror[2]);
  printf("%8.4e, %4.4e, %8.4e\n",derror[3],derror[4],derror[5]);
  DSDPFunctionReturn(0);
}


static char dsdpoptions[]="\
                  -gaptol <1e-6> stop when relative duality gap less than  \n\
                  -r0 <-1> if nonnegative, initialize S by adding this multiple of the identity matrix \n\
                  -penalty <1e10>< penalize dual infeasibility \n\
                  -boundy <1e7> bound for variables y \n\
                  -maxit <200> set maximum iterates \n\
                  -zbar <1e10> Upper bound for dual solution \n\
                  -mu0 <-1> if positive, set initial barrier parameter \n\
                  -rho <3> Potential parameter as multiple of dimension \n\
                  -drho <1> Use dynamic rho strategy \n\
                  -pnormtol <1e30> stop only if pnorm less than \n\
                  -reuse <4> Reuse the Schur Matrix this many times\n\
                  -dobjmin <> apply a known lower bound for the objective at solution as a constraint. \n\
                  -bigM <0> if positive, modify algorithm to make dual \n\
                     infeasibility positive with a large associated cost \n\
                  -dloginfo <0> - print more information for higher numbers \n\
                  -params <filename> to read selected options from a file \n\
                  -help  for this help message\n";

/*!
\fn DSDPPrintOptions();
\brief Print runtime options;
\sa DSDPSetOptions()
\sa DSDPReadOptions()
\sa DSDPView()
\ingroup DSDPSolver
 */
int DSDPPrintOptions(){
  DSDPFunctionBegin;
  printf("%s",dsdpoptions); 
  DSDPFunctionReturn(0);
}
