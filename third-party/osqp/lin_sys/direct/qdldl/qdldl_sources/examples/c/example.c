#include "qdldl.h"
#include <stdlib.h>
#include <stdio.h>

void print_arrayi(const QDLDL_int* data, QDLDL_int n,char* varName);
void print_arrayf(const QDLDL_float* data, QDLDL_int n, char* varName);
void print_line(void);

const QDLDL_int   An   = 10;
const QDLDL_int   Ap[] = {0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 17};
const QDLDL_int   Ai[] = {0, 1, 1, 2, 3, 4, 1, 5, 0, 6, 3, 7, 6, 8, 1, 2, 9};
const QDLDL_float Ax[] = {1.0, 0.460641, -0.121189, 0.417928, 0.177828, 0.1,
                       -0.0290058, -1.0, 0.350321, -0.441092, -0.0845395,
                       -0.316228, 0.178663, -0.299077, 0.182452, -1.56506, -0.1};
const QDLDL_float  b[] = {1,2,3,4,5,6,7,8,9,10};

int main()

{
  QDLDL_int i; // Counter

  //data for L and D factors
  QDLDL_int Ln = An;
  QDLDL_int *Lp;
  QDLDL_int *Li;
  QDLDL_float *Lx;
  QDLDL_float *D;
  QDLDL_float *Dinv;

  //data for elim tree calculation
  QDLDL_int *etree;
  QDLDL_int *Lnz;
  QDLDL_int  sumLnz;

  //working data for factorisation
  QDLDL_int   *iwork;
  QDLDL_bool  *bwork;
  QDLDL_float *fwork;

  //Data for results of A\b
  QDLDL_float *x;


  /*--------------------------------
   * pre-factorisation memory allocations
   *---------------------------------*/

  //These can happen *before* the etree is calculated
  //since the sizes are not sparsity pattern specific

  //For the elimination tree
  etree = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);
  Lnz   = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);

  //For the L factors.   Li and Lx are sparsity dependent
  //so must be done after the etree is constructed
  Lp    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(An+1));
  D     = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);
  Dinv  = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);

  //Working memory.  Note that both the etree and factor
  //calls requires a working vector of QDLDL_int, with
  //the factor function requiring 3*An elements and the
  //etree only An elements.   Just allocate the larger
  //amount here and use it in both places
  iwork = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(3*An));
  bwork = (QDLDL_bool*)malloc(sizeof(QDLDL_bool)*An);
  fwork = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);

  /*--------------------------------
   * elimination tree calculation
   *---------------------------------*/
  sumLnz = QDLDL_etree(An,Ap,Ai,iwork,Lnz,etree);

  /*--------------------------------
   * LDL factorisation
   *---------------------------------*/

  //First allocate memory for Li and Lx
  Li    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*sumLnz);
  Lx    = (QDLDL_float*)malloc(sizeof(QDLDL_float)*sumLnz);

  //now factor
  QDLDL_factor(An,Ap,Ai,Ax,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork);

  /*--------------------------------
   * solve
   *---------------------------------*/
  x = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);

  //when solving A\b, start with x = b
  for(i=0;i < Ln; i++) x[i] = b[i];
  QDLDL_solve(Ln,Lp,Li,Lx,Dinv,x);

  /*--------------------------------
   * print factors and solution
   *---------------------------------*/
  printf("\n");
  printf("A (CSC format):\n");
  print_line();
  print_arrayi(Ap, An + 1, "A.p");
  print_arrayi(Ai, Ap[An], "A.i");
  print_arrayf(Ax, Ap[An], "A.x");
  printf("\n\n");

  printf("elimination tree:\n");
  print_line();
  print_arrayi(etree, Ln, "etree");
  print_arrayi(Lnz, Ln, "Lnz");
  printf("\n\n");

  printf("L (CSC format):\n");
  print_line();
  print_arrayi(Lp, Ln + 1, "L.p");
  print_arrayi(Li, Lp[Ln], "L.i");
  print_arrayf(Lx, Lp[Ln], "L.x");
  printf("\n\n");

  printf("D:\n");
  print_line();
  print_arrayf(D, An,    "diag(D)     ");
  print_arrayf(Dinv, An, "diag(D^{-1})");
  printf("\n\n");

  printf("solve results:\n");
  print_line();
  print_arrayf(b, An, "b");
  print_arrayf(x, An, "A\\b");
  printf("\n\n");


  /*--------------------------------
   * clean up
   *---------------------------------*/
  free(Lp);
  free(Li);
  free(Lx);
  free(D);
  free(Dinv);
  free(etree);
  free(Lnz);
  free(iwork);
  free(bwork);
  free(fwork);
  free(x);

return 0 ;


}


void print_line(void){
  printf("--------------------------\n");
}

void print_arrayi(const QDLDL_int* data, QDLDL_int n,char* varName){

  QDLDL_int i;
  printf("%s = [",varName);
  for(i=0; i< n; i++){
    printf("%i,",(int)data[i]);
  }
  printf("]\n");

}

void print_arrayf(const QDLDL_float* data, QDLDL_int n, char* varName){

  QDLDL_int i;
  printf("%s = [",varName);
  for(i=0; i< n; i++){
    printf("%.3g,",data[i]);
  }
  printf("]\n");

}
