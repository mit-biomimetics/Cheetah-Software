#include "dsdplapack.h"

typedef long int integer;
typedef long int logical;

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmax(a,b) (double)max(a,b)

static int lsame2(const char c1[], const char c2[]){
  if (c1[0]==c2[0]) return 1;
  else return 0;
}

/* Subroutine */ int dtrsm2(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, double *alpha, double *a, integer *
	lda, double *b, integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i1, i2, i3;
    /* Local variables */
    static integer info;
    static double temp,alpha2,aref;
    static integer i, j, k;
    static logical lside;
    static integer nrowa;
    static logical upper;
    static logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    /* Function Body */
    lside = lsame2(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame2(diag, "N");
    upper = lsame2(uplo, "U");
    info = 0;
    if (! lside && ! lsame2(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame2(uplo, "L")) {
	info = 2;
    } else if (! lsame2(transa, "N") && ! lsame2(transa,
	     "T") && ! lsame2(transa, "C")) {
	info = 3;
    } else if (! lsame2(diag, "U") && ! lsame2(diag, 
	    "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	return info;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.) {
	i1 = *n;
	for (j = 1; j <= i1; ++j) {
	    i2 = *m;
	    for (i = 1; i <= i2; ++i) {
		b_ref(i, j) = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }
/*     Start the operations. */
    if (lside) {
      if (lsame2(transa, "N")) {
	/*           Form  B := alpha*inv( A )*B. */
	if (upper) {
	  i1 = *n;
	  for (j = 1; j <= i1; ++j) {
	    if (*alpha != 1.) {
	      i2 = *m;
	      alpha2=*alpha;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, j) *= alpha2;
		/* L30: */
	      }
	    }
	    for (k = *m; k >= 1; --k) {
	      if (b_ref(k, j) != 0.) {
		if (nounit) {
		  b_ref(k, j) /=  a_ref(k, k);
		}
		i2 = k - 1;
		aref=-b_ref(k, j);
		for (i = 1; i <= i2; ++i) {
		  b_ref(i, j) += aref * a_ref(i, k);
/* L40: */
		}
			}
	      /* L50: */
	    }
	    /* L60: */
	  }
	} else {
	  i1 = *n;
	  for (j = 1; j <= i1; ++j) {
	    if (*alpha != 1.) {
	      i2 = *m;
	      alpha2=*alpha;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, j) *= alpha2;
		/* L70: */
	      }
	    }
	    i2 = *m;
	    for (k = 1; k <= i2; ++k) {
	      if (b_ref(k, j) != 0.) {
		if (nounit) {
		  b_ref(k, j) /=  a_ref(k, k);
		}
		i3 = *m;
		aref= -b_ref(k, j);
		for (i = k + 1; i <= i3; ++i) {
		  b_ref(i, j) += aref * a_ref(i, k);
		  /* L80: */
		}
	      }
	      /* L90: */
	    }
	    /* L100: */
	  }
	}
      } else {
	/*           Form  B := alpha*inv( A' )*B. */
	if (upper) {
	  i1 = *n;
	  for (j = 1; j <= i1; ++j) {
	    i2 = *m;
	    alpha2=*alpha;
	    for (i = 1; i <= i2; ++i) {
	      temp = alpha2 * b_ref(i, j);
	      i3 = i - 1;
	      for (k = 1; k <= i3; ++k) {
		temp -= a_ref(k, i) * b_ref(k, j);
		/* L110: */
	      }
	      if (nounit) {
		temp /= a_ref(i, i);
	      }
	      b_ref(i, j) = temp;
/* L120: */
	    }
	    /* L130: */
	  }
	} else {
	  i1 = *n;
	  for (j = 1; j <= i1; ++j) {
	    for (i = *m; i >= 1; --i) {
	      temp = alpha2 * b_ref(i, j);
	      i2 = *m;
	      for (k = i + 1; k <= i2; ++k) {
		temp -= a_ref(k, i) * b_ref(k, j);
		/* L140: */
	      }
	      if (nounit) {
		temp /= a_ref(i, i);
	      }
	      b_ref(i, j) = temp;
	      /* L150: */
	    }
/* L160: */
	  }
	}
	}
    } else {
      if (lsame2(transa, "N")) {
	/*           Form  B := alpha*B*inv( A ). */
	if (upper) {
	  i1 = *n;
	  for (j = 1; j <= i1; ++j) {
	    if (*alpha != 1.) {
	      i2 = *m;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, j) *= alpha2;
		/* L170: */
	      }
	    }
	    i2 = j - 1;
	    for (k = 1; k <= i2; ++k) {
	      if (a_ref(k, j) != 0.) {
		i3 = *m;
		aref=-a_ref(k, j);
		for (i = 1; i <= i3; ++i) {
		  b_ref(i, j) += aref *  b_ref(i, k);
		  /* L180: */
		}
	      }
	      /* L190: */
	    }
	    if (nounit) {
	      temp = 1. / a_ref(j, j);
	      i2 = *m;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, j) *= temp;
		/* L200: */
	      }
	    }
	    /* L210: */
	  }
	} else {
	  for (j = *n; j >= 1; --j) {
	    if (*alpha != 1.) {
	      i1 = *m;
	      for (i = 1; i <= i1; ++i) {
		b_ref(i, j) *= alpha2;
		/* L220: */
	      }
		    }
	    i1 = *n;
	    for (k = j + 1; k <= i1; ++k) {
	      if (a_ref(k, j) != 0.) {
		i2 = *m;
		aref=-a_ref(k, j);
		for (i = 1; i <= i2; ++i) {
		  b_ref(i, j) += aref *  b_ref(i, k);
		  /* L230: */
		}
	      }
	      /* L240: */
	    }
	    if (nounit) {
	      temp = 1. / a_ref(j, j);
	      i1 = *m;
	      for (i = 1; i <= i1; ++i) {
		b_ref(i, j) *= temp;
		/* L250: */
	      }
	    }
	    /* L260: */
	  }
	}
      } else {
	/*           Form  B := alpha*B*inv( A' ). */
	if (upper) {
	  for (k = *n; k >= 1; --k) {
	    if (nounit) {
	      temp = 1. / a_ref(k, k);
	      i1 = *m;
	      for (i = 1; i <= i1; ++i) {
		b_ref(i, k) *= temp;
		/* L270: */
	      }
	    }
	    i1 = k - 1;
	    for (j = 1; j <= i1; ++j) {
	      if (a_ref(j, k) != 0.) {
		temp = a_ref(j, k);
		i2 = *m;
		for (i = 1; i <= i2; ++i) {
		  b_ref(i, j) -= temp * b_ref(i, k);
		  /* L280: */
		}
	      }
	      /* L290: */
	    }
	    if (*alpha != 1.) {
	      i1 = *m;
	      for (i = 1; i <= i1; ++i) {
		b_ref(i, k) *= alpha2;
		/* L300: */
	      }
	    }
	    /* L310: */
	  }
	} else {
	  i1 = *n;
	  for (k = 1; k <= i1; ++k) {
	    if (nounit) {
	      temp = 1. / a_ref(k, k);
	      i2 = *m;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, k) *= temp;
		/* L320: */
	      }
	    }
	    i2 = *n;
	    for (j = k + 1; j <= i2; ++j) {
	      if (a_ref(j, k) != 0.) {
		temp = a_ref(j, k);
		i3 = *m;
		for (i = 1; i <= i3; ++i) {
		  b_ref(i, j) -= temp * b_ref(i, k);
		  /* L330: */
		}
	      }
	      /* L340: */
		    }
	    if (*alpha != 1.) {
	      i2 = *m;
	      for (i = 1; i <= i2; ++i) {
		b_ref(i, k) *= alpha2;
		/* L350: */
	      }
	    }
	    /* L360: */
	  }
	}
      }
    }
    return 0;
    /*     End of DTRSM . */
} /* dtrsm_ */
#undef b_ref
#undef a_ref

