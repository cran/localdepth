#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>


/* 2D DIAMETRO */
void twoDdiam(double *x, double *y, int *nx, double *diameter) {
  double d1, d2, d3;
  int ind1, i1, i2, i3; 
  ind1 = -1;
  for (i1 = 0 ; i1 < (*nx-2) ; i1++) {
    for (i2 = (i1+1) ; i2 < (*nx-1) ; i2++) {
      for (i3 = (i2+1); i3 < *nx ; i3++) {
        ind1 += 1;
        d1 = sqrt(R_pow((x[i1]-x[i2]), 2.0)+R_pow((y[i1]-y[i2]), 2.0));
        d2 = sqrt(R_pow((x[i2]-x[i3]), 2.0)+R_pow((y[i2]-y[i3]), 2.0));
        d3 = sqrt(R_pow((x[i1]-x[i3]), 2.0)+R_pow((y[i1]-y[i3]), 2.0));
        diameter[ind1] = d3;
        if (d1 > d2) {
          if (d1 > d3) {
            diameter[ind1] = d1;
          }
        } else {
          if (d2 > d3) {
            diameter[ind1] = d2;
          }
        }
      }
    }
  }
}

/* 2D AREA */
void twoDarea(double *x, double *y, int *nx, double *area) {
  double dsum, d1, d2, d3, dtemp;
  int ind1, i1, i2, i3; 
  ind1 = -1;
  for (i1 = 0 ; i1 < (*nx-2) ; i1++) {
    for (i2 = (i1+1) ; i2 < (*nx-1) ; i2++) {
      for (i3 = (i2+1); i3 < *nx ; i3++) {
        ind1 += 1;
        d1 = sqrt(R_pow((x[i1]-x[i2]), 2.0)+R_pow((y[i1]-y[i2]), 2.0));
        d2 = sqrt(R_pow((x[i2]-x[i3]), 2.0)+R_pow((y[i2]-y[i3]), 2.0));
        d3 = sqrt(R_pow((x[i1]-x[i3]), 2.0)+R_pow((y[i1]-y[i3]), 2.0));
        dsum = (d1+d2+d3)/2.0;
        dtemp = dsum*(dsum-d1)*(dsum-d2)*(dsum-d3);
	area[ind1] = sqrt(fmax(0.0, dtemp));
      }
    }
  }
}



