
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void oneDdiam(double *x, int *nx, double *diameter) {
  int ind1, i1, i2; 
  ind1 = -1;
  for (i1 = 0 ; i1 < *nx-1 ; i1++) {
    for (i2 = (i1+1) ; i2 < *nx ; i2++) {
      ind1 += 1;
      diameter[ind1] = fabs(x[i2]-x[i1]);
    }
  }
}

