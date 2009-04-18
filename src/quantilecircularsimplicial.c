
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void circdiam(double *x, int *nx, double *diameter) {

  int ind1, i1, i2; 
  double d12, d21;
  ind1 = -1;
  for (i1 = 0 ; i1 < *nx-1 ; i1++) {
    for (i2 = (i1+1) ; i2 < *nx ; i2++) {
      ind1 += 1;
      d21 = fabs(fmod((x[i2]-x[i1]+2.0*M_PI), (2.0*M_PI)));
      d12 = fabs(fmod((x[i1]-x[i2]+2.0*M_PI), (2.0*M_PI)));
      if (d21 >= d12)
        diameter[ind1] = d12;
      else
        diameter[ind1] = d21;
    }
  }
}
