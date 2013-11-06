
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>
#include <stdio.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif


#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))
#ifdef R_160_and_older
#define both_non_NA both_FINITE
#else
#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#endif


void ldcircsimp(double *x, double *y, int *nx, int *ny, double *tau,
      int *nuse, double *depth, double *ldepth, double *diameter) {

  int amb, ind1, ind2, i1, i2; 
  double z, x1, x2, d12, d21, spherical;

  for (ind2 = 0 ; ind2 < *ny ; ind2++) {
    depth[ind2] = 0.0;
    ldepth[ind2] = 0.0;
  }

  ind1 = -1;
  for (i1 = 0 ; i1 < *nx-1 ; i1++) {
    for (i2 = (i1+1) ; i2 < *nx ; i2++) {
      x1 = x[i1];
      x2 = x[i2];
      ind1 += 1;
      d21 = fabs(fmod((x2-x1+2.0*M_PI), (2.0*M_PI)));
      d12 = fabs(fmod((x1-x2+2.0*M_PI), (2.0*M_PI)));
      amb = 0;
      if (d21 == d12 && d12 == M_PI) {
        diameter[ind1] = d21;
        amb = 1;
      } else if (d21 < d12) {
        diameter[ind1] = d21;
      } else {
        diameter[ind1] = d12;
        z = x1;
        x1 = x2;
        x2 = z;
      }
      for (ind2 = 0 ; ind2 < *ny ; ind2++) {
        if (amb == 1) {
          depth[ind2] += 0.5;
	  if (*tau >= M_PI) {
            ldepth[ind2] += 0.5;
	  }
        } else { 
          if (fabs(fmod((y[ind2]-x1+2.0*M_PI), (2.0*M_PI))) <= diameter[ind1]) {
	    depth[ind2] += 1.0;
            if (*nuse == 1 | *nuse == 2) {
	      if (diameter[ind1] <= *tau) {
	        ldepth[ind2] += 1.0;
              }
            } else {
              spherical = fmax(fabs(fmod((y[ind2]-x1+2.0*M_PI), (2.0*M_PI))), fabs(fmod((x2-y[ind2]+2.0*M_PI), (2.0*M_PI))));
              if (spherical <= *tau) {
		ldepth[ind2] += 1.0;
              }
            }
          }
        }
      }
    }
  }
}
