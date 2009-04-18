#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

/* per il diametro */
void ld2Ddiamsimpsim(double *x1, double *x2, double *y1, double *y2, double *tau, int *nx, int *ny, int *nt, double *depth, double *ldepth) {
  double d1, d2, d3, dy1, dy2, dy3, a1, a2, a3, b1, b2, b3;
  int ind1, ind2, ind3, i, i1, i2, i3; 
  ind1 = -1;
  ind3 = -1;
  for (i1 = 0 ; i1 < (*nx-2) ; i1++) {
    for (i2 = (i1+1) ; i2 < (*nx-1) ; i2++) {
      for (i3 = (i2+1); i3 < *nx ; i3++) {
        ind1 += 1;
        d1 = sqrt(R_pow((x1[i1]-x1[i2]), 2.0)+R_pow((x2[i1]-x2[i2]), 2.0));
        d2 = sqrt(R_pow((x1[i2]-x1[i3]), 2.0)+R_pow((x2[i2]-x2[i3]), 2.0));
        d3 = sqrt(R_pow((x1[i1]-x1[i3]), 2.0)+R_pow((x2[i1]-x2[i3]), 2.0));
        for (ind2 = 0 ; ind2 < *ny ; ind2++) {
	  ind3 +=1;
          dy1 = sqrt(R_pow((x1[i1]-y1[ind2]), 2.0)+R_pow((x2[i1]-y2[ind2]), 2.0));
          dy2 = sqrt(R_pow((x1[i2]-y1[ind2]), 2.0)+R_pow((x2[i2]-y2[ind2]), 2.0));
          dy3 = sqrt(R_pow((x1[i3]-y1[ind2]), 2.0)+R_pow((x2[i3]-y2[ind2]), 2.0));

          if (dy1==0.0 | dy2==0.0 | dy3==0.0) {
            depth[ind3] = 1.0;
            if (d1 <= *tau & d2 <= *tau & d3 <= *tau) {
              ldepth[ind3] = 1.0;
            }
          } else {
            // printf("dy %6.3f %6.3f %6.3f \n", dy1, dy2, dy3);

	    a1 = atan2((x1[i1]-y1[ind2]), (x2[i1]-y2[ind2]));
	    a2 = atan2((x1[i2]-y1[ind2]), (x2[i2]-y2[ind2]));
	    a3 = atan2((x1[i3]-y1[ind2]), (x2[i3]-y2[ind2]));

          // printf("atan %6.3f %6.3f %6.3f \n", a1, a2, a3);

            if (a1 < 0.0)  a1 += 2.0*M_PI;
            if (a2 < 0.0)  a2 += 2.0*M_PI;
            if (a3 < 0.0)  a3 += 2.0*M_PI;

          // printf("a %6.3f %6.3f %6.3f \n", a1, a2, a3);

            if (a1 > a2 & a1 > a3) {
	      b3 = a1;
              if (a2 > a3) {
	        b1 = a3;
                b2 = a2;
              } else {
  	        b1 = a2;
                b2 = a3;             
              }
	    } else if (a2 > a3) {
              b3 = a2;
              if (a1 > a3) {
                b1 = a3;
                b2 = a1;
              } else {
                b1 = a1;
                b2 = a3;
              }
            } else {
              b3 = a3;
              if (a1 > a2) {
                b1 = a2;
                b2 = a1;
              } else {
                b1 = a1;
                b2 = a2;
              }
            }
          // printf("b %6.3f %6.3f %6.3f \n", b1, b2, b3);

            a1 = b3-b1;
            a2 = b1-b2;
            a3 = b2-b3;
            if (a1 < 0.0)  a1 += 2.0*M_PI;
            if (a2 < 0.0)  a2 += 2.0*M_PI;
            if (a3 < 0.0)  a3 += 2.0*M_PI;

          // printf('angle %6.3f %6.3f %6.3f \n', a1, a2, a3);

            if (a1 >= M_PI & a2 >= M_PI & a3 >= M_PI) {
            // printf("%6.3f %6.3f %6.3f \n", a1, a2, a3);
	      depth[ind3] = 1.0; 
               if (d1 <= *tau & d2 <= *tau & d3 <= *tau) {
                 ldepth[ind3] = 1.0;
               }             
            }
          }
        }
      }
    }
  }
  //  ind3 = -1;
  //for (ind2 = 0 ; ind2 < (*ny-2) ; ind2++) {
  //for (ind1 = 0 ; (ind2+1) < (*ny-1) ; ind2++) {
  //   ind3 += 1;
  //   simdepth[ind3] = 0.0;
  //   simldepth[ind3] = 0.0;
  //   for (int i=0;i< *nt;i++) {
  //	 simdepth[ind3] += depth[ind2][i] * depth[ind1][i];
  //	 simldepth[ind3] += ldepth[ind2][i] * ldepth[ind1][i];
  //   }
  //}
  //}
}


/* per l'area dei simplessi */
void ld2Dareasimpsim(double *x1, double *x2, double *y1, double *y2, double *tau, int *nx, int *ny, int *nt, double *depth, double *ldepth) {
  double dsum, area, d1, d2, d3, dy1, dy2, dy3, a1, a2, a3, b1, b2, b3;
  int ind1, ind2, ind3, i1, i2, i3; 
  ind1 = -1;
  ind3 = -1;
  for (i1 = 0 ; i1 < (*nx-2) ; i1++) {
    for (i2 = (i1+1) ; i2 < (*nx-1) ; i2++) {
      for (i3 = (i2+1); i3 < *nx ; i3++) {
        ind1 += 1;
        d1 = sqrt(R_pow((x1[i1]-x1[i2]), 2.0)+R_pow((x2[i1]-x2[i2]), 2.0));
        d2 = sqrt(R_pow((x1[i2]-x1[i3]), 2.0)+R_pow((x2[i2]-x2[i3]), 2.0));
        d3 = sqrt(R_pow((x1[i1]-x1[i3]), 2.0)+R_pow((x2[i1]-x2[i3]), 2.0));
        dsum = (d1+d2+d3)/2.0;
        area = sqrt(dsum*(dsum-d1)*(dsum-d2)*(dsum-d3));

        for (ind2 = 0 ; ind2 < *ny ; ind2++) {
	  ind3 += 1;
          dy1 = sqrt(R_pow((x1[i1]-y1[ind2]), 2.0)+R_pow((x2[i1]-y2[ind2]), 2.0));
          dy2 = sqrt(R_pow((x1[i2]-y1[ind2]), 2.0)+R_pow((x2[i2]-y2[ind2]), 2.0));
          dy3 = sqrt(R_pow((x1[i3]-y1[ind2]), 2.0)+R_pow((x2[i3]-y2[ind2]), 2.0));

          if (dy1==0.0 | dy2==0.0 | dy3==0.0) {
            depth[ind3] = 1.0;
            if (area <= *tau) {
              ldepth[ind3] = 1.0;
            }
	    // printf("spigolo %6.3f %6.3f %6.3f %6.3f \n", area, *tau, depth[ind3], ldepth[ind3]);
          } else {
            // printf("dy %6.3f %6.3f %6.3f \n", dy1, dy2, dy3);

	    a1 = atan2((x1[i1]-y1[ind2]), (x2[i1]-y2[ind2]));
	    a2 = atan2((x1[i2]-y1[ind2]), (x2[i2]-y2[ind2]));
	    a3 = atan2((x1[i3]-y1[ind2]), (x2[i3]-y2[ind2]));

          // printf("atan %6.3f %6.3f %6.3f \n", a1, a2, a3);

            if (a1 < 0.0)  a1 += 2.0*M_PI;
            if (a2 < 0.0)  a2 += 2.0*M_PI;
            if (a3 < 0.0)  a3 += 2.0*M_PI;

          // printf("a %6.3f %6.3f %6.3f \n", a1, a2, a3);

            if (a1 > a2 & a1 > a3) {
	      b3 = a1;
              if (a2 > a3) {
	        b1 = a3;
                b2 = a2;
              } else {
  	        b1 = a2;
                b2 = a3;             
              }
	    } else if (a2 > a3) {
              b3 = a2;
              if (a1 > a3) {
                b1 = a3;
                b2 = a1;
              } else {
                b1 = a1;
                b2 = a3;
              }
            } else {
              b3 = a3;
              if (a1 > a2) {
                b1 = a2;
                b2 = a1;
              } else {
                b1 = a1;
                b2 = a2;
              }
            }
          // printf("b %6.3f %6.3f %6.3f \n", b1, b2, b3);

            a1 = b3-b1;
            a2 = b1-b2;
            a3 = b2-b3;
            if (a1 < 0.0)  a1 += 2.0*M_PI;
            if (a2 < 0.0)  a2 += 2.0*M_PI;
            if (a3 < 0.0)  a3 += 2.0*M_PI;

          // printf('angle %6.3f %6.3f %6.3f \n', a1, a2, a3);

            if (a1 >= M_PI & a2 >= M_PI & a3 >= M_PI) {
            // printf("%6.3f %6.3f %6.3f \n", a1, a2, a3);
	      depth[ind3] = 1.0; 
               if (area <= *tau) {
                 ldepth[ind3] = 1.0;
               }
	       // printf("no spi. %6.3f %6.3f %6.3f %6.3f \n", area, *tau, depth[ind3], ldepth[ind3]);
             
            }
          }
        }
      }
    }
  }
}
