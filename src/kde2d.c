#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>



/*  Compute the density of the normal distribution. (fast dnorm) */
double fast_dnorm( double x ){

  const double sigma = 1.;
  //const double M_1_SQRT_2PI = 1./sqrt(2.* M_PI);

  return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma ;

}



/* Two-Dimensional Kernel Density Estimation
 

Description: Two dimensional kernel density estimation with an axis-aligned
bivariate normal kernel, evaluated on a square grid

*/


void kde2d ( int nxy, int n, double *x, double *y, 
                      double h_in, double *lims,
                      double *gx, double *gy, double *z){


   double deltax = (lims[1] - lims[0]) / (n - 1.);
   double deltay = (lims[3] - lims[2]) / (n - 1.);

   int i, j, k;
   double h = h_in/4.;
   double sum;

   for (i = 0; i < n; i ++){

      gx[i] = lims[0] + i * deltax;
      gy[i] = lims[2] + i * deltay;
   }
   gx[n-1] = lims[1];
   gy[n-1] = lims[3];
   

   for(j = 0 ; j < n; j++){   // number of elements in gy
      for (i = 0; i < n; i++){      // number of elements in gx

         sum = 0;
         for(k = 0; k < nxy; k++){   // number of elements in x znd y

            //printf("%d %d %d %12.8f  %12.8f\n", i, j, k , fast_dnorm((gy[j]-y[k])/h) ,fast_dnorm( (gx[i] - x[k])/h)); 
            sum += fast_dnorm( (gy[j] - y[k])/h ) * fast_dnorm( (gx[i] - x[k])/h);
         }

         z[ j*n + i ]  = sum  / (nxy * h * h);


      }
   }

   return;


}


SEXP kde2d_wrap_( SEXP in_n, SEXP in_x, SEXP in_y,
                  SEXP in_h, SEXP in_lims){


   int nxy, n;
   double *x, *y, *lims, h;

   const char *names[] = {"x", "y", "z", ""};
   SEXP sgx, sgy, sz;
   double *gx, *gy, *z;

   nxy = (LENGTH(in_x) < LENGTH(in_y)) ? LENGTH(in_x) : LENGTH(in_y);
   n = INTEGER(in_n)[0];


   if (n < 2 || nxy < 2) return(R_NilValue);
   x = REAL(in_x);
   y = REAL(in_y);

   if (LENGTH(in_lims) < 4) return(R_NilValue);
   lims = REAL(in_lims);

   h = REAL(in_h)[0];



   PROTECT(sgx = allocVector(REALSXP, n));
   gx = REAL(sgx);

   PROTECT(sgy = allocVector(REALSXP, n));
   gy = REAL(sgy);

   PROTECT(sz = allocMatrix(REALSXP, n, n ));
   z = REAL(sz);

   kde2d( nxy, n, x, y, h, lims, gx, gy, z);


   SEXP res = PROTECT(mkNamed(VECSXP, names));
   SET_VECTOR_ELT(res, 0, sgx);
   SET_VECTOR_ELT(res, 1, sgy);
   SET_VECTOR_ELT(res, 2, sz);

   UNPROTECT(4);
   return res;

}


