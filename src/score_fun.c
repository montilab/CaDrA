#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>              



// Define a structure, so we could pass it back to R
struct ks_score
{
  double score;
  double pvalue;
};


// -------------------------------------------- //
//
// swap
// Swap two elements in an array
//
// INPUT & OUTPUT:
//   int *a
//   int *b
// --------------------------------------------- //
static
  void swap(int* a, int* b) 
  { 
    int t = *a; 
    *a = *b; 
    *b = t; 
  } 


// -------------------------------------------- //
//
// cadra_pkstwo
// Simplified implementation of pkstwo function
// 
// INPUT:
//    x
//    tol - tolerance (default: = 1e-6)
// 
// OUTPUT:
//    \sum_{-\infty}^\infty (-1)^k e^{-2k^2x^2}
// -------------------------------------------- //

static double
  cadra_pkstwo(double x, double tol)
  {
    /* This is a simple case of the original pkstwo function for n = 1*/
    double new, old, s, w, z;
    int k, k_max;
    double res;
    
    k_max = (int) sqrt(2 - log(tol));
    
    if(x < 1) {
      z = - (M_PI_2 * M_PI_4) / (x * x);
      w = log(x);
      s = 0;
      for(k = 1; k < k_max; k += 2) {
        s += exp(k * k * z - w);
      }
      res = s / M_1_SQRT_2PI;
    }
    else {
      z = -2 * x * x;
      s = -1;
      k = 1;
      old = 0;
      new = 1;
      while(fabs(old - new) > tol) {
        old = new;
        new += 2 * s * exp(z * k * k);
        s *= -1;
        k++;
      }
      res = new;
    }
    return(res);
  }


// -------------------------------------------- //
// 
// partition
// Swap elements if the left element is larger than the right element
// 
// INPUT:
//    arr - an array to be sorted
//    order - an array with the order of the eleents
//    low - left position
//    high - right position
//
// OUTPUT:
//    arr - an array to be sorted
//    order - an array with the order of the eleents
// -------------------------------------------- //

static
  int partition (int *arr, int *order, int low, int high) 
  { 
    int pivot = arr[high]; // pivot 
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
    int j; 
    
    for (j = low; j <= high - 1; j++) 
    { 
      // If current element is smaller than the pivot 
      if (arr[j] < pivot) 
      { 
        i++; // increment index of smaller element 
        swap(&arr[i], &arr[j]); 
        swap(&order[i], &order[j]);
      } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    swap(&order[i+1], &order[high]);
    return (i + 1); 
  } 


// -------------------------------------------- //
// 
// quickSort
// Quick Sort of an integer array
//
// INPUT:
//    arr - input array
//    order - initial order
//    low - left position within array
//    high - right position within array
//
// OUTPUT:
//    array - sorted array
//    order - order after the sorting
// -------------------------------------------- //
static
  void quickSort(int *arr, int *order, int low, int high) 
  { 
    if (low < high) 
    { 
      /* pi is partitioning index, arr[p] is now 
       at right place */
      int pi = partition(arr, order, low, high); 
      
      // Separately sort elements before 
      // partition and after partition 
      quickSort(arr, order, low, pi - 1); 
      quickSort(arr, order, pi + 1, high); 
    } 
  } 


// ---------------------------------------------- //
//                                                //
// match_subarray
//
//  Find position of array2 within array 1
//
// INPUT:
//    N - number of elements in array 1
//    array1 
//    array2
//
// OUTPUT:
//   res - resulting array with positions
// ------------------------------------------------ //
static
  void match_subarray( int N, int* array1, int* array2, int* res)
  {
    int* ptr = array1;
    int* ptr2 = array2;
    
    int i=0;
    
    while(ptr < array1 + N) {
      if(*ptr == *ptr2) {
        ptr++;
        res[i++] = ptr2 -array2;
      }
      ptr2++;
    }
    return;
  }


// ---------------------------------------------- //
//                                                //
// make_Y_array
//
// Prepare Y array for computing KS score
//
// INPUT:
//   N - number of elements in the input y array
//   y - array with integers
//
// OUTPUT:
//   y2 - output array
// ----------------------------------------------- //

static
  int make_Y_array( int N, int* y, int* y2 ) 
  {
    int i, j=2;
    
    y2[0] = y[0] - 1;
    y2[1] = y[0];
    for ( i = 1; i < N; i++)
    {
      if(y[i] - 1 == y2[j-1])
        y2[j++] = y[i];
      else
      {
        y2[j++] = y[i] - 1;
        y2[j++] = y[i];
      }   
    }
    
    return(j);
  } 

// ---------------------------------------------- //
//                                                //
// ks_test_simple
//
// Simple implementation of Kolmogorov-Smirnov test
// score computation
//
//
static
  double ks_test_simple( int n_x, int n_y, int* y, int alt)
  {
    
    double N = 1.0 * n_x * n_y/(n_x + n_y);
    int *w;
    int *order;
    int i, temp;
    double *z, minz, maxz, maxabsz;
    double result;
    
    w = (int*) malloc( (n_x + n_y) * sizeof(int) );
    order = (int*) malloc( (n_x + n_y) * sizeof(int) );
    z = (double*) malloc( (n_x + n_y) * sizeof(double) );
    
    for (i = 0; i < n_x; i++) {
      w[i] = i+1;
      order[i] = i+1;
    }
    for (i = n_x; i < n_x + n_y; i++){
      w[i] = y[i - n_x];
      order[i] = i+1;
    }
    
    // Calculate order of the elements in w array
    quickSort( w, order, 0, n_x + n_y - 1);
    
    for (i = 1; i < n_x + n_y; i++){
      if (w[i] == w[i-1]) {
        temp=order[i];
        order[i]=order[i-1];
        order[i-1]=temp;
      }
    }
    
    if ( order[0] <= n_x )z[0] = 1.0/n_x; else z[0] = -1.0/n_y;
    for (i = 1; i < n_x + n_y; i++)
    {
      if ( order[i] <= n_x ) z[i] = z[i-1] + 1.0/n_x; else z[i] = z[i-1] - 1.0/n_y;
    }
    
    
    minz =  1.0 * n_x + n_y;
    maxz = -999999.;
    maxabsz=-1;
    for (i = 0; i < n_x + n_y -1; i++){
      if ( w[i] != w[i+1] ) {
        minz = (z[i] < minz) ? z[i] : minz;
        maxz = (z[i] > maxz) ? z[i] : maxz;
        maxabsz = ( fabs(z[i]) > maxabsz ) ? fabs(z[i]) : maxabsz;
      }
    }
    minz = (z[n_x+n_y-1] < minz)? z[n_x+n_y-1] : minz;
    maxz = (z[n_x+n_y-1] > maxz)? z[n_x+n_y-1] : maxz;
    maxabsz = ( fabs(z[n_x+n_y-1]) > maxabsz) ? fabs(z[n_x+n_y-1]) : maxabsz;
    
    free(w);
    free(order);
    free(z);
    
    //printf("minz=%12.9f\n", minz);
    switch(alt){
    case 1:
      result = exp(-2.0 * N * minz*minz);
      break;
    case -1:
      result = exp(-2.0 * N * maxz*maxz);
      break;
    case 0:
      result = 1.0 - cadra_pkstwo(sqrt(1.0*N) * maxabsz, 1E-06);
      break;
    default:
      result = exp(-2.0 *N * minz*minz);
    }
    return (result);
  }

// ---------------------------------------------- //
//                                                //
// ks_genscore
//
// INPUT:
// int n_x
// int n_y - number of values in y array
// int* y  
// int n_w - number of elements in weight 
// int* weight - an integer array of weights
//
// OUTPUT:
//   ks_score structure: ks score and p-value
//
// int alt:   0: "two.sided", +1: "less"; -1: "greater"
// ----------------------------------------------------
static
  struct ks_score 
  ks_genscore( int n_x, int n_y, int* y, int n_w, double* weight, int alt)
  /* +1 = left-sided, 0 - two-sided, -1 = right-sided */
  {
    struct ks_score res;
    
    //double N = 1.0 * n_x * n_y/(n_x + n_y);
    double hit = 1.0/n_y;
    double mis = 1.0/n_x;
    
    
    int *y2;
    int y2_N;
    
    double *Phit, *Pmis;
    
    int* D;
    int* y_match;
    int i, i_save=0;
    double zmax = -999999.9;
    
    
    if (n_w == n_x)
    {
      Pmis = (double*) malloc (n_x * sizeof(double) );
      Phit = (double*) malloc (n_x * sizeof(double) );
      
      for (i = 0; i < n_x; i++) {
        Pmis[i] = 1.;
        Phit[i] = 0.;
      }


      for (i = 0; i < n_y; i++){
        Pmis[ y[i] -1 ] = 0.0;
        Phit[ y[i] -1 ] = (weight[ y[i] -1 ] >0) ? weight[ y[i] -1 ]: -weight[ y[i] -1 ];
      }
      for (i = 1; i < n_x; i++) {
        Pmis[i] += Pmis[i-1];
        Phit[i] += Phit[i-1];

      }
      for (i = 0; i < n_x; i++){
        
        if ( fabs( Pmis[i]/(n_x - n_y) - Phit[i]/Phit[ n_x -1]  ) > zmax) {
          zmax = fabs( Pmis[i]/(n_x - n_y) - Phit[i]/Phit[ n_x -1 ] ) ;
        }
      }
      res.score = zmax;
 
      free(Pmis);
      free(Phit);
      
    } else
    {
      
      y2 = (int*) malloc( 2 * n_y * sizeof(int) );
      y_match = (int*) malloc( n_y * sizeof(int) );
      
      y2_N = make_Y_array (n_y, y, y2);
      match_subarray( n_y, y, y2, y_match);
      
      D = (int*) malloc ( y2_N * sizeof(int) );
      for (i = 0; i < y2_N; i++) D[i] = 0;
      for (i = 0; i < n_y; i++) D[ y_match[i] ] = i + 1;
      for (i = 1; i < y2_N; i++) if (D[i] == 0) D[i] = D[i-1];
      for (i = 0; i < y2_N; i++) {
        if( fabs(D[i] * hit - y2[i] * mis) > zmax ) {
          zmax = fabs(D[i] * hit - y2[i]*mis);
          i_save = i;
        }
      }
      
      res.score =D[i_save] * hit - y2[i_save]*mis ;
      
      free(y2);
      free(y_match);
      free(D);
    }
    
    res.pvalue =  ks_test_simple( n_x, n_y, y, alt);
    
    return(res);
  }

// ---------------------------------------------------------------------- //
//                                                                        //
//  Wrapper for ks_genescore function                                     //
//
// ---------------------------------------------------------------------- //
SEXP ks_genescore_wrap_(SEXP in_n_x, SEXP in_y, SEXP in_w, SEXP alternative)
{
  struct ks_score ks;
  int *y;
  double *w;
  int alt;
  int n_x, n_y, n_w;
  double *rans;
  SEXP ans;

  n_x = INTEGER(in_n_x)[0];
  n_y = LENGTH(in_y);
  n_w = LENGTH(in_w);
  if (LENGTH(in_w)) w = REAL(in_w);
  if ( n_y < 1 || n_x < 1 ) return(R_NilValue);
  y = INTEGER(in_y);
  
  if (LENGTH(alternative) > 0) 
  {
    alt= INTEGER(alternative)[0]; /* +1 = left-sided, 0 = double sided, -1 = right sided */
  } else 
  {
    alt = 1;
  }

  PROTECT(ans = allocVector(REALSXP, 2 ));
  rans = REAL(ans);

 
  ks = ks_genscore( n_x, n_y, y, n_w, w, alt);
  rans[0] = ks.score;
  rans[1] = ks.pvalue;

  UNPROTECT(1);
  return (ans);


  }



// ---------------------------------------------------------------------- //
//                                                                        //
//  Compute directional KS scores for each row of a given binary matrix   //
//
// INPUT:
//  mat - matrix of binary features to compute row-wise scores for 
//        based on the  Kolmogorov-Smirnov test
//  alt - an integer value specifying the alternative hypothesis, 
//        must be one of 0: "two.sided", +1: "less"; -1: "greater". 
//  w   - a vector of weights to use if performing a weighted-KS test. 
//        Default is NULL.
//
// OUTPUT:
// A 2 x N matrix (N = number of rows in matrix) with 
//         the first and second rows corresponding to 
//         the KS scores and p-values, respectively 
// 
// ------------------------------------------------------------------------ //
SEXP ks_genescore_mat_(SEXP mat, SEXP w, SEXP alternative)
{
  struct ks_score ks;
  double *ymat, *rans;
  double *weight;
  int N, ncol, nrow, nw, i, j;
  int *yarray;
  int alt;
  SEXP Rdim, ans;
  
  ymat = REAL(mat);
  
  Rdim = getAttrib(mat, R_DimSymbol);
  nrow = INTEGER(Rdim)[0];
  ncol = INTEGER(Rdim)[1];
  nw = LENGTH(w);
  //printf("length of w: %d\n", nw);
  
  if (LENGTH(alternative) > 0) 
  {
    alt= INTEGER(alternative)[0]; /* +1 = left-sided, 0 = double sided, -1 = right sided */
  } else 
  {
    alt = 1;
  }
  
  if (nw > 0 && nw != ncol) 
  {
    return(R_NilValue);
  }
  
  if (nw > 0) weight = REAL(w);
  
  yarray = (int*) malloc( ncol * sizeof(int) );
  PROTECT(ans = allocMatrix(REALSXP, 2, nrow ));
  rans = REAL(ans);
  
  for (i = 0; i < nrow; i++)
  {
    
    for (j = 0, N = 0; j < ncol; j++)
    {
      if (ymat[ j * nrow + i] >0.1 ) { 
        yarray[N++] = j+1;
      }
    }
    ks = ks_genscore(ncol, N, yarray, nw, weight, alt);
    rans[i * 2] = ks.score;
    rans[i * 2 + 1] = ks.pvalue;
  }
  free(yarray);
  
  UNPROTECT(1);
  return (ans);
}


