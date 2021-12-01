#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>              


struct ks_score{
   double score;
   double pvalue;
};

static
void swap(int* a, int* b) 
{ 
    int t = *a; 
    *a = *b; 
    *b = t; 
} 

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

static
void match_subarray( int N, int* array1, int* array2, int* res)
{
    int* ptr = array1;
    int* ptr2 = array2;

    int i=0;

    while(ptr < array1 + N) {
        if(*ptr == *ptr2) {
            ptr++;
            /* printf("%d\n", ptr2 - array2); */
            res[i++] = ptr2 -array2;
        }
        ptr2++;
    }
   return;
}


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

static
double ks_test_simple( int n_x, int n_y, int* y)
{

  double N = 1.0 * n_x * n_y/(n_x + n_y);
  int *w;
  int *order;
  int i, j, temp;
  double *z, minz;

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
  for (i = 0; i < n_x + n_y -1; i++){
     if ( w[i] != w[i+1] ) minz = (z[i] < minz) ? z[i] : minz;
  }   
  minz = (z[n_x+n_y-1] < minz)? z[n_x+n_y-1] : minz;

  free(w);
  free(order);
  free(z);

  //printf("minz=%12.9f\n", minz);

  return (exp(-2.0 * N * minz*minz));
}

static
struct ks_score ks_genscore( int n_x, int n_y, int* y)
{
   struct ks_score res;

   double N = 1.0 * n_x * n_y/(n_x + n_y);
   double hit = 1.0/n_y;
   double mis = 1.0/n_x;


   int *y2;
   int y2_N;

   int* D;
   int* y_match;
   int i, i_save=0;
   double zmax = -999999.9;

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
   res.pvalue =  ks_test_simple( n_x, n_y, y);

   free(y2);
   free(y_match);
   free(D); 

   return(res);
}


SEXP ks_genescore_c(SEXP n_x, SEXP y)
{
    int N = LENGTH(y);
    int nx = asInteger(n_x);
    struct ks_score ks;
    SEXP ans;
    int *yarray;
    y = PROTECT(coerceVector(y, INTSXP));
    yarray = INTEGER(y);
    ks = ks_genscore(nx, N, yarray);
    PROTECT(ans = allocVector(REALSXP, 2 ));
    REAL(ans)[0] = ks.score;
    REAL(ans)[1] = ks.pvalue;
    UNPROTECT(2);
    return (ans);
}

SEXP ks_genescore_mat_(SEXP mat)
{
    struct ks_score ks;
    SEXP Rdim, ans;
    double *ymat, *rans;
    int N, ncol, nrow, i, j;
    int *yarray;

    ymat = REAL(mat);
    Rdim = getAttrib(mat, R_DimSymbol);
    nrow = INTEGER(Rdim)[0];
    ncol = INTEGER(Rdim)[1];
    yarray = (int*) malloc( ncol * sizeof(int) );
    PROTECT(ans = allocMatrix(REALSXP, 2, nrow ));
    rans = REAL(ans);

    //printf("number of rows = %d, number of columns = %d\n", nrow, ncol);
    for (i = 0; i < nrow; i++)
    {
      //printf("i = %d\n", i);
      for (j = 0, N=0; j < ncol; j++){
        if (ymat[ j * nrow + i] >0.1 ) { 
            //printf("%d, ", j);
            yarray[N++] = j+1;
         }
       }
       ks = ks_genscore(ncol, N, yarray);
       rans[i * 2] = ks.score;
       rans[i * 2 + 1] = ks.pvalue;
       //printf( "\n %12.8f %12.8f\n", ks.score, ks.pvalue);
    }
    free(yarray);

    //REAL(ans)[0] = ks.score;
    //REAL(ans)[1] = ks.pvalue;
    UNPROTECT(1);
    return (ans);
}


