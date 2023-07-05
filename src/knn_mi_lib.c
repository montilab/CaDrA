#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

extern int mutual_inf_cc_vec( double *input_x,  double *input_y,  int n_elems, int k, double *mi) ;
extern int mutual_inf_cd_vec( double *input_x,  int *input_y,  int n_elems,  int k, double *mi ) ;
extern int cond_mutual_inf_vec( double *input_x,   double *input_y,  double *input_z,  int n_elems,  int k, double *mi) ;

/* Functions defined here:
 *  _mutual_inf_cc - MI(x;y) c-c case where y is an NxM matrix. Returns an MI vector of length M.
 *  _mutual_inf_cd  - mutual information MI(x;y) c-d case where y is an NxM matrix.
 *  _cond_mutual_inf_ccc  - CMI(x;y|z), where x,y,z are continuous values. x is of length N, y and z are vectors
 *    of size N or matrices of size NxM.
 *  _cond_mutual_inf_cdd  - CMI(x;y|z), where x is continuous and both y and z are discrete. x is of length N,
 *      y and z are integer vectors of size N or integer matrices of size NxM.
 *
 *  The 1d functions all return a single numeric value. The 2d functions return a numeric vector of length M.
 *
 *  The R code that calls these functions will implement checks to make sure that vector/matrix sizes are
 *  correct.
 *
 */


SEXP _mutual_inf_cc( SEXP r_input_x, SEXP r_input_y, SEXP k) {
  /* R C wrapper for:
   * int mutual_inf_cc_vec(const double *input_x, const double *input_y, const int n_elems,  
   *                       const int k, double *mi)
   * r_input_x is the target vector. It is an R vector of length N.
   * r_input_y is the features matrix. It can be an R vector or length N or matrix of size MxN.
   * Assume that the size double-checking happened in R.
   * k is the number of nearest neighbors.
   * The mutual information is returned as a numeric vector of length 1.
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi, *p_y, *p_x, *p_buf_y ;
  int n_rows, n_cols,   k_value, i, j  ;
  
  /* Assume x is a vector - this is checked on the R side.*/
  n_cols = LENGTH(r_input_x) ;

  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Check if r_input_y is a vector or a matrix. */
  n_rows = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    n_rows = Rf_nrows(r_input_y);
  } 

  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, n_rows));
  
  /* get pointers to the inputs and outputs */
  p_mi = REAL(mi);
  p_y = REAL(r_input_y) ;
  p_x = REAL(r_input_x) ;
  
  /* If y is a matrix, allocate a buffer */
  if (n_rows > 1) {
    p_buf_y = (double *) malloc(n_cols * sizeof(double)) ;
  } else {
    /* Point the buffer to the y vector */
    p_buf_y = p_y ;
  }
  
  for (i = 0 ; i < n_rows ; i++) {
    /* Copy the buffer if there are more than 1 row */
    if (n_rows > 1) {
      for (j = 0 ; j < n_cols ; j++) {
        p_buf_y[j] = p_y[j * n_rows + i] ;
      }
    }
    mutual_inf_cc_vec(p_x, p_buf_y, n_cols, k_value, &p_mi[i]);
  }
  
  /* free the buffer if needed */
  if (n_rows > 1) {
    free(p_buf_y);  
  } 
  
  /* Lift R garbage protection */
  UNPROTECT(1);  
  /* Return the mi vector */
  return( mi );
}

SEXP _mutual_inf_cd( SEXP r_input_x, SEXP r_input_y, SEXP k) {
  /* R C wrapper for:
   * int mutual_inf_cd_vec(const double *input_x, const int *input_y, const int n_elems, 
   *                       const int k, double *mi )
   *  r_input_x - the vector of continuous data. size N.
   *  r_input_y - input matrix. size NxM
   *  k - number of nearest neighbors.
   *  The mutual information is returned as a vector of length M.
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi;
  int n_rows, n_cols, k_value ;
  int *p_y, *p_buf_y ;
  double *p_x ;
  int i, j ;
  
  /* Assume x is a vector - this is checked on the R side.*/
  n_cols = LENGTH(r_input_x) ;
  
  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Check if r_input_y is a vector or a matrix. */
  n_rows = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    n_rows = Rf_nrows(r_input_y);
  } 
  
  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, n_rows));
  
  /* get pointers to the inputs and outputs */
  p_mi = REAL(mi);
  p_y = INTEGER(r_input_y) ;
  p_x = REAL(r_input_x) ;
  
  /* If y is a matrix, allocate a buffer */
  if (n_rows > 1) {
    p_buf_y = (int *) malloc(n_cols * sizeof(int)) ;
  } else {
    /* Point the buffer to the y vector */
    p_buf_y = p_y ;
  }
  
  for (i = 0 ; i < n_rows ; i++) {
    /* Copy the buffer if there are more than 1 row */
    if (n_rows > 1) {
      for (j = 0 ; j < n_cols ; j++) {
        p_buf_y[j] = p_y[j * n_rows + i] ;
      }
    }
    mutual_inf_cd_vec(p_x, p_buf_y, n_cols, k_value, &p_mi[i]);
  }
  
  /* free the buffer if needed */
  if (n_rows > 1) {
    free(p_buf_y);  
  } 
  
  /* Lift R garbage protection */
  UNPROTECT(1);  
  /* Return the mi vector */
  return( mi );
}


SEXP _cond_mutual_inf( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k, SEXP case_) {
  /* R C wrapper for:
   * int cond_mutual_inf( double *input_x,  int x_elems,  double *input_y,  double *input_z,  int nrows,  int ncols,  int k, double *mi) ;
   *  r_input_x - the vector of continuous data. size N.
   *  r_input_y - integer vector or matrix. size N or size NxM
   *  r_input_z - integer vector or matrix. size N or size NxM
   *  k - number of nearest neighbors.
   *
   *  The mutual information is returned as a vector of length M.
   *
   *  in CaDrA  terms  input_score = x, input_mat=y, expression_score=z
   */
  
  /* declare the output mutual information */
  SEXP mi ;
  double *p_mi;
  int len_x, num_rows, k_value, case_value ;
  double *p_buf_y, *p_buf_z ;
  int i, j ;
  double *p_x, *p_y, *p_z ;
  char err_msg[256] ;
  
  /* In the R wrapper make sure k is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  /* Length of x vector */
  len_x =  Rf_length(r_input_x);
  
  /* 4 cases: 
   *  0: y and z are vectors, 
   *  1: y is a vector, z is a matrix
   *  2: y is a matrix, z if a vector
   *  3: y and z are both matrices
   *  
   *  With matrices the number of columns matches the
   *  length of x. As they come from R they're stored in
   *  column-major format, so they'll be copied row-by-row
   *  to temporary buffers where needed. 
   * In the R wrapper make sure case_ is passed as an integer. */
  case_value = INTEGER(case_)[0] ;
  /* Check to make sure this is in a valid range */
  if (case_value < 0 || case_value > 3) {
    sprintf(err_msg, "value of case argument must be in the range 0-3. Value given: %d", case_value) ;
    error(err_msg) ;
  };

  /* Get C pointers to the R data */
  p_y = REAL(r_input_y) ;
  p_z = REAL(r_input_z) ;
  p_x = REAL(r_input_x) ;

  /* and....go! */
  switch (case_value) {
    case 0:
      /* R memory allocation for the return vector */
      mi = PROTECT(allocVector(REALSXP, 1));
      p_mi = REAL(mi);
      cond_mutual_inf_vec(p_x, p_y, p_z, len_x, k_value, &p_mi[0]) ;
      break ;
    case 1:
      /* Get the number of rows in the Z matrix */    
      num_rows = Rf_nrows(r_input_z);
      /* R memory allocation for the return vector */
      mi = PROTECT(allocVector(REALSXP, num_rows));
      p_mi = REAL(mi);
      /* Get a row buffer */
      p_buf_z = (double *) malloc(len_x * sizeof(double)) ;
      /* Loop over each row of Z */
      for (i = 0 ; i < num_rows ; i++) {
        for (j = 0 ; j < len_x ; j++) {
          /* Move along the row when the matrix is column-major */
          p_buf_z[j] = p_z[j * num_rows + i] ;
        }
        cond_mutual_inf_vec(p_x, p_y, p_buf_z, len_x, k_value, &p_mi[i]) ; 
      }
      free(p_buf_z) ;
      break ;
    case 2:
      /* Get the number of rows in the Y matrix */    
      num_rows = Rf_nrows(r_input_y);
      /* R memory allocation for the return vector */
      mi = PROTECT(allocVector(REALSXP, num_rows));
      p_mi = REAL(mi);
      /* Get a row buffer */
      p_buf_y = (double *) malloc(len_x * sizeof(double)) ;
      /* Loop over each row of Z */
      for (i = 0 ; i < num_rows ; i++) {
        for (j = 0 ; j < len_x ; j++) {
          /* Move along the row when the matrix is column-major */
          p_buf_y[j] = p_y[j * num_rows + i] ;
        }
        cond_mutual_inf_vec(p_x, p_buf_y, p_z, len_x, k_value, &p_mi[i]) ; 
      }
      free(p_buf_y);
      break ;
    case 3:
      /* Get the number of rows in one of the matrices */    
      num_rows = Rf_nrows(r_input_y);
      
      /* R memory allocation for the return vector */
      mi = PROTECT(allocVector(REALSXP, num_rows));
      p_mi = REAL(mi);
      /* Get row buffers */
      p_buf_y = (double *) malloc(len_x * sizeof(double)) ;
      p_buf_z = (double *) malloc(len_x * sizeof(double)) ;
      /* Loop over each row of Z */
      for (i = 0 ; i < num_rows ; i++) {
        for (j = 0 ; j < len_x ; j++) {
          /* Move along the row when the matrix is column-major */
          p_buf_y[j] = p_y[j * num_rows + i] ;
          p_buf_z[j] = p_z[j * num_rows + i] ;
        }
        cond_mutual_inf_vec(p_x, p_buf_y, p_buf_z, len_x, k_value, &p_mi[i]) ; 
      }
      free(p_buf_y) ;
      free(p_buf_z) ;
      break ;    
  }
  UNPROTECT(1);
  return( mi );
} 

