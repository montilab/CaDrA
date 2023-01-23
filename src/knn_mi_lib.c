#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

extern int mutual_inf_cc_vec( double *input_x,  double *input_y,  int n_elems,   int k, double *mi) ;
extern int mutual_inf_cd_vec( double *input_x,  int *input_y,  int n_elems,  int k, double *mi ) ;
extern int cond_mutual_inf_vec( double *input_x,   double *input_y,  double *input_z,  int n_elems,  int k, double *mi) ;

/* Functions defined here:
 *  _mutual_inf_cc_1d - mutual information MI(x;y) continuous-continuous case. x and y are vectors of length N.
 *  _mutual_inf_cc_2d - MI(x;y) c-c case where y is an NxM matrix. Returns an MI vector of length M.
 *  _mutual_inf_cd_1d - MI(x;y) continuous-discrete case. x and y are vectors of length N. y is an integer vector of
 *      discrete values. Accepts a flag to internally use the continuous-continuous algorithm instead of a faster c-d algorithm.
 *  _mutual_inf_cd_2d - mutual information MI(x;y) c-d case where y is an NxM matrix.
 *      Also accepts the flag to switch the algorithm.
 *  _cond_mutual_inf_ccc_1d - conditional mutual information CMI(x;y|z), where x,y,z are continuous values and are
 *      vectors of length N.
 *  _cond_mutual_inf_ccc_2d - CMI(x;y|z), where x,y,z are continuous values. x is of length N, y and z are matrices of
 *      size NxM.
 *  _cond_mutual_inf_cdd_1d - conditional mutual information CMI(x;y|z), where x is continuous. y and z are
 *      integer vectors of discrete values. All are vectors of length N.
 *  _cond_mutual_inf_ccc_2d - CMI(x;y|z), where x is continuous and both y and z are discrete. x is of length N,
 *      y and z are integer matrices of size NxM.
 *
 *  The 1d functions all return a single numeric value. The 2d functions return a numeric vector of length M.
 *
 *  The R code that calls these functions will implement checks to make sure that vector/matrix sizes are
 *  correct.
 *
 */




SEXP _mutual_inf_cc_1d( SEXP r_input_x, SEXP r_input_y, SEXP k) {
    /* R C wrapper for:
     * int mutual_inf_cc_vec(double *input_x, double *input_y, int n_elems, double *mi,  int k)
     *
     * r_input_x and r_input_y are vectors of the same length. k is the number of nearest neighbors.
     * The mutual information is returned as a numeric vector of length 1.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int n_elems, k_value ;

    n_elems = LENGTH(r_input_x);
    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, 1));

    p_mi = REAL(mi);

    mutual_inf_cc_vec(REAL(r_input_x), REAL(r_input_y), n_elems,k_value, p_mi);

    /* Lift R garbage protection */
    UNPROTECT(1);

    return( mi );
}

SEXP _mutual_inf_cc_2d( SEXP r_input_x, SEXP r_input_y, SEXP k) {
    /* R C wrapper for:
     * int mutual_inf_cc_vec(double *input_x, double *input_y, int n_elems, double *mi,  int k)
     *  r_input_x - the vector of continuous data. size N.
     *  r_input_y - input matrix. size NxM
     *  k - number of nearest neighbors.
     *
     *  The mutual information is returned as a vector of length M.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi, *p_y, *p_x ;
    int y_nrows, y_ncols, k_value ;
    int i ;

    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* Get the dimensions of r_input_y */
    y_ncols = Rf_ncols(r_input_y);
    y_nrows = Rf_nrows(r_input_y);

    /* Verify we've loaded the dimensions correctly. */
    /* printf("%d %d %d %d\n", x_elems, k_value, y_ncols, y_nrows); */

    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, y_ncols));

    /* get pointers */
    p_mi = REAL(mi);
    p_y = REAL(r_input_y) ;
    p_x = REAL(r_input_x) ;
    for (i = 0 ; i < y_ncols ; i++) {
        mutual_inf_cc_vec(p_x, &p_y[y_nrows * i], y_nrows, k_value, &p_mi[i]);
    }
    /* Lift R garbage protection */
    UNPROTECT(1);

    return( mi );
}

SEXP _mutual_inf_cd_1d( SEXP r_input_x, SEXP r_input_y, SEXP k, SEXP use_cc) {
    /* R C wrapper for:
     * int mutual_inf_cc_vec( double *input_x,  int *input_y,  int n_elems,  int k, double *mi ) ;
     *
     * r_input_x and r_input_y are vectors of the same length. k is the number of nearest neighbors.
     * r_input_y is a vector of integers.
     *
     * The mutual information is returned as a vector of length 1.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int n_elems, k_value ;
    int prot_ctr = 0 ;
    int i ;
    int *p_y ;
    SEXP dbl_y ;
    double *p_dbl_y, *p_x ;
    int use_cc_val = asLogical(use_cc);

    n_elems = LENGTH(r_input_x);
    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, 1));
    prot_ctr++ ;

    p_mi = REAL(mi);

    p_y = INTEGER(r_input_y) ;
    p_x = REAL(r_input_x) ;
    /* There is an implementation provided for the specific continuous-discrete case
    * for calculating MI.  It returns the same result as the continuous-continuous
    * implementation (while running 2x faster) for larger sets of data (N > 100?)
    * but sometimes returns a negative value for very small sets.  The R code will
    * be set up to use the c-c calculation by default but the flag can be used
    * to switch it to the c-d algorithm if desired. */
    if (!use_cc_val) {
        mutual_inf_cd_vec(p_x, p_y, n_elems, k_value, p_mi);
    } else {
        dbl_y = PROTECT(allocVector(REALSXP, n_elems));
        prot_ctr++ ;
        p_dbl_y = REAL(dbl_y);
        for (i = 0 ; i < n_elems ; i++) {
            p_dbl_y[i] = (double) p_y[i] ;
        }
        mutual_inf_cc_vec(p_x, p_dbl_y, n_elems,k_value, p_mi);
    }

    /* Lift R garbage protection */
    UNPROTECT(prot_ctr);

    return( mi );
}

SEXP _mutual_inf_cd_2d( SEXP r_input_x, SEXP r_input_y, SEXP k, SEXP use_cc) {
    /* R C wrapper for:
     * int mutual_inf_cd(double *input_x, int x_elems, int *input_y, int y_nrows, int y_ncols, int k, double *mi) ;
     *  r_input_x - the vector of continuous data. size N.
     *  r_input_y - input matrix. size NxM
     *  k - number of nearest neighbors.
     *  use_cc - logical flag, TRUE falls back to the continuous-continuous algorithm.
     *  The mutual information is returned as a vector of length M.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int y_nrows, y_ncols, k_value ;
    int *p_y ;
    SEXP dbl_y ;
    double *p_dbl_y, *p_x ;
    int prot_ctr = 0 ;
    int use_cc_val = asLogical(use_cc);
    int i,j ;

    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* Get the dimensions of r_input_y */
    y_ncols = Rf_ncols(r_input_y);
    y_nrows = Rf_nrows(r_input_y);

    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, y_ncols));
    prot_ctr++ ;

    p_mi = REAL(mi);
    p_x = REAL(r_input_x);
    p_y = INTEGER(r_input_y) ;

    /* As with the 1D case try using the c-d algorithm unless a negative MI is returned.
     * If it is...re-do the entire calculation with the c-c algorithm, doing type conversions
     * of the integer columns as we go. */
    if (!use_cc_val) {
        for (j = 0 ; j < y_ncols ; j++) {
            mutual_inf_cd_vec(p_x, &p_y[j * y_nrows], y_nrows, k_value, &p_mi[j]);
        }
    } else {
        dbl_y = PROTECT(allocVector(REALSXP, y_nrows));
        prot_ctr++ ;
        p_dbl_y = REAL(dbl_y) ;
        for (j = 0 ; j < y_ncols ; j++) {
            for (i = 0 ; i < y_nrows ; i++) {
                p_dbl_y[i] = (double) p_y[j * y_nrows + i] ;
            }
            mutual_inf_cc_vec(p_x, p_dbl_y, y_nrows, k_value, &p_mi[j]);
        }
    }

    /* Lift R garbage protection */
    UNPROTECT(prot_ctr);

    return( mi );
}

SEXP cond_mutual_inf_ccc_1d_( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
    /* R C wrapper for:
     * int cond_mutual_inf_vec( double *input_x,   double *input_y,  double *input_z,  int n_elems,  int k, double *mi) ;
     *
     * r_input_x, r_input_z, and r_input_y are vectors of the same length. k is the number of nearest neighbors.
     * The mutual information is returned as a vector of length 1.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int n_elems, k_value ;

    n_elems = LENGTH(r_input_x) ;
    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, 1)) ;

    p_mi = REAL(mi) ;

    cond_mutual_inf_vec(REAL(r_input_x), REAL(r_input_y), REAL(r_input_z), n_elems, k_value, p_mi) ;

    /* Lift R garbage protection */
    UNPROTECT(1) ;

    return( mi ) ;
}

SEXP _cond_mutual_inf_cdd_1d( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
    /* R C wrapper for:
     * int cond_mutual_inf_vec( double *input_x,   double *input_y,  double *input_z,  int n_elems,  int k, double *mi) ;
     *
     * r_input_x, r_input_z, and r_input_y are vectors of the same length. k is the number of nearest neighbors.
     *
     * Here r_input_y and r_input_z are discrete vectors so they're integers.
     *
     * The mutual information is returned as a vector of length 1.
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int n_elems, k_value ;
    int i ;
    int *p_y, *p_z ;
    SEXP dbl_y, dbl_z ;
    double *p_dbl_y, *p_dbl_z ;

    n_elems = LENGTH(r_input_x);
    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, 1));

    p_mi = REAL(mi);

    /* There is only 1 conditional mutual information function and it doesn't deal with integers.
     * Allocate storage for the y & z inputs as double arrays, convert them, and pass those. */
    dbl_y = PROTECT(allocVector(REALSXP, n_elems));
    dbl_z = PROTECT(allocVector(REALSXP, n_elems));
    p_dbl_y = REAL(dbl_y) ;
    p_dbl_z = REAL(dbl_z) ;

    p_y = INTEGER(r_input_y) ;
    p_z = INTEGER(r_input_z) ;

    for (i = 0 ; i < n_elems ; i++) {
        p_dbl_y[i] = (double) p_y[i] ;
        p_dbl_z[i] = (double) p_z[i] ;
    }

    cond_mutual_inf_vec(REAL(r_input_x), p_dbl_y, p_dbl_z, n_elems, k_value, p_mi);

    /* Lift R garbage protection */
    UNPROTECT(3);

    return( mi );
}


SEXP _cond_mutual_inf_ccc_2d( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
    /* R C wrapper for:
     * int cond_mutual_inf( double *input_x,  int x_elems,  double *input_y,  double *input_z,  int nrows,  int ncols,  int k, double *mi) ;
     *  r_input_x - the vector of continuous data. size N.
     *  r_input_y - input matrix. size NxM
     *  r_input_z - input matrix. size NxM
     *  k - number of nearest neighbors.
     *
     *  The mutual information is returned as a vector of length M.
     *
     *  in CaDrA  terms  input_score = x, input_mat=y, expression_score=z
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int  yz_nrows, yz_ncols, k_value ;
    int j ;
    double *p_x, *p_y, *p_z ;

    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* Get the dimensions of r_input_y */
    /* r_input_z must have the same dimensions, that should all be checked
     * on the R side. */
    yz_ncols = Rf_ncols(r_input_y);
    yz_nrows = Rf_nrows(r_input_y);

    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, yz_ncols));
    p_mi = REAL(mi);

    p_x = REAL(r_input_x) ;
    p_y = REAL(r_input_y) ;
    p_z = REAL(r_input_z) ;
    for (j = 0 ; j < yz_ncols ; j++) {
        cond_mutual_inf_vec(p_x, &p_y[j * yz_nrows], &p_z[j * yz_nrows], yz_nrows, k_value, &p_mi[j]);
    }

    /* Lift R garbage protection */
    UNPROTECT(1);

    return( mi );
}


SEXP _cond_mutual_inf_cdd_2d( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
    /* R C wrapper for:
     * int cond_mutual_inf( double *input_x,  int x_elems,  double *input_y,  double *input_z,  int nrows,  int ncols,  int k, double *mi) ;
     *  r_input_x - the vector of continuous data. size N.
     *  r_input_y - integer input matrix. size NxM
     *  r_input_z - integer input matrix. size NxM
     *  k - number of nearest neighbors.
     *
     *  The mutual information is returned as a vector of length M.
     *
     *  in CaDrA  terms  input_score = x, input_mat=y, expression_score=z
    */

    /* declare the output mutual information */
    SEXP mi;
    double *p_mi;
    int yz_nrows, yz_ncols, k_value ;
    double *p_dbl_y, *p_dbl_z ;
    SEXP dbl_y, dbl_z ;
    int i, j ;
    int *p_y, *p_z ;
    double *p_x ;

    /* In the R wrapper make sure this is passed as an integer. */
    k_value = INTEGER(k)[0] ;
    /* Get the dimensions of r_input_y */
    /* r_input_z must have the same dimensions, that should all be checked
     * on the R side. */
    yz_ncols = Rf_ncols(r_input_y);
    yz_nrows = Rf_nrows(r_input_y);

    /* R memory allocation */
    mi = PROTECT(allocVector(REALSXP, yz_ncols));
    p_mi = REAL(mi);

    /* To avoid excess memory usage, loop here and call cond_mutual_inf_vec, not cond_mutual_inf,
     * while doing integer -> double conversions */
    dbl_y = PROTECT(allocVector(REALSXP, yz_nrows));
    dbl_z = PROTECT(allocVector(REALSXP, yz_nrows));
    p_dbl_y = REAL(dbl_y) ;
    p_dbl_z = REAL(dbl_z) ;
    p_y = INTEGER(r_input_y) ;
    p_z = INTEGER(r_input_z) ;
    p_x = REAL(r_input_x) ;
    for (j = 0 ; j < yz_ncols ; j++) {
        /* Convert the jth column of the input y & z to double */
        for (i = 0 ; i < yz_nrows ; i++) {
            p_dbl_y[i] = (double) p_y[i + j * yz_nrows] ;
            p_dbl_z[i] = (double) p_z[i + j * yz_nrows] ;
        }
        cond_mutual_inf_vec(p_x, p_dbl_y, p_dbl_z, yz_nrows, k_value, &p_mi[j]);
    }

    /* Lift R garbage protection */
    UNPROTECT(3);

    return( mi );
}


