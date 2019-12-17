#include <stdlib.h>
/*  
 * irow:    row index
 * pcol:    col pointer
 *  val:    matrix non-zero value
 * diag:    diag value of matrix
 */ 
typedef struct matrix{
    int *irow, *pcol;
    double *val, *diag;
} Matrix;

double tolerance(double *x, int n);
void smv(int *irow, int *pcol, double *val, double *x, double *b, int n);
void smv_diag(int *irow, int *pcol, double *val, double *x, double *b, int n);
/*double Norm(int *irow, int *pcol, double *val, double *x, double *b, int n);
*/
