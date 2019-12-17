 /* b = A*x
  * A: lower triangular of matrix A
  * irow: row index
  * pcol: col pointer
  * val:  matrix val
  * diag: diag value of A
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void smv(int *irow, int *pcol, double *val, double *x, double *b, int n){
    int i, j;
    double sum;
    
    for(i=0; i<n; i++){
        b[i] = 0.0;
    }
    
    for(i=0; i<n; i++){
        sum  = val[pcol[i]]*x[i];
        for(j=pcol[i]+1; j<pcol[i+1]; j++){
            b[irow[j]] += val[j]*x[i];
            sum += val[j]*x[irow[j]];
        }
        b[i] += sum;
    }
}