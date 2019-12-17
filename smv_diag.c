 /* b = A*x
  * A: symatric matrix 
  * irow: row index
  * pcol: col pointer
  * val:  matrix val
  * diag: diag array of A
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void smv_diag(int *irow, int *pcol, double *val, double *x, double *b, int n){
    int i, j;
    for(i=0; i<n; i++){
        b[i] = 0.0;
        for(j=pcol[i]; j<pcol[i+1]; j++){
            b[i] += val[j]*x[irow[j]];
        }
    }
}