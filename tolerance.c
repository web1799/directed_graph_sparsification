/* norm-2 of array x */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double tolerance(double *x, int n){
    int i;
    double beta, tol;
    beta = 0.0;
    for(i=0; i<n; i++){
        beta = beta + fabs(x[i]*x[i]);
    }
    tol = sqrt(beta);
    return tol;
}