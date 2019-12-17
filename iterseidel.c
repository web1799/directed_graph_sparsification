/*
 * iterjacobi: iterative method using Gauss-Seidel for calculating Ax = b
 * Since the for loop in matlab is slow, so we use the c function to implement it.
*/

#include "mex.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "solver.h"

#define x_out plhs[0]
#define C prhs[0]
#define x0_in prhs[1]
#define b_in prhs[2]
#define iter_in prhs[3]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int n;
    int i, j;
    Matrix A;
    double *x, *b, *x0;
    double *array;
    double *y;
	int iter, it;
    double tol, TOL, tol0;
    double *beta, gamma;
    int  flag;
    clock_t tstart, tend;
    double telapsed;
 
    
    
    /* check input and output arguments*/
    if(nlhs != 1){
        mexErrMsgIdAndTxt("iterseidel:nlhs", "one output required.");
    }
    if(nrhs != 4){
        mexErrMsgIdAndTxt("iterseidel:nrhs", "four input required.");
    } 
    if(!mxIsStruct(prhs[0]) || mxGetNumberOfElements(prhs[0]) == 0){
        mexErrMsgIdAndTxt("iterseidel:NotStruct", "The first input must be a non-empty struct");
    }
    if(mxGetN(prhs[1]) != 1){
        mexErrMsgIdAndTxt("iterseidel:NotColVector", "Second input must be a column vector");
    }
    /*
    if(mxIsDouble(prhs[2] != 1)){
        mexErrMsgIdAndTxt("iterseidel:NotDouble", "Third input must be double type");
    }*/

    n = (int) mxGetScalar(mxGetField(C, 0, "n"));
    A.irow = (int *) mxGetPr(mxGetField(mxGetField(C, 0, "A"), 0, "irow"));
    A.pcol = (int *) mxGetPr(mxGetField(mxGetField(C, 0, "A"), 0, "pcol"));
    A.val = (double *) mxGetPr(mxGetField(mxGetField(C, 0, "A"), 0, "val"));        
    A.diag = (double *) mxGetPr(mxGetField(C, 0, "diag")); 
    
    x_out = mxCreateDoubleMatrix(n, 1, mxREAL); /* create rhs vector */
    x = (double *) mxGetPr(x_out);
    b = (double *) mxGetPr(b_in);
	x0 = (double*) mxGetPr(x0_in);
    iter = (int) mxGetScalar(iter_in); 
    array = (double*)mxMalloc(n*sizeof(double));
    beta = (double*)mxMalloc(n*sizeof(double));
    y = (double*)mxMalloc(n*sizeof(double));
     
    /* initial x=x0 */
    for(i=0; i<n; i++){
        x[i] = x0[i];
    }
    
    
    
    /* full */
    if(1){    
       /* calculate norm(b-Ax) */
   /*  smv_diag(A.irow, A.pcol, A.val, x, array, n); array = A*x
    for(i=0; i<n; i++){
        array[i] = b[i] - array[i];
    }
    tol = tolerance(array, n);  
     */
    tstart = clock();
    it = 0;
    flag = 0;

    while(it < iter){
        for(i=0; i<n; i++){
            gamma = 0;
            for(j=A.pcol[i]; j<A.pcol[i+1]; j++){
                if(i != A.irow[j])
                    gamma += A.val[j]*x[A.irow[j]];
            }
            x[i] = (b[i]-gamma)/(A.diag[i]);
        }
        
        
        it = it+1;
        
    }   
    }
    
    
    
    
    
    
    
    
    if(0){
    /* calculate norm(b-Ax) */
    smv(A.irow, A.pcol, A.val, x, beta, n); /* array = A*x */
    for(i=0; i<n; i++){
       array[i] = b[i] - beta[i];
    }
    tol0 = tolerance(array, n); ;  /* norm(b-Ax) */ 
    tol = tol0;
    tstart = clock();
    it = 0;
    flag = 0;

    while(it < iter){
        for(i=0; i<n; i++){
            array[i] = beta[i] - A.diag[i]*x[i]; /* array = (A-D)*x */
            y[i] = (b[i]-array[i])/(A.diag[i]);
        }
        
        for(i=0; i<n; i++){
            x[i] = y[i];
        }
        
        it = it+1;
        
        /* update norm(b-Ax) */
        smv(A.irow, A.pcol, A.val, x, beta, n); /* beta = A*x */
        for(i=0; i<n; i++){
            array[i] = b[i] - beta[i]; /* array = b-A*x */
        }
       /* mexPrintf("tol:%f\n", tol);*/
        tol = tolerance(array, n); 
    }
    }
    
    
    
    
    tend = clock();
    telapsed = (double)(tend - tstart)/CLOCKS_PER_SEC;
    
/*
    mexPrintf("\nSeidel iterative method converges at %d-th iteration\n", iter);
    mexPrintf("%d Error accuracy is\n", iter);*/

 /*   mexPrintf("\nSeidel Iteration Time: %f\n", telapsed);  */
    mxFree(array);
    mxFree(beta);
    mxFree(y);
}
