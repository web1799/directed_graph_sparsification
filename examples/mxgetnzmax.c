/*=================================================================
 * mxgetnzmax.c
 *
 * This example illustrates how to use mxGetNzMax.  It takes a sparse
 * matrix as an input argument and it displays the number of nonzero
 * elements in the input argument and the maximum number of nonzero
 * elements that can be stored.  The maximum number of elements that
 * can be stored is the value of nzmax.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2011 The MathWorks, Inc.  
 * All rights reserved.
 *=================================================================*/
/* $Revision: 1.5.6.3 $ */
#include "mex.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mwSize nzmax; 
    mwSize columns;
    mwIndex nnz;
    
    (void) plhs;    /* unused parameter */
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:mxgetnzmax:invalidNumInputs",
                "One input argument required.");
    }
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:mxgetnzmax:maxlhs",
                "Too many output arguments.");
    }
    
    if (!mxIsSparse(prhs[0]))  {
        mexErrMsgIdAndTxt( "MATLAB:mxgetnzmax:invalidInputSparsity",
                "Input argument must be a sparse array.");
    }
    nzmax = mxGetNzmax(prhs[0]);
    columns = mxGetN(prhs[0]);
    
    /* NOTE: nnz is the actual number of nonzeros and is stored as the
     * last element of the jc array where the size of the jc array is the
     * number of columns + 1 */
    nnz = *(mxGetJc(prhs[0]) + columns);
    
    mexPrintf("Contains %d nonzero elements.\n", nnz);
    mexPrintf("Can store up to %d nonzero elements.\n", nzmax);
}
