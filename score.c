#include <mex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define H_in prhs[0]
#define r_off_in prhs[1]
#define c_off_in prhs[2]
#define ht_in prhs[3]
#define score_out plhs[0]

void check_argument(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int nL, nA, i, len, k; 
	int roff, coff;
	mwSize j, n;
	int *irow, *pcol, *r_off, *c_off;
	double *val, *Aval, *s, *ht;
	double lsep_ht, w_epq_ht;

	check_argument(nlhs, plhs, nrhs, prhs);

	/* extract field element of cell H */	
	nL = (int)mxGetScalar(mxGetField(mxGetField(H_in, 0, "Ls"), 0, "n"));
	nA = (int)mxGetScalar(mxGetField(mxGetField(H_in, 0, "A"), 0, "n"));
	irow = (int*) mxGetData(mxGetField(mxGetField(H_in, 0, "Ls"), 0, "irow"));
	pcol = (int*) mxGetData(mxGetField(mxGetField(H_in, 0, "Ls"), 0, "pcol"));
	val = (double*) mxGetPr(mxGetField(mxGetField(H_in, 0, "Ls"), 0, "val"));

	Aval = (double*) mxGetPr(mxGetField(mxGetField(H_in, 0, "A"), 0, "val"));

	r_off = (int*)mxGetPr(r_off_in);
	c_off = (int*)mxGetPr(c_off_in);
	ht = (double*)mxGetPr(ht_in);

	n = (mwSize) nL;
	len = (int)mxGetM(r_off_in);
	score_out = mxCreateNumericMatrix(len, 1, mxDOUBLE_CLASS, mxREAL);
	s = (double*) mxGetData(score_out);	

	for(i = 0; i < len; i++){
		roff = r_off[i];
		coff = c_off[i];

		lsep_ht = 0;
		for(k=pcol[roff]; k<pcol[roff+1]; k++)
			lsep_ht += val[k]*ht[irow[k]];

		w_epq_ht = Aval[i]*(ht[roff]-ht[coff]);
      		s[i] = 2*w_epq_ht*lsep_ht; 
//            
//            s[i] = 2*w_epq_ht*lsep_ht+w_epq_ht*w_epq_ht;
	}
}


void check_argument(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if(nrhs != 4 || nlhs != 1){
		mexErrMsgTxt("score.c: Wrong input or output argument number");
	}
}
