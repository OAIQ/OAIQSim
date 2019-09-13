// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: lumpy_mex_cpu.cc
// purpose: This is a mex CPU method to evaluate a lumpy background with specified lump centers 
//
// Author:  Nick Henscheid
// Date:    9-2016
// Contact: nhenscheid@math.arizona.edu
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include "mex.h"
#include "objects/lumpy/lumpy_kernel_cpu.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    int dim         = (int)mxGetScalar(prhs[0]);
	int K           = (int)mxGetScalar(prhs[1]); 
	int nEval       = (int)mxGetScalar(prhs[2]);   // vector of size dim
	double* centers = (double*)mxGetData(prhs[3]);
	double* evalpts = (double*)mxGetData(prhs[4]); 
	double B0       = (double)mxGetScalar(prhs[5]);
	double* b0       = (double*)mxGetData(prhs[6]);
	double* rb2      = (double*)mxGetData(prhs[7]);
    int nTotal = 1;
	plhs[0]         = mxCreateNumericMatrix(nEval,1,mxDOUBLE_CLASS,mxREAL);
	double* result  = (double*)mxGetData(plhs[0]);
	
    //Launch kernel
	LumpyKernel2D<double>(centers,evalpts,result,K,nEval,B0,b0,rb2);
}

