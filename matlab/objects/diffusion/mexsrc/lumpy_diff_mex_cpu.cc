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
#include "objects/diffusion/lumpy_diff_kernel_cpu.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    int dim        = (int)mxGetScalar(prhs[0]);    //Dimension (e.g. 2,3)
	int K          = (int)mxGetScalar(prhs[1]);    //Number of lumps
	int nEval      = (int)mxGetScalar(prhs[2]);    //Number of eval points'
    int auc        = (int)mxGetScalar(prhs[3]);    //Whether to evaluate "AUC" lump or not
	float *centers = (float*)mxGetData(prhs[4]);   //Lump centers
	float *evalpts = (float*)mxGetData(prhs[5]);   //Field evaluation points
    float *timepts = (float*)mxGetData(prhs[6]);   //Vector of times
    float evaltime = (float)mxGetScalar(prhs[7]);  //Evaluation time
	float C0       = (float)mxGetScalar(prhs[8]);  //Constant/"DC offset"
    float D        = (float)mxGetScalar(prhs[9]); //Diffusion coefficient
	
	plhs[0]        = mxCreateNumericMatrix(nEval,1,mxSINGLE_CLASS,mxREAL);
	float *result  = (float*)mxGetData(plhs[0]);

    //Launch kernel
	LumpyDiffKernel<float>(dim,K,nEval,centers,evalpts,timepts,evaltime,result,C0,D);
}




// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
//     int dim        = (int)mxGetScalar(prhs[0]);    //Dimension (e.g. 2,3)
// 	int K          = (int)mxGetScalar(prhs[1]);    //Number of lumps
// 	int nEval      = (int)mxGetScalar(prhs[2]);    //Number of eval points'
//     int auc        = (int)mxGetScalar(prhs[3]);    //Whether to evaluate "AUC" lump or not
// 	float *centers = (float*)mxGetData(prhs[4]);   //Lump centers
// 	float *evalpts = (float*)mxGetData(prhs[5]);   //Field evaluation points
//     float *timepts = (float*)mxGetData(prhs[6]);   //Vector of times
//     float evaltime = (float)mxGetScalar(prhs[7]);  //Evaluation time
// 	float C0       = (float)mxGetScalar(prhs[8]);  //Constant/"DC offset"
//     float D        = (float)mxGetScalar(prhs[9]); //Diffusion coefficient
// 	
// 	plhs[0]        = mxCreateNumericMatrix(nEval,1,mxSINGLE_CLASS,mxREAL);
// 	float *result  = (float*)mxGetData(plhs[0]);
//     float *centers_d,*evalpts_d,*timepts_d,*result_d;
//     
//     // For debug
//     mexPrintf("(dim,K,nEval,auc) = (%i,%i,%i,%i)\n",dim,K,nEval,auc);
//     //for(int i=0;i<12;++i){
//     //    mexPrintf("(centers[%i],evalpts[%i],theta[%i]) = (%f,%f,%f)\n",i,i,i,centers[i],evalpts[i],theta[i]);
//     //}
// 
// 	MEXCHECK(cudaMalloc((void**)&(centers_d),dim*K*sizeof(float)));
//     MEXCHECK(cudaMalloc((void**)&(evalpts_d),dim*nEval*sizeof(float)));
// 	MEXCHECK(cudaMalloc((void**)&(timepts_d),K*sizeof(float)));
//     MEXCHECK(cudaMalloc((void**)&(result_d),nEval*sizeof(float)));
//     
// 	MEXCHECK(cudaMemcpy(centers_d,centers,dim*K*sizeof(float),cudaMemcpyHostToDevice));
// 	MEXCHECK(cudaMemcpy(evalpts_d,evalpts,dim*nEval*sizeof(float),cudaMemcpyHostToDevice));
//     MEXCHECK(cudaMemcpy(timepts_d,timepts,K*sizeof(float),cudaMemcpyHostToDevice));
// 	MEXCHECK(cudaMemcpy(result_d,result,nEval*sizeof(float),cudaMemcpyHostToDevice));
// 
// 	// Cuda parameters (1d grid of 1d blocks)
// 	dim3 block(512,1);
// 	dim3 grid((nEval+block.x-1)/block.x,1);
//     mexPrintf("grid.x,grid.y = (%i,%i)\n",grid.x,grid.y);
// 	//Lauch kernel
//     if(auc==0){
//         LumpyDiffKernel<float><<<grid,block>>>(dim,K,nEval,centers_d,evalpts_d,timepts_d,evaltime,result_d,C0,D);
//     } else if(auc==1){
//         LumpyAUCKernel<float><<<grid,block>>>(dim,K,nEval,centers_d,evalpts_d,timepts_d,evaltime,result_d,C0,D);
//     } else{
//         mexErrMsgTxt("Invalid AUC parameter value");
//     }
//     MEXCHECK(cudaPeekAtLastError());
// 	MEXCHECK(cudaDeviceSynchronize());
// 	MEXCHECK(cudaMemcpy(result,result_d,nEval*sizeof(float),cudaMemcpyDeviceToHost));
// 	MEXCHECK(cudaDeviceReset());
// }