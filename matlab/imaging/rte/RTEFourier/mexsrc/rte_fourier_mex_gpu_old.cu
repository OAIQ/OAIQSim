// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_fourier_cpu_mex.h
// purpose: Contains the mex interface to compute ballistic RTE solutions  
// Author:  Nick Henscheid
// Date:    4-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include<mex.h>
#include<iostream>
#include "imaging/rte/gpu/rte_gpu.cuh"  // RTE Class file 
#include "imaging/rte/gpu/domain_gpu.cuh"  // Domain definitions
#include "rte_parameters_gpu.cuh"  // Custom parameters
#include "compute/util_gpu.cuh"

// This is a specialized mex interface to compute attenuated ray transforms of Fourier basis elements 
// Phi_k(r) = exp(2*pi*i*k\cdot r)chi(r) 

__global__ void ComputeRTE(double* R,double* S,GPU::complexnumtype* w_res,int neval,double k1,double k2,double hmin,double cm)
{
    DiscDomain D; //Standard unit disc domain
    MuFun mu; mu.A = 0.0; //No attenuation
    XiFun xi(k1,k2); 
    RTEBallistic<DiscDomain,MuFun,XiFun> w(mu,xi,D);
    w.SetStep(hmin);
    w.SetLightSpeed(cm);
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    GPU::numtype r[2];
    GPU::numtype s[2];
    if(idx<neval){
        // This is 2D.  Should make it generic.
        r[0] = R[2*idx];
        r[1] = R[2*idx+1];
        s[0] = S[2*idx];
        s[1] = S[2*idx+1];
        w_res[idx] = w(r,s,0.0);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    double* R  = (double*)mxGetData(prhs[0]);  // positions  (2xneval, col major)
    double* S  = (double*)mxGetData(prhs[1]);  // directions (2xneval, col major)
    double* k  = (double*)mxGetData(prhs[2]);  // Wavenumber 
    int neval  = (int)mxGetScalar(prhs[3]);    // Number of evaluation points
    double cm  = (double)mxGetScalar(prhs[4]); // Speed of light in the medium
    double hmin  = (double)mxGetScalar(prhs[5]);
    
    printf("neval, cm, hmin, kx, ky = %i,%f,%f,%f,%f\n",neval,cm,hmin,k[0],k[1]);

    plhs[0]    = mxCreateDoubleMatrix(neval, 1, mxCOMPLEX);
    double* Wr  = (double*)mxGetPr(plhs[0]);
    double* Wi  = (double*)mxGetPi(plhs[0]);

    double *R_d,*S_d;
    GPU::complexnumtype *w_res_d,*w_res;
    w_res = (GPU::complexnumtype *)malloc(neval*sizeof(GPU::complexnumtype));
    CHECK(cudaMalloc((void**)&(w_res_d),neval*sizeof(GPU::complexnumtype)));
    CHECK(cudaMalloc((void**)&(R_d),2*neval*sizeof(double)));
    CHECK(cudaMalloc((void**)&(S_d),2*neval*sizeof(double)));
    CHECK(cudaMemcpy(R_d,R,2*neval*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(S_d,S,2*neval*sizeof(double),cudaMemcpyHostToDevice));
    
    for(int i=0;i<neval;++i)
    {
        mexPrintf("R[%i] = (%f,%f), ",i,R[2*i],R[2*i+1]);
        mexPrintf("S[%i] = (%f,%f)\n",i,S[2*i],S[2*i+1]);
    }
    
    dim3 block(512,1);
	dim3 grid((neval+block.x-1)/block.x,1);
    mexPrintf("(grid,block) = ((%i,%i),(%i,%i))\n",grid.x,grid.y,block.x,block.y);
    //ComputeRTE<<<grid,block>>>(R_d,S_d,w_res_d,neval,k[0],k[1],hmin,cm);
    CHECK(cudaMemcpy(w_res,w_res_d,neval*sizeof(GPU::complexnumtype),cudaMemcpyDeviceToHost));
    CHECK(cudaDeviceSynchronize());

    for(int i=0;i<neval;++i){
         Wr[i] = w_res[i].real();
         Wi[i] = w_res[i].imag();
         //printf("(w_res[i].real(),w_res[i].imag()) = (%2.3e,%2.3e)\n",Wr[i],Wi[i]);
    }

}

