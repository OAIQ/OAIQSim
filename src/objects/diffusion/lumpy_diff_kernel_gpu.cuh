// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: lumpy_kernel_gpu.cu
// purpose: LumpyKernel computes generalized lumpy background functions 
//          u(r) = B0 + sum_{j=1}^K phi(r-r_j;theta_j)
//          The standard example is the "Gaussian" lumpy background where 
//          phi(r) = (b0/pi*rb^2)*exp(-|x|^2/rb^2)
//          The "lump function" phi is defined as a device function 
// Inputs:  dim is the ambient dimension (any dim>0 possible, in theory)
//          K is the number of lumps
//          nEval is the number of evaluation points 
//          nParam is the number of parameters (number of columns in theta)
//          centers is an nEval-by-dim (row-major) array of centers (r_j)
//          evalpts is an nEval-by-dim (row-major) array of eval points (r)
//          result is an nEval-by-1 array to store the result 
//          theta is an K-by-nParam (row-major) array of parameters (theta_j)
//          The template parameter T is typically either float or double

// Author:  Nick Henscheid
// Date:    9-2016
// Contact: nhenscheid@math.arizona.edu
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#ifndef __CGRI_OBJECTS_LUMPY_DIFF_KERNEL_GPU_H__
#define __CGRI_OBJECTS_LUMPY_DIFF_KERNEL_GPU_H__
#include <math.h>
//#include "util.h"
//#define POW2 2.0

template<typename T>
__device__ T IsoDiffusionLump(int dim, T* R, T* R0, T tau, T D)
{
    // Evaluates a single (Isotropic) Diffusion lump
    // G(r,t;r0,t0) = exp(-||r-r0||^2/4D(t-t0))/(4*pi*D*(t-t0))^(d/2)
    // R is the evaluation point, R0 is the center, tau = t-t0

    T exponent = 0;
    T pow2 = 2.0;
    T powd = (T)dim/2.0;
    
    if((D>0)&&(tau>0)){
        T A = 1/powf(4.0*M_PI*D*tau,powd);  
        for(int idim = 0;idim<dim;++idim){
            exponent += powf(R[idim]-R0[idim],pow2);
        }
        return A*exp(-exponent/(4.0*D*tau)); 
    } else{
        return 0;
    }
}

template<typename T>
__global__ void LumpyDiffKernel(int dim,int K,int nEval,T* centers,T* evalpts,T* timepts,T evaltime,T* result,T C0,T D)
{
    unsigned int tId = blockIdx.x*blockDim.x + threadIdx.x; //1d thread structure
    unsigned int idx = tId*dim;
    
	if(tId<nEval){
		result[tId] = C0;  // Constant DC shift
		for(int iLump=0;iLump<K;++iLump)
        {
            result[tId] += IsoDiffusionLump(dim,&evalpts[idx],&centers[dim*iLump],evaltime - timepts[iLump],D);
        }
	}
}


#endif