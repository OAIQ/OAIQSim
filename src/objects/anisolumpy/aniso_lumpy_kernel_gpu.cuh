// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: lumpy_kernel_cpu.h
// purpose: This header contains the CPU compute kernel for the 
//          standard lumpy background 
// Author:  Nick Henscheid
// Date:    9-2016
// Contact: nhenscheid@math.arizona.edu
// References: "Factors Influencing Lesion Detection in Medical Imaging" J.P. Rolland, 1990
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#ifndef __CGRI_OBJECTS_ANISO_LUMPY_KERNEL_GPU_H__
#define __CGRI_OBJECTS_ANISO_LUMPY_KERNEL_GPU_H__
#include <math.h>

//  LumpyKernel computes the standard d-dimensional lumpy background u(r) = B0 + sum_{j=1}^K G(r-r_j;b0,rb^2)
//	with Gaussian lump function given by G(x) = (b0/pi*rb^2)*exp(-|x|^2/rb^2)
//  the parameters K, B0, b0 and rb^2 are respectively the number of lumps,
//  DC offset, lump scale, and width (rb^2 = 2*var).
//  dim is the dimension of the set, nEval is the number of evaluation points   
//  T is typically either float or double
//  centers is K-by-dim, stored row-major; evalpts is nEval-by-dim, stored row-major
//  result will be nEval-by-1.

template<typename T>
__device__ T GaussianLump(int d, T* R, T w, T* mu, T* Linv)
{
    // Evaluates a single weighted Gaussian Lump 
    // y = w*exp(-0.5*(R-mu)'*inv(S)*(x-mu))/(det(2*pi*S))^(1/2).
    // where S = LL', i.e. L is a cholesky factor of S
    // Inputs: 
    //  d is the evaluation dimension
    //  R is the evaluation point (d-by-1)
    //  w is the weight
    //  mu is the mean vector (d-by-1)
    //  Linv is the inverse of the cholesky factor, stored in lower 
    //  triangular row-major i.e. [l11,l21,l22,l31,l32,l33,...]
    if(w==0){
        return 0.0;
    }else{ 
        T powd2 = (T)d/2.0; 
        T pow2 = 2.0;
        T A = w/powf(2*M_PI,powd2);
        T exponent = 0;
        T temp = 0;
        int idx = 0;
        // Compute determinant
        for(int i=0;i<d;++i){ 
            A *= Linv[idx];
            idx += (i+2);
        }
        idx = 0;
        // Compute exponent 
        for(int i=0;i<d;++i){
            for(int j=0;j<(i+1);++j){
                temp = temp+Linv[idx]*(R[j]-mu[j]);
                ++idx;
            }
            exponent += powf(temp,pow2);
            temp = 0.0;
        }
        return A*exp(-0.5*exponent);
    }    
}


template<typename T>
__global__ void LumpyKernel(int d,int nLump,int nEval,T B0,T* b0,T* mu,T* Linv,T* evalpts,T* result)
{
    unsigned int tId = blockIdx.x*blockDim.x + threadIdx.x; //1d thread structure
    unsigned int idx = tId*d;
    int Lsize = d*(d+1)/2;
    
	if(tId<nEval){
		result[tId] = B0;  // Constant DC shift
		for(int iLump=0;iLump<nLump;++iLump)
        {
            result[tId] += GaussianLump(d,&evalpts[idx],b0[iLump],&mu[iLump*d],&Linv[iLump*Lsize]);
                           
        }
	}
}


#endif