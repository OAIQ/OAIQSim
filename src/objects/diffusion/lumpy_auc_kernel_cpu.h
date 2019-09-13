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
#ifndef __CGRI_OBJECTS_LUMPY_KERNEL_CPU_H__
#define __CGRI_OBJECTS_LUMPY_KERNEL_CPU_H__
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
void LumpyKernel2D(T* centers,T* evalpts, T* result,int K,int nEval,T B0,T* b0,T* rb2)
{
    T exponent = 0;
    T power2 = 2.0; //Must convert to T else strange bugs appear
   
    T xEval,yEval,A;

    for(int ix=0;ix<nEval;++ix){
        xEval = evalpts[2*ix];
        yEval = evalpts[2*ix+1];
        result[ix] = B0;
        for(int iLump=0;iLump<K;++iLump){
            if((rb2[iLump]>0)&&(b0[iLump]>0)){
                A = (b0[iLump]/(M_PI*rb2[iLump])); 
                exponent = pow(xEval-centers[2*iLump],power2)+pow(yEval-centers[2*iLump+1],power2); 
                result[ix] += A*exp(-exponent/rb2[iLump]); 
            } else{
                //printf("Invalid input\n");
                result[ix] += 0.0;
            }
            
        }
    }
}


// for(int iEval=0;iEval<nEval;++iEval){
//         result[iEval] = B0;
//         for(int iLump=0;iLump<K;++iLump){
//             exponent = 0;
//             for(int iDim=0;iDim<dim;++iDim){
//                 exponent += pow(evalpts[dim*iEval+iDim]-centers[dim*iLump+iDim],power2); 
//             }
//             result[iEval] += A*exp(-exponent/rb2); 
//         }
//     }


#endif