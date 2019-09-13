// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_parameters_gpu.h
// purpose: Contains the class definitions for the CPU version of the ballistic 
//          RTE solver 
// Author:  Nick Henscheid
// Date:    4-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#ifndef __CGRI_IMAGING_RTE_PARAMETERS_GPU_H__
#define __CGRI_IMAGING_RTE_PARAMETERS_GPU_H__

#include "compute/util_gpu.cuh"
#include "compute/linalg/linear_algebra_gpu.cuh"

using namespace GPU;

struct MuFun
{
    GPU::numtype var = 0.05;
    GPU::numtype mx  = 0.0;
    GPU::numtype my  = 0.0;
    GPU::numtype A   = 0.0;
    __dh__ GPU::numtype operator()(GPU::numtype* r,GPU::numtype* s,GPU::numtype e){
        return 0.0;
        //return A*exp(-(pow(r[0]-mx,2.0)+pow(r[1]-my,2.0))/(2*var));
    }
};

// template<typename T>
// struct XiFun
// {
//     T* k;
//     __dh__ XiFun(T k1,T k2){
//         T temp[2]; temp[0]=k1;temp[1]=k2;
//         k = temp;
//     }
//     __dh__ XiFun(){};
//     __dh__ T operator()(T* r,T* s,T e){
//         return cos(2*M_PI*(k[0]*r[0]+k[1]*r[1]))/(4.0*M_PI);
//     }
// };

// 2-dimensional complex Fourier mode with wavevector k

struct XiFun
{
    GPU::numtype k[2];  // Mode
    __dh__ XiFun(GPU::numtype k1,GPU::numtype k2){
        k[0]=k1;k[1]=k2;
    }
    __dh__ XiFun(){};
    __dh__ GPU::complexnumtype operator()(GPU::numtype* r,GPU::numtype* s,GPU::numtype e){
        //GPU::complexnumtype I(0.0,1.0);
        return exp((numtype)(2*M_PI*(k[0]*r[0]+k[1]*r[1]))*(GPU::complexnumtype(0.0,1.0)))/((numtype)(4.0*M_PI));
    }
};
#endif