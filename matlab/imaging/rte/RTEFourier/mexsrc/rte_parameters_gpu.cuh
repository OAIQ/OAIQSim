// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_cpu.h
// purpose: Contains the class definitions for the CPU version of the ballistic 
//          RTE solver 
// Author:  Nick Henscheid
// Date:    4-2017, 3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#ifndef __OAIQSIM_IMAGING_RTE_PARAMETERS_GPU_H__
#define __OAIQSIM_IMAGING_RTE_PARAMETERS_GPU_H__

#include "compute/util_gpu.cuh"
#include "compute/linalg/linear_algebra_gpu.cuh"

using namespace GPU;

// D-dimensional isotropic Gaussian attenuation function
template<int D>
struct MuFun: public Managed // For unified memory.  
{
    GPU::numtype var = 0.1;
    GPU::numtype* m;
    GPU::numtype A   = 0.0;
    MuFun(){
        CHECK(cudaMallocManaged(&m,D*sizeof(GPU::numtype)));
        for(int i=0;i<D;++i){
            m[i] = 0.5;
        }
    }
    __dh__ GPU::numtype operator()(GPU::numtype* r,GPU::numtype* s,GPU::numtype e){
        numtype exponent = 0.0;
        for(int i=0;i<D;++i){
            exponent +=pow(r[i]-m[i],(GPU::numtype)2.0);
        }
        return A*exp(-exponent/(2*var));
    }
};

// d-dimensional complex Fourier mode with wavevector k
template<int D>
struct XiFun : public Managed  // For unified memory
{
    numtype* k;  // Fourier Mode
    XiFun(numtype* k_i){
        CHECK(cudaMallocManaged(&k,D*sizeof(GPU::numtype)));
        for(int i=0;i<D;++i){k[i] = k_i[i];}
    }
    XiFun(){
        CHECK(cudaMallocManaged(&k,D*sizeof(GPU::numtype)));
        for(int i=0;i<D;++i){k[i] = 1.0;}
    }
    __dh__ GPU::complexnumtype operator()(GPU::numtype* r,GPU::numtype* s,GPU::numtype e){
        numtype exponent = 0.0;
        for(int i=0;i<D;++i){
            exponent += k[i]*r[i];
        }
        return exp((numtype)(2.0*PI)*complexnumtype(0.0,1.0)*exponent)/((numtype)(4.0*PI));
    }
};


#endif