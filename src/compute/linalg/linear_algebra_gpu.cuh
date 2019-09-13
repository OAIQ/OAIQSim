// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// file: linear_algebra_gpu.h
// purpose: Contains various typedefs, includes and utility functions for 
//          linear algebra operations within the toolbox.
// Author:  Nick Henscheid
// Date:    5-2017, 10-2018
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:      
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_COMPUTE_LIN_ALG_GPU_H__
#define __OAIQSIM_COMPUTE_LIN_ALG_GPU_H__

// Linear algebra includes
//#include<vector>
//#include "compute/extern/cuda_complex/cuda_complex.hpp"
//#include <thrust/complex.h>
#include <compute/util_gpu.cuh>
#include <compute/extern/eigen3/Eigen/Dense>

namespace GPU{
    typedef Eigen::Matrix<GPU::numtype,2,1> vec2;
    typedef Eigen::Matrix<GPU::numtype,3,1> vec3;
    typedef Eigen::Matrix<GPU::numtype,4,1> vec4;
    typedef Eigen::Matrix<GPU::complexnumtype,2,1> complexvec2;
    typedef Eigen::Matrix<GPU::complexnumtype,3,1> complexvec3;
    typedef Eigen::Matrix<GPU::complexnumtype,4,1> complexvec4;
    typedef std::vector<vec2> vec2vec;
    //typedef thrust::device_vector<vec2> vec2vec_d;
    typedef std::vector<vec3> vec3vec;
    typedef std::vector<vec4> vec4vec;
    typedef std::vector<complexvec2> complexvec2vec;
    typedef std::vector<complexvec3> complexvec3vec;
    typedef std::vector<complexvec4> complexvec4vec;
    typedef std::vector<GPU::numtype> numvec; 
    typedef std::vector<GPU::complexnumtype> complexnumvec; 

    template<typename V>
    GPU::numtype* vec2GPU(V Xvec_h,int dim)
    {
        // Converts a vecnvec from host to GPU & returns the pointer to the device array
        int N = Xvec_h.size();
        GPU::numtype* X;
        CHECK(cudaMallocManaged((void**)&X,dim*N*sizeof(GPU::numtype)));
        if(dim>1){
            for(int i=0;i<N;++i){
                for(int j=0;j<dim;++j){
                    X[dim*i+j] = Xvec_h[i](j);
                }
            }
        }
        //GPU::numtype* X_d;
        
        //CHECK(cudaMemcpy(X_d,X_h,dim*N*sizeof(GPU::numtype),cudaMemcpyHostToDevice));
        return X;
    }
}

//WHY IS ALL THIS COMMENTED OUT...?

// Utility functions 
// 
// //Function to convert a C-style array to an std::vector of Eigen::Vectorxd.
// template<int d> //Dimension of vector
// std::vector<Eigen::Matrix<numtype,d,1> > ArrayToVec(const numtype* x_in,int n)
// {
//     Eigen::Matrix<numtype,d,1> v_temp;
//     std::vector<Eigen::Matrix<numtype,d,1> > x_out;
//     for(int i=0;i<n;++i){
//         for(int j=0;j<d;++j){
//             v_temp(j) = x_in[d*i+j];
//         }
//         x_out.push_back(v_temp);
//     }
//     return x_out;
// }
// // Function to convert std::vector of Eigen::Vectorxd to C-style array
// template<int d> //Dimension of vector
// numtype* VecToArray(const std::vector<Eigen::Matrix<numtype,d,1> >& x_in)
// {
//     int n = x_in.size();
//     numtype* x_out = (numtype*) malloc(n*d*sizeof(numtype));
//     
//     for(int i=0;i<n;++i){
//         for(int j=0;j<d;++j){
//             x_out[d*i+j] = x_in[i](j);
//         }
//     }
//     return x_out;
// }
// // Function to convert a C-style array to an std::vector of numtypes
// std::vector<numtype> ArrayToNumVec(const numtype* x_in,int n)
// {
//     std::vector<numtype> x_out;
//     for(int i=0;i<n;++i){
//         x_out.push_back(x_in[i]);
//     }
//     return x_out;
// }
// // Function to convert std::vector of numtypes to a C-style array
// void NumVecToArray(const std::vector<numtype>& x_in,numtype* x_out)
// {
//     // Warning! Assumes that x_out is already allocated! bad design. fix.
//     int n = x_in.size();
//     for(int i=0;i<n;++i){
//         x_out[i] = x_in[i];
//     }
// }



#endif