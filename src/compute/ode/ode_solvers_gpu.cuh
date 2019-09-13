// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// file: ode_solvers.h
// purpose: Contains various ODE solvers.  Used mainly by the RTE 
//          solver in the imaging module.
// Author:  Nick Henscheid
// Date:    4-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_COMPUTE_ODE_SOLVERS_GPU_H__
#define __OAIQSIM_COMPUTE_ODE_SOLVERS_GPU_H__
#include <compute/util_gpu.cuh>

#include "compute/linalg/linear_algebra_gpu.cuh"
//using namespace GPU;

//Notes: this function only returns the final y value, not the entire vector.  
//       this is obviously a special case

template<class OdeRhs> //This is what allows polymorphism without abstract base classes, which CUDA can't handle  
__dh__ GPU::complexnumtype rk4kernel(GPU::numtype tinit,GPU::complexnumtype yinit,int ntime,GPU::numtype h,OdeRhs* F)
{
    GPU::complexnumtype y = yinit;
    GPU::numtype t = tinit;
    GPU::complexnumtype k1,k2,k3,k4;
    for(int i=0;i<ntime-1;++i) 
    {
        k1 = (*F)(t,y);
        k2 = (*F)(t+ h/(GPU::numtype)2.0,y + h*k1/(GPU::numtype)2.0);
        k3 = (*F)(t+ h/(GPU::numtype)2.0,y + h*k2/(GPU::numtype)2.0);
        k4 = (*F)(t+ h,y + h*k3);
        y = y + (h/(GPU::numtype)6.0)*(k1 + (GPU::numtype)2.0*k2 + (GPU::numtype)2.0*k3 + k4);   
        t = t + h;     
        //printf("(k1,t,y) = (%f,%f,%4.16f+%4.16fi)\n",k1,t,y.real(),y.imag());
    }
    return y;
}

//  Might need deep copy, thrust, or unified memory to implement this for more complicated OdeRhs 
//  (Specifically if OdeRhs requires an array in its definition)
template<class OdeRhs>
__global__ void rk4(GPU::numtype* tinit,GPU::complexnumtype* yi,GPU::complexnumtype* yf,int n_ode,int* n_time,GPU::numtype* h,OdeRhs** F)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    
    if(idx<n_ode){
        //printf("Hello! (idx,tinit[idx],yi[idx]) = (%i,%8.8f,%f+%fi)\n",idx,tinit[idx],yi[idx].real(),yi[idx].imag());
        yf[idx] = rk4kernel(tinit[idx],yi[idx],n_time[idx],h[idx],F[idx]);
        //printf("Result[idx] = %8.16f+%8.16fi\n",yf[idx].real(),yf[idx].imag());
    }
}



// namespace GPU{
//     class RK4Integrator{
//         public:
//             RK4Integrator(void){n_step = 1000;}
//             void SetStep(int n){n_step = n;}
// 
//             template<typename rhs>
//             GPU::complexnumvec operator()(std::vector<rhs>& rhs_vec,GPU::vec2vec& Tinterval,GPU::complexnumvec& y0){
//                 int N = rhs_vec.size();
//                 GPU::complexnumvec result;
//                 GPU::numtype ti = Tinterval[0](0,0);
//                 GPU::numtype tf = Tinterval[0](1,0);
//                 GPU::numtype h = (tf-ti)/n_step;
//                 printf("(N,ti,tf,y0,h) = (%i,%f,%f,%f,%f)\n",N,ti,tf,y0[0].real(),h);
//                 // Need to copy rhs vector to device
//                 // then launch __global__ rk4, not the device integrator
//                 rhs* rhs_vec_h = (rhs*)malloc(N*sizeof(rhs));
//                 rhs* rhs_vec_d;
//                 GPU::numtype* ti_h  = (GPU::numtype*)malloc(N*sizeof(GPU::numtype));
//                 GPU::numtype* ti_d;
//                 GPU::complexnumtype *yi_h = (GPU::complexnumtype*)malloc(N*sizeof(GPU::complexnumtype));
//                 GPU::complexnumtype *yf_h = (GPU::complexnumtype*)malloc(N*sizeof(GPU::complexnumtype));
//                 GPU::complexnumtype *yi_d, *yf_d;
//                     
//                 for(int i=0;i<N;++i){
//                     rhs_vec_h[i] = rhs_vec[i];
//                     ti_h[i] = Tinterval[i](0,0);
//                     yi_h[i] = y0[i];
//                 }
//                 CHECK(cudaMallocManaged((void**)&rhs_vec_d,N*sizeof(rhs)));
//                 CHECK(cudaMalloc((void**)&ti_d,N*sizeof(GPU::numtype)));
//                 CHECK(cudaMalloc((void**)&yi_d,N*sizeof(GPU::complexnumtype)));
//                 CHECK(cudaMalloc((void**)&yf_d,N*sizeof(GPU::complexnumtype)));
//                 // This assumes that no deep copying is necessary!
//                 CHECK(cudaMemcpy(rhs_vec_d,rhs_vec.data(),N*sizeof(rhs),cudaMemcpyHostToDevice));
//                 CHECK(cudaMemcpy(ti_d,rhs_vec_h,N*sizeof(GPU::numtype),cudaMemcpyHostToDevice));
//                 CHECK(cudaMemcpy(yi_d,yi_h,N*sizeof(GPU::complexnumtype),cudaMemcpyHostToDevice));
//                
//                 rk4<<<1,N>>>(ti_d,yi_d,yf_d,N,n_step,h,rhs_vec_d);
//                 CHECK(cudaMemcpy(yf_h,yf_d,N*sizeof(GPU::complexnumtype),cudaMemcpyDeviceToHost));
//                 printf("yf_h[0] = %f\n",yf_h[0].real());
//                 CHECK(cudaDeviceSynchronize());
//                 //GPU::complexnumtype res0 = rk4kernel(Tinterval[0](0,0),y0[0],n_step,h,rhs_vec_d[0]);
//                 for(int i=0;i<N;++i){
//                     result.push_back(yf_h[i]);
//                 }
//                 printf("result[0].real() = %f\n",result[0].real());
//                 CHECK(cudaDeviceSynchronize());
//                 //rk4kernel(rhs_vec,result,Tinterval,y0,n_step);
//                 return result;
//             }
//         private:
//             int n_step;
//     };
// }


#endif  