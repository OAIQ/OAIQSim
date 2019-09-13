// ------------------------------------------------------------------------
// ************************ OAIQSIM TOOLBOX *******************************
// file: rte_gpu.cuh
// purpose: Contains the class definitions for the CPU version of the 
//          ballistic RTE solver 
// Author:  Nick Henscheid
// Date:    4-2017, 10-2018
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:   - 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_IMAGING_RTE_GPU_H__
#define __OAIQSIM_IMAGING_RTE_GPU_H__

#include "compute/ode/ode_solvers_gpu.cuh"   
#include "compute/linalg/linear_algebra_gpu.cuh"   
#include "compute/util_gpu.cuh"

template<typename Domain,typename Mutype,typename Xitype> class RTEBallistic;
template<typename Domain,typename Mutype,typename Xitype> class RTEOdeFun;
//template<typename Domain,typename Mutype,typename Xitype> void RTEKernel(RTEBallistic* W, double* R,double* S,double* E,complexnumtype* w_res,int neval);

// template<typename Domain,typename Mutype,typename Xitype>
// __global__ void RTEKernel(RTEOdeFun<Domain,Mutype,Xitype>* F, GPU::complexnumtype* w_res,int neval)
// {
//     unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
//     if(idx<neval)
//     {
//         //SHOULD MAKE A LOCAL COPY BUT IT ISN'T WORKING FOR SOME REASON??
//         //RTEOdeFun<Domain,Mutype,Xitype> F_t = F[idx];
//         int nstep = ceil((F[idx].taum-F[idx].taup)/(F[idx].W->hmin));
//         //int nstep = 1000;
// //             printf("(taup,taum,nstep = (%2.5e,%2.5e,%i)\n",exittimes[0],exittimes[1],nstep);
//         if(nstep>0){
//             GPU::numtype h = (F[idx].taum-F[idx].taup)/((GPU::numtype)nstep-1.0); // Step size
// //                 printf("h = %2.5e\n",h);
//             w_res[idx] = rk4kernel(F[idx].taup,0.0,nstep,h,F[idx])/(F[idx].W->cm);  // Solve the ode
//         } else{
//             w_res[idx] = GPU::complexnumtype(0.0);
//         }
//     }
// }

template<typename Domain,typename Mutype,typename Xitype>
__dh__ GPU::complexnumtype RTERHSGPU(GPU::numtype t,GPU::complexnumtype y,GPU::numtype* r0,GPU::numtype* s0,GPU::numtype* tau,RTEBallistic<Domain,Mutype,Xitype>* W){   
    // tau = [taum,taup]
        GPU::numtype tempr[2];
        tempr[0] = r0[0] + (t-tau[1])*s0[0];
        tempr[1] = r0[1] + (t-tau[1])*s0[1];
        
        return -W->mu(tempr,s0,0.0)*y + W->xi(tempr,s0,0.0);
}

template<typename Domain,typename Mutype,typename Xitype>
__dh__ GPU::complexnumtype RTErk4kernel(GPU::numtype* tau,GPU::numtype* r0,GPU::numtype* s0,int ntime,GPU::numtype h,RTEBallistic<Domain,Mutype,Xitype>* W)
{
    GPU::complexnumtype y = 0.0;
    GPU::numtype t = 0.0;  //Should be 0
    GPU::complexnumtype k1,k2,k3,k4;
    //RTEBallistic* W_local = W;
    
    for(int i=0;i<ntime-1;++i) 
    {
        //printf("y.real() = %f\n",y.real());
        k1 = RTERHSGPU(t,y,r0,s0,tau,W);
        k2 = RTERHSGPU(t+ h/(GPU::numtype)2.0,y + h*k1/(GPU::numtype)2.0,r0,s0,tau,W);
        k3 = RTERHSGPU(t+ h/(GPU::numtype)2.0,y + h*k2/(GPU::numtype)2.0,r0,s0,tau,W);
        k4 = RTERHSGPU(t+ h,y + h*k3,r0,s0,tau,W);
        y = y + (h/(GPU::numtype)6.0)*(k1 + (GPU::numtype)2.0*k2 + (GPU::numtype)2.0*k3 + k4);   
        t = t + h;     
    }
   
    return y;
}

template<typename Domain,typename Mutype,typename Xitype>
__global__ void RTEKernel(RTEBallistic<Domain,Mutype,Xitype>* W,GPU::numtype* R,GPU::numtype* S, GPU::complexnumtype* w_res,int neval)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx<neval)
    {
        RTEBallistic<Domain,Mutype,Xitype>* W_local = W;
        GPU::numtype* r0 = R+2*idx;
        GPU::numtype* s0 = S+2*idx;  
        GPU::numtype tau[2];
        W_local->D.GetExitTimes(tau,r0,s0);
        //SHOULD MAKE A LOCAL COPY BUT IT ISN'T WORKING FOR SOME REASON??
        //RTEOdeFun<Domain,Mutype,Xitype> F_t = F[idx];
        int nstep = ceil((tau[1]-tau[0])/W_local->hmin);
        //printf("r0[0],r0[1],tau[0],tau[1],nstep = (%f,%f,%f,%f,%i)\n",r0[0],r0[1],tau[0],tau[1],nstep);
        //int nstep = 1000;
//             printf("(taup,taum,nstep = (%2.5e,%2.5e,%i)\n",exittimes[0],exittimes[1],nstep);
        if(nstep>0){
            GPU::numtype h = (tau[1]-tau[0])/((GPU::numtype)nstep-1.0); // Step size
//                 printf("h = %2.5e\n",h);
            w_res[idx] = RTErk4kernel(tau,r0,s0,nstep,h,W_local)/(W_local->cm);  // Solve the ode
            //printf("w_res[%i] = %f\n",idx,w_res[idx]);
        } else{
            w_res[idx] = GPU::complexnumtype(0.0);
        }
    }
}

// template<class MuFun,class XiFun>
// struct RTEOdeFun:public Managed   // For unified memory
//         {
//             GPU::numtype taup;   // Domain entry time
//             GPU::numtype taum;   // Domain exit time
//             GPU::numtype* r0;    // Initial point
//             GPU::numtype* s0;    // Initial direction
//             GPU::numtype e0;     // Energy
//             GPU::numtype cm;     // Speed of light
//             MuFun mu;
//             XiFun xi;
//             __dh__ RTEOdeFun(MuFun &mu_in,XiFun &xi_in,GPU::numtype cm_in){
//                 mu = mu_in;xi = xi_in;cm = cm_in;
//                 CHECK(cudaMallocManaged(&r0,2*sizeof(GPU::numtype)));
//                 CHECK(cudaMallocManaged(&s0,2*sizeof(GPU::numtype)));
//                 //memcpy(r0,r0_in,2*sizeof(numtype));   //Set elsewhere
//                 //memcpy(s0,s0_in,2*sizeof(numtype));
//             }
//             __dh__ RTEOdeFun(){};
//             __dh__ RTEOdeFun(RTEOdeFun &F_i){
//                 CHECK(cudaMallocManaged(&r0,2*sizeof(GPU::numtype)));
//                 CHECK(cudaMallocManaged(&s0,2*sizeof(GPU::numtype)));
//                 memcpy(r0,F_i.r0,2*sizeof(GPU::numtype));   //Set elsewhere
//                 memcpy(s0,F_i.s0,2*sizeof(GPU::numtype));
//             }
//             __dh__ GPU::complexnumtype operator()(GPU::numtype t,GPU::complexnumtype y){     
//                 GPU::numtype tempr[2];
//                 tempr[0] = r0[0] + (t-taup-taum)*s0[0];
//                 tempr[1] = r0[1] + (t-taup-taum)*s0[1];
//                 return -mu(tempr,s0,e0)*y + (xi(tempr,s0,e0));
//             }
//         };
        
template<class Domain,class Mutype,class Xitype>
struct RTEOdeFun : public Managed   // For unified memory
{
    GPU::numtype taup = 0.0;   // Domain entry time
    GPU::numtype taum = 0.0;   // Domain exit time
    GPU::numtype* r0;    // Initial point
    GPU::numtype* s0;    // Initial direction
    GPU::numtype e0=0.0;     // Energy
    //GPU::numtype cm;     // Speed of light
    RTEBallistic<Domain,Mutype,Xitype>* W;
    //MuFun mu;
    //XiFun xi;
    __dh__ RTEOdeFun(RTEBallistic<Domain,Mutype,Xitype>* W_in)
    {
        //CHECK(cudaMallocManaged(&W,sizeof(RTEBallistic<Domain,MuFun,XiFun>)));
        W = W_in;
        //mu = mu_in;xi = xi_in;cm = cm_in;
        CHECK(cudaMallocManaged(&r0,2*sizeof(GPU::numtype)));
        CHECK(cudaMallocManaged(&s0,2*sizeof(GPU::numtype)));
        //memcpy(r0,r0_in,2*sizeof(numtype));   //Set elsewhere
        //memcpy(s0,s0_in,2*sizeof(numtype));
    }
    __dh__ RTEOdeFun(){};
    __dh__ RTEOdeFun(RTEOdeFun &F_i)
    {
        W = F_i.W;
        taup = F_i.taup;
        taum = F_i.taum;
        e0   = F_i.e0;
        CHECK(cudaMallocManaged(&r0,2*sizeof(GPU::numtype)));
        CHECK(cudaMallocManaged(&s0,2*sizeof(GPU::numtype)));
        memcpy(r0,F_i.r0,2*sizeof(GPU::numtype));   
        memcpy(s0,F_i.s0,2*sizeof(GPU::numtype));
    }
    __dh__ GPU::complexnumtype operator()(GPU::numtype t,GPU::complexnumtype y)
    {     
        GPU::numtype tempr[2];
        tempr[0] = r0[0] + (t-taup-taum)*s0[0];
        tempr[1] = r0[1] + (t-taup-taum)*s0[1];
        return -W->mu(tempr,s0,e0)*y + (W->xi(tempr,s0,e0));
    }
};


template<class Domain,class Mutype,class Xitype> // Polymorphism sans abstract base class
class RTEBallistic : public Managed
{
    public:
        // This is the ode right hand side that will get passed to the ode integrator
        
        // Constructor 
        __dh__ RTEBallistic(Mutype &mu_in, Xitype &xi_in, Domain &D_in){
            CHECK(cudaMallocManaged((void**)&mu,sizeof(Mutype)));
            CHECK(cudaMallocManaged((void**)&xi,sizeof(Xitype)));
            CHECK(cudaMallocManaged((void**)&D,sizeof(Domain)));
            //CHECK(cudaMallocManaged((void**)&F,sizeof(RTEOdeFun<MuFun,XiFun>)));
            mu = mu_in;xi = xi_in;D = D_in;
            //F = RTEOdeFun<Domain,MuFun,XiFun>(mu,xi,cm);
        }
        __dh__ void SetStep(GPU::numtype h){hmin = h;}
        __dh__ void SetLightSpeed(GPU::numtype c){cm = c;}//F.cm=c;}
        // NEED TO HAVE A VERSION THAT ALLOWS EVALUATION AT A VECTOR OF INPUTS
        // Evaluation function.  Allows for simple w(r,s,e) evaluation. 
        __dh__ GPU::complexnumtype operator()(GPU::numtype* r,GPU::numtype* s,GPU::numtype e){
            GPU::numtype exittimes[2];
            D.GetExitTimes(exittimes,r,s); // Get exit times
            //RTEOdeFun<MuFun,XiFun> *F_loc = new RTEOdeFun<MuFun,XiFun>(F);
            RTEOdeFun<Domain,Mutype,Xitype>* F = new RTEOdeFun<Domain,Mutype,Xitype>(this);
            F->r0 = r; // Set up ode RHS
            F->s0 = s; // Set up ode RHS
            F->taup = exittimes[0]; // Entry time 
            F->taum = exittimes[1]; // Exit time 
            //printf("tau = (%f,%f)\n",exittimes[0],exittimes[1]);
            int nstep = ceil((exittimes[1]-exittimes[0])/hmin);
            //int nstep = 1000;
//             printf("(taup,taum,nstep = (%2.5e,%2.5e,%i)\n",exittimes[0],exittimes[1],nstep);
            if(nstep>0){
                GPU::numtype h = (exittimes[1]-exittimes[0])/((GPU::numtype)nstep-1.0); // Step size
//                 printf("h = %2.5e\n",h);
                return rk4kernel(exittimes[0],0.0,nstep,h,*F)/cm;  // Solve the ode
            } else{
                return GPU::complexnumtype(0.0);
            }
        }
        
        GPU::complexnumvec operator()(GPU::vec2vec r,GPU::vec2vec s,GPU::numvec e)
        {
            // Set up kernel launch params
            int neval = r.size();
            dim3 grid((neval+blksz.x-1)/blksz.x,1);
            
            GPU::numtype* r_d = GPU::vec2GPU(r,2);
            GPU::numtype* s_d = GPU::vec2GPU(s,2);
            GPU::numtype* e_h = e.data();
            GPU::numtype* e_d; 
            CHECK(cudaMalloc(&e_d,neval*sizeof(GPU::numtype)));
            CHECK(cudaMemcpy(e_d,e_h,neval*sizeof(GPU::numtype),cudaMemcpyHostToDevice));
            GPU::complexnumtype* res_d;
            CHECK(cudaMallocManaged(&res_d,neval*sizeof(GPU::complexnumtype)));
            
//             for(int i=0;i<32;++i){
//                 printf("r_d[2*i],r_d[2*i+1] = (%f,%f)\n",r_d[2*i],r_d[2*i+1]);
//             }
            // Lanuch kernel
            
            //GPU::numtype exittimes[2];
//             RTEOdeFun<Domain,Mutype,Xitype> *F;
//             CHECK(cudaMallocManaged((void**)&F,neval*sizeof(RTEOdeFun<Domain,Mutype,Xitype>)));
//             printf("Sizeof(RTEOdeFun) = %i\n", sizeof(RTEOdeFun<Domain,Mutype,Xitype>));
//             for(int i=0;i<neval;++i){
//                 F[i] = RTEOdeFun<Domain,Mutype,Xitype>(this);
//                 GPU::numtype* r = r_d + 2*i;
//                 GPU::numtype* s = s_d + 2*i;
//                 D.GetExitTimes(exittimes,r,s);
//                 F[i].r0 = r; // Set up ode RHS
//                 F[i].s0 = s; // Set up ode RHS
//                 F[i].taup = exittimes[0]; // Entry time 
//                 F[i].taum = exittimes[1]; // Exit time 
//                 //printf("tau = (%f,%f)\n",exittimes[0],exittimes[1]);
//             }
            
            RTEKernel<<<grid,blksz>>>(this,r_d,s_d,res_d,neval);
            CHECK(cudaDeviceSynchronize());
            //GPU::numtype* r_temp;
            //GPU::numtype* s_temp;
            GPU::complexnumvec res;
            for(int i=0;i<neval;++i){
//                 r_temp = r_d+2*i;
//                 s_temp = s_d+2*i;
//                 printf("computing %i/%i\n",i+1,neval);
//                 //printf("res_d.real() = %f\n",res_d[i].real());
//                 res.push_back(this->operator()(r_temp,s_temp,0.0));
                   res.push_back(res_d[i]);
            }
            return res;
        }
 
        Mutype mu;     
        Xitype xi;    
        //RTEOdeFun<MuFun,XiFun> F; // Don't think this will be needed.
        Domain D; 
        GPU::numtype cm = 299792458;   // Unless otherwise specified
        GPU::numtype hmin  = 0.0001;   // Unless otherwise specified
        dim3 blksz = dim3(128,1);        // Unless otherwise specified
        //int nstepmin = 10;      // Unless otherwise specified 
        
};






#endif