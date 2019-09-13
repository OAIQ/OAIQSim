// ------------------------------------------------------------------------
// ************************ OAIQSIM TOOLBOX *******************************
// file: rte_gpu.cuh
// purpose: Contains the class definitions for the GPU version of the 
//          ballistic RTE solver 
// Author:  Nick Henscheid
// Date:    4-2017, 10-2018, 3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:   - 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ------------------------------------------------------------------------
#ifndef __CGRI_IMAGING_RTE_GPU_H__
#define __CGRI_IMAGING_RTE_GPU_H__
#include "compute/ode/ode_solvers_gpu.cuh"   //Update to correct path once moved
#include "compute/linalg/linear_algebra_gpu.cuh"   
#include "compute/util_gpu.cuh"

namespace GPU{
template<int dim,class Domain,class MuFun,class XiFun> // Polymorphism sans abstract base class
class RTEBallistic
{
    public:
        // This is the ode right hand side that will get passed to the ode integrator
        struct OdeRHS : public Managed //Must inherit from managed for the host-side initialization!
        {
            numtype taup,taum;  // Domain entry and exit times
            numtype *r0,*s0;    // Initial point and direction
            numtype e0;         // Energy
            numtype cm = CVAC;  // Speed of light (default vacuum)
            MuFun* mu;          // Pointer to the attenuation function (small mem footprint!)
            XiFun* xi;          // Pointer to the source function      (small mem footprint!)
            OdeRHS(MuFun* mu_i,XiFun* xi_i,numtype* r0_i,numtype* s0_i,numtype e0_i,
                          numtype cm_i,numtype taup_i,numtype taum_i):
            mu(mu_i),xi(xi_i),r0(r0_i),s0(s0_i),e0(e0_i),cm(cm_i),taup(taup_i),taum(taum_i){}
            OdeRHS(){};  
            __dh__ complexnumtype operator()(numtype t,complexnumtype y){    
                numtype tempr[dim];
                for(int i=0;i<dim;++i){
                    tempr[i] = r0[i] + (t-taum)*s0[i];
                }
                return -(*mu)(tempr,s0,e0)*y + ((*xi)(tempr,s0,e0))/cm;
            }
        };
        // Constructor 
        RTEBallistic(MuFun* mu_i, XiFun* xi_i, Domain* D_i):mu(mu_i),xi(xi_i),D(D_i){}
        __dh__ void setStep(numtype h){hmin = h;}
        __dh__ void setLightSpeed(numtype c){cm = c;}

        // Single point evaluation function.  Allows for w(r,s,e) evaluation for single inputs.
        __dh__ complexnumtype operator()(numtype* r,numtype* s,numtype e){
            numtype* exittimes = (numtype*)malloc(2*sizeof(numtype));
            D->GetExitTimes(exittimes,r,s); // Get exit times
            OdeRHS* F = new OdeRHS(mu,xi,r,s,e,cm,exittimes[0],exittimes[1]);
            int nstep = max((int)ceil((exittimes[1]-exittimes[0])/hmin),nstepmin);
            if(nstep>0){
                numtype h = (exittimes[1]-exittimes[0])/((numtype)nstep-1.0); // Step size
                complexnumtype result = rk4kernel(exittimes[0],0.0,nstep,h,F);
                return result;  // Solve the ode
            } else{
                return complexnumtype(0.0);
            }
        }
        template<int d> // Dimension of input 
        complexnumvec operator()(std::vector<Eigen::Matrix<GPU::numtype,d,1> > r,
                                 std::vector<Eigen::Matrix<GPU::numtype,d,1> > s, numvec e){
            // This is a host-side function that evaluates the ray transform for a vector of
            // inputs (r,s,e).  It sets up the vector of OdeRHS, sets up the ODE integrator, 
            // then launches it.
            int neval = r.size(); // Assumes r.size() = s.size() = e.size()!
            // Set up thread grid.  
            // Note: should write a better grid selection function that 
            // selects block size based on neval.
            dim3 block(128,1);
            dim3 grid((neval+block.x-1)/block.x,1);
    
            // Get exit times
            numtype* exittimes = (numtype*)malloc(2*neval*sizeof(numtype));
            D->GetExitTimes(exittimes,r,s); 
            
            // Transfer r and s to the GPU
            numtype* r_d = vec2GPU(r,dim);
            numtype* s_d = vec2GPU(s,dim);
            
            // Set up right-hand-side vector
            OdeRHS** F_vec;
            CHECK(cudaMallocManaged(&F_vec,neval*sizeof(OdeRHS*)));
            for(int i=0;i<neval;++i){
                F_vec[i] = new OdeRHS(mu,xi,r_d+dim*i,s_d+dim*i,0.0,cm,exittimes[2*i],exittimes[2*i+1]);
            }
            
            // Set up ODE initial conditions
            // Note: since all the RTE ODEs start at (t,y) = (0,0), the ODE
            // integrator could be set up to accept a "common" initial cond
            // for a set of ODE, which would save on memory.
            complexnumtype* yi;
            CHECK(cudaMalloc(&yi,neval*sizeof(complexnumtype)));
            numtype* tinit;
            CHECK(cudaMallocManaged(&tinit,neval*sizeof(numtype)));
            initArray<<<grid,block>>>(neval,tinit,(numtype)0.0); // defined in util_gpu.cuh
            initArray<<<grid,block>>>(neval,yi,(complexnumtype)0.0);
            CHECK(cudaDeviceSynchronize());

            // Set up result arrays 
            complexnumtype *result, *result_d; 
            result = (complexnumtype*)malloc(neval*sizeof(complexnumtype));
            CHECK(cudaMalloc(&result_d,neval*sizeof(complexnumtype)));
            
            // Construct step vectors
            // Note: using mallocManaged for convenience.  Could probably be
            // made faster by doing this on the device.
            int* nstep; 
            numtype* h;  
            CHECK(cudaMallocManaged(&nstep,neval*sizeof(int)));
            CHECK(cudaMallocManaged(&h,neval*sizeof(numtype)));

            for(int i=0;i<neval;++i){
                nstep[i] = max((int)ceil(exittimes[2*i+1]/hmin),nstepmin);
                if(nstep[i]>1){
                    h[i] = exittimes[2*i+1]/((numtype)nstep[i]-1);
                }else{ 
                    h[i] = 0.0; // Shouldn't happen if nstepmin > 1.
                }
            }
            // Call the ODE integrator
            rk4<<<grid,block>>>(tinit,yi,result_d,neval,nstep,h,F_vec);
            CHECK(cudaMemcpy(result,result_d,neval*sizeof(complexnumtype),cudaMemcpyDeviceToHost));
            CHECK(cudaDeviceSynchronize());
         
            complexnumvec result_vec;
            for(int i=0;i<neval;++i){
                result_vec.push_back(result[i]);
            }
            return result_vec;
        }
        
    private: 
        MuFun* mu;     
        XiFun* xi;    
        Domain* D; 
        numtype cm = CVAC; // Vac. light speed, unless otherwise specified
        // Note: these are for a one-shot integrator; adaptive integration 
        // (e.g. RKCK) should be used w/ a convergence tolerance!!
        numtype hmin  = 0.0001;  // Unless otherwise specified
        int nstepmin = 10;       // Unless otherwise specified   
};

}// namespace GPU
#endif