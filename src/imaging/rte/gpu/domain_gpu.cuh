// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: domain_gpu.cuh
// purpose: Contains the class definition for the domain component of  
//          the RTE simulation module
// Author:  Nick Henscheid
// Date:    9-2016, 10-2018, 3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:   - 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#ifndef __CGRI_IMAGING_DOMAIN_GPU_H__
#define __CGRI_IMAGING_DOMAIN_GPU_H__

#include <compute/util_gpu.cuh>
#include <compute/linalg/linear_algebra_gpu.cuh>

namespace GPU{
// forward declarations
template<int dim>
__global__ void ComputeSphereExitTimes(numtype*,numtype*,numtype*,numtype*,int,numtype); 
template<int dim>
__dh__ numtype* ComputeRectExitTimes(numtype* L,numtype* r,numtype* s);
template<int dim>
__global__ void ComputeRectExitTimesKernel(numtype*,numtype*,numtype*,numtype*,int);

//d-dimensional disc/sphere domain. 
template<int dim> 
class SphereDomain //: public Managed // Don't think I need to inherit from managed - 
                                     // domain objects will be purely host-side?
{
    public:
        SphereDomain(numtype R_in){R = R_in;
            for(int i=0;i<dim;++i){r0[i] = 0.0;}
        }
        SphereDomain(void){R = 1.0;
            for(int i=0;i<dim;++i){r0[i] = 0.0;}
        }
        numtype radius(void){return R;}
        void set_radius(numtype R_i){R = R_i;}
        Eigen::Matrix<numtype,dim,1> center(void) {return r0;}
        void set_center(Eigen::Matrix<numtype,dim,1> r0_i){r0 = r0_i;} 
        
        void GetExitTimes(numtype* tau,numtype* r,numtype* s){
            //  WRITTEN FOR 2D - UPDATE!
            
            // Returns the forward and backward exit times for a single evaluation point 
            // tau(0) = forward, tau(1) = backward
            
//             numtype a = r[0]*s[0]+r[1]*s[1]; //dot prod;
//             numtype b = r[0]*r[0]+r[1]*r[1];
//             numtype d = pow(a,(numtype)2.0) + pow(R,(numtype)2.0) - b; // Discriminant
//             if(a>=0&&d>=0.0){
//                 tau[0] = a - sqrt(d);
//                 tau[1] = a + sqrt(d);
//             } else{
//                 tau[0] = 0.0; tau[1] = 0.0;
//             }
        }
        // Function to get exit times (for host call only!)
        void GetExitTimes(numtype* tau,const std::vector<Eigen::Matrix<numtype,dim,1> > &r,
                                       const std::vector<Eigen::Matrix<numtype,dim,1> > &s){
            // Returns the forward and backward domain exit times. 
            // tau(0) = forward, tau(1) = backward
            // If (r,s) is an outflow boundary point, tau(0) = 0.0.
            // 1) Transfer to gpu
            int n = r.size(); 
            numtype* r_d = vec2GPU(r,dim);   // Defined in linear_algebra_gpu.cuh
            numtype* s_d = vec2GPU(s,dim);
            numtype* tau_d,*r0_d; 
            CHECK(cudaMalloc(&tau_d,2*n*sizeof(numtype)));
            CHECK(cudaMalloc(&r0_d,dim*sizeof(numtype)));
            CHECK(cudaMemcpy(r0_d,r0.data(),dim*sizeof(numtype),cudaMemcpyHostToDevice));
            // 2) compute exit times on gpu
            dim3 block(128,1);// NEED n-DEPENDENT GRID SELECTION
            dim3 grid((n+block.x-1)/block.x,1);
            ComputeSphereExitTimes<dim><<<grid,block>>>(tau_d,r_d,r0_d,s_d,n,R);
            // 3) return to cpu
            CHECK(cudaMemcpy(tau,tau_d,2*n*sizeof(numtype),cudaMemcpyDeviceToHost));
        }
        
        bool isInside(Eigen::Matrix<numtype,dim,1> r){
            return ((r-r0).dot(r-r0)<R*R)? true:false;
        }

    private:
        numtype R;  //Radius of the sphere
        Eigen::Matrix<numtype,dim,1> r0;  // Center of sphere      
};  // SphereDomain

template<int dim>
__global__ void ComputeSphereExitTimes(numtype* tau,numtype* r,numtype* r0,
                                       numtype* s,int n,numtype R){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx<n){
        // Making local copies increases speed
        numtype* r_loc  = r + dim*idx;
        numtype* s_loc  = s + dim*idx;
        numtype* tau_loc = tau + 2*idx;
        numtype* r0_loc = r0;  
        numtype a=0.0,b=0.0,c=0.0,d=0.0; 
        for(int i=0;i<dim;++i){
            a+= s_loc[i]*s_loc[i];
            b+= 2.0*s_loc[i]*(r_loc[i]-r0_loc[i]); 
            c+= (r_loc[i]-r0_loc[i])*(r_loc[i]-r0_loc[i]);
        }
        c -= R*R; 
        d = b*b - 4.0*a*c;
        //printf("(a,b,c,d) = (%f,%f,%f,%f)\n",a,b,c,d);
        // Note: +b instead of -b below - tau_pm are computed using the reverse direction!
        if(a>0&&d>=0.0){
            tau[2*idx]   = (b - sqrt(d))/(2*a); 
            tau[2*idx+1] = (b + sqrt(d))/(2*a);
        }else{
            tau[2*idx] = 0.0; tau[2*idx+1] = 0.0;
        }
    }
}

/**********************************************/
/****** d-DIMENSIONAL RECTANGULAR DOMAIN ******/
/**********************************************/
template<int d>
class RectDomain
{
    public:
        RectDomain(Eigen::Matrix<numtype,d,2> L_in){
            L = L_in;
        }
        RectDomain(void){
            for(int i=0;i<d;++i){
                L(i,0) = 0.0;
                L(i,1) = 1.0; 
            }
        }
        Eigen::Matrix<numtype,d,2> dims(void){return L;}
        Eigen::Matrix<numtype,d,1> getCenter(void){
            Eigen::Matrix<numtype,d,1> result;
            for(int i=0;i<d;++i){
                result[i] = (L(i,0)+L(i,1))/2.0;
            }
            return result;
        }
        void updateDims(Eigen::Matrix<numtype,d,2> L_i){
            L = L_i;
        }
        vec2    GetExitTimes(Eigen::Matrix<numtype,d,1> &r,
                             Eigen::Matrix<numtype,d,1> &s); //Defined below
        vec2vec GetExitTimes(std::vector<Eigen::Matrix<numtype,d,1> > &r,
                             std::vector<Eigen::Matrix<numtype,d,1> > &s); //Defined below
        bool isInside(Eigen::Matrix<numtype,d,1> r){
            bool x = 1;
            for(int i=0;i<d;++i){
                x*= (r(i)>=L(i,0))&&(r(i)<=L(i,1));
            }
            return x;
        }
    private:
        Eigen::Matrix<numtype,d,2> L;       // Rectangle dimensions
};


// Exit time implementation (3D, single point evaluation)
template<int d>
vec2 RectDomain<d>::GetExitTimes(Eigen::Matrix<numtype,d,1> &r,
                                 Eigen::Matrix<numtype,d,1> &s) {
    vec2 result;
    numtype* L_arr = L.data();
    numtype* r_arr = r.data();
    numtype* s_arr = s.data();
    numtype* result_arr = ComputeRectExitTimes<d>(L_arr,r_arr,s_arr);
    return vec2(result_arr[0],result_arr[1]);
}

// Exit time implementation (3D, vector input)
template<int d>
vec2vec RectDomain<d>::GetExitTimes(std::vector<Eigen::Matrix<numtype,d,1> > &r,
                                    std::vector<Eigen::Matrix<numtype,d,1> > &s) {
    
    // 1) Transfer to gpu
    // 2) compute exit times on gpu
    // 3) return to cpu
    
    int n = r.size(); 
    vec2vec tau;
    
    printf("n = %i\n",n);
    numtype* r_d = vec2GPU(r,d);   // Defined in linear_algebra_gpu.cuh
    numtype* s_d = vec2GPU(s,d);
    numtype  *tau_h,*tau_d,*L_d; 
    tau_h = (numtype*)malloc(2*n*sizeof(numtype));
    CHECK(cudaMalloc(&tau_d,2*n*sizeof(numtype)));
    CHECK(cudaMalloc(&L_d,2*d*sizeof(numtype)));
    CHECK(cudaMemcpy(L_d,L.data(),2*d*sizeof(numtype),cudaMemcpyHostToDevice));
    // 2) compute exit times on gpu
    dim3 block(128,1);// NEED n-DEPENDENT GRID SELECTION
    dim3 grid((n+block.x-1)/block.x,1);
    ComputeRectExitTimesKernel<d><<<grid,block>>>(tau_d,L_d,r_d,s_d,n);
    // 3) return to cpu
    CHECK(cudaMemcpy(tau_h,tau_d,2*n*sizeof(numtype),cudaMemcpyDeviceToHost));
    
    for(int i=0;i<r.size();++i){
        tau.push_back(vec2(tau_h[2*i],tau_h[2*i+1]));
    }
    return tau; 
}

template<int dim>
__global__ void ComputeRectExitTimesKernel(numtype* tau,numtype* L,
                                           numtype* r,numtype* s,int n_ray){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx<n_ray){
        numtype* r_loc   = r + dim*idx;
        numtype* s_loc   = s + dim*idx;
        numtype* tau_loc = tau + 2*idx;
        tau_loc      = ComputeRectExitTimes<dim>(L,r_loc,s_loc);
        tau[2*idx]   = tau_loc[0];
        tau[2*idx+1] = tau_loc[1];
    }
}


// This is the rect domain exit time implementation.
// L is the (column-major) vector defining the box
// (r,s) is the ray
template<int dim>
__dh__ numtype* ComputeRectExitTimes(numtype* L,numtype* r,numtype* s) {
    numtype tau[2];  
    numtype inv_dir[dim];
    numtype tmin[dim];
    numtype tmax[dim];
    int sign[dim];
    for(int i=0;i<dim;++i){
        inv_dir[i] = -1./s[i]; // Negative, since computing "reverse" ray intersection
        sign[i] = (inv_dir[i] < 0);
    }
    
    for(int i=0;i<dim;++i){
        tmin[i] = (L[i+dim*sign[i]] - r[i]) * inv_dir[i];
        tmax[i] = (L[i+dim*(1-sign[i])] - r[i]) * inv_dir[i];
        if(isnan(tmin[i]))
            tmin[i] = -1.0/0.0;  // -inf
        if(isnan(tmax[i]))
            tmax[i] = 1.0/0.0;   // +inf
        if ( (tmin[0] > tmax[i]) || (tmin[i] > tmax[0]) ) 
            i=dim; // Done, exit
        if (tmin[0] > tmin[0])
            tmin[0] = tmin[i];
        if (tmax[i] < tmax[0])
            tmax[0] = tmax[i];
    }
    
    tau[0] = tmin[0];
    tau[1] = tmax[0];
    return tau;
}


} // namespace GPU

#endif