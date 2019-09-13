// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// file: domain_cpu.h
// purpose: Contains the class definition for the disc domain component of  
//          the RTE simulation module
// Author:  Nick Henscheid
// Date:    9-2016, 11-2018, 3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_IMAGING_DOMAIN_CPU_H__
#define __OAIQSIM_IMAGING_DOMAIN_CPU_H__

#include "compute/util_cpu.h"
#include "compute/linalg/linear_algebra_cpu.h"

namespace CPU{
/**********************************************/
/****** d-DIMENSIONAL SPHERE/DISC DOMAIN ******/
/**********************************************/
template<int dim>
class SphereDomain
{
    public:
        SphereDomain(numtype R_in){R = R_in;
            for(int i=0;i<dim;++i){r0[i] = 0.0;}
        }
        SphereDomain(void){R = 1.0;
            for(int i=0;i<dim;++i){r0[i] = 0.0;}
        }
        numtype radius(void)  {return R;}
        void set_radius(numtype r){R = r;}
        Eigen::Matrix<numtype,dim,1> center(void) {return r0;}
        void set_center(Eigen::Matrix<numtype,dim,1> r0_i) {r0 = r0_i;}
        
        // Exit time functions, defined below. 
        vec2 GetExitTimes(const Eigen::Matrix<numtype,dim,1> &r,
                          const Eigen::Matrix<numtype,dim,1> &s); // Single eval version
        vec2vec GetExitTimes(const std::vector<Eigen::Matrix<numtype,dim,1> > &r,
                             const std::vector<Eigen::Matrix<numtype,dim,1> > &s); // Multiple eval
        // Domain indicator function
        bool IsInside(vec2 r){
            return (r.dot(r)<R*R)? true:false;
        }
    private:
         numtype R;  //Radius of the circle
         Eigen::Matrix<numtype,dim,1> r0;  // Center of sphere   
};
// Exit time implementation (single eval point)
template<int dim>
vec2 SphereDomain<dim>::GetExitTimes(const Eigen::Matrix<numtype,dim,1> &r,
                                     const Eigen::Matrix<numtype,dim,1> &s)
{
    // Returns the forward and backward exit times for a single 
    // phase point (r,s)
    // tau(0) = forward, tau(1) = backward
    vec2 tau;
    numtype a=0.0,b=0.0,c=0.0,d=0.0; 
    Eigen::Matrix<numtype,dim,1> rtemp = r-r0;
    a = s.dot(s);
    b = 2.0*s.dot(rtemp);
    c = rtemp.dot(rtemp) - R*R;
    d = b*b - 4.0*a*c;
    if(a>0&&d>=0.0){
            tau[0] = (b - sqrt(d))/(2*a); 
            tau[1] = (b + sqrt(d))/(2*a);
    }else{
            tau[0] = 0.0; tau[1] = 0.0;
    }
    return tau;
}
// Exit time implementation (vector of eval points) 
template<int dim>
vec2vec SphereDomain<dim>::GetExitTimes(const std::vector<Eigen::Matrix<numtype,dim,1> > &r,
                                        const std::vector<Eigen::Matrix<numtype,dim,1> > &s)
{
    // Returns the forward and backward exit times for a vector 
    // of phase space points (r,s)
    // tau[i](0) = forward, tau[i](1) = backward
    vec2vec tau;
    for(int i=0;i<r.size();++i){
        tau.push_back(GetExitTimes(r[i],s[i]));
    }
    return tau;
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
        vec2    GetExitTimes(const Eigen::Matrix<numtype,d,1> &r,
                             const Eigen::Matrix<numtype,d,1> &s) const; //Defined below
        vec2vec GetExitTimes(const std::vector<Eigen::Matrix<numtype,d,1> > &r,
                             const std::vector<Eigen::Matrix<numtype,d,1> > &s) const; //Defined below
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

// Exit time implementation (Vector input)
template<int dim>
vec2vec RectDomain<dim>::GetExitTimes(const std::vector<Eigen::Matrix<numtype,dim,1> > &r,
                                      const std::vector<Eigen::Matrix<numtype,dim,1> > &s) const {
    vec2vec tau;
    for(int i=0;i<r.size();++i){
        tau.push_back(GetExitTimes(r[i],s[i]));
    }
    return tau;
}

// Exit time implementation (3D)
template<int dim>
vec2 RectDomain<dim>::GetExitTimes(const Eigen::Matrix<numtype,dim,1> &r,
                                   const Eigen::Matrix<numtype,dim,1> &s) const {
    
    vec2 tau;  
    numtype inv_dir[dim];
    numtype tmin[dim];
    numtype tmax[dim];
    int sign[dim];
    for(int i=0;i<dim;++i){
        inv_dir[i] = -1./s[i]; // Negative, since computing "reverse" ray intersection
        sign[i] = (inv_dir[i] < 0);
    }
    
    for(int i=0;i<dim;++i){
        tmin[i] = (L(i,sign[i]) - r[i]) * inv_dir[i];
        tmax[i] = (L(i,(1-sign[i])) - r[i]) * inv_dir[i];
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
    
    tau(0) = tmin[0];
    tau(1) = tmax[0];
    return tau;
}


}// namespace CPU

#endif