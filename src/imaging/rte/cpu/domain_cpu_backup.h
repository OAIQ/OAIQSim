// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// file: domain_cpu.h
// purpose: Contains the class definition for the disc domain component of  
//          the RTE simulation module
// Author:  Nick Henscheid
// Date:    9-2016, 11-2018
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

/*2D DISC DOMAIN */
class DiscDomain
{
    public:
        DiscDomain(double R_in){R = R_in;}
        DiscDomain(void){R = 1.0;}
        void UpdateRadius(double r){R = r;}
        // Exit time functions, defined below. 
        vec2 GetExitTimes(vec2 r,vec2 s); // Single eval version
        vec2vec GetExitTimes(vec2vec r,vec2vec s); // Multiple eval
        // Domain indicator function
        bool IsInside(vec2 r){
            return (r.dot(r)<R*R)? true:false;
        }
    private:
         double R;  //Radius of the circle
};
// Exit time implementations
vec2 DiscDomain::GetExitTimes(vec2 r,vec2 s)
{
    // Returns the forward and backward exit times for a single 
    // phase point (r,s)
    // tau(0) = forward, tau(1) = backward
    vec2 tau;
    double a = r.dot(s);
    tau(0) = fabs(a - sqrt(pow(a,2.0) - (r.dot(r)-pow(R,2.0))));
    tau(1) = fabs(a + sqrt(pow(a,2.0) - (r.dot(r)-pow(R,2.0))));
    return tau;
}
vec2vec DiscDomain::GetExitTimes(vec2vec r,vec2vec s)
{
    // Returns the forward and backward exit times for a vector 
    // of phase space points (r,s)
    // tau[i](0) = forward, tau[i](1) = backward
    vec2vec tau;
    vec2 tau_temp,r_temp,s_temp;
    double a;
    for(int i=0;i<r.size();++i){
        r_temp = r[i]; s_temp = s[i];
        a = r_temp.dot(s_temp);
        tau_temp(0) = fabs(a - sqrt(pow(a,2.0) - (r_temp.dot(r_temp)-pow(R,2.0))));
        tau_temp(1) = fabs(a + sqrt(pow(a,2.0) - (r_temp.dot(r_temp)-pow(R,2.0))));
        tau.push_back(tau_temp);
    }
    return tau;
}

/*2D RECTANGULAR DOMAIN */
template<int d>
class RectDomain
{
    public:
        RectDomain(Eigen::Matrix<numtype,d,2> L_in){
            L = L_in;
            for(int i=0;i<d;++i){
                center(i)= (L(i,1)-L(i,0))/2;
            }
        }
        RectDomain(void){
            for(int i=0;i<d;++i){
                L(i,0) = 0.0;
                L(i,1) = 1.0; 
                center(i) = 0.5;
            }
        }
        void updateDims(Eigen::Matrix<numtype,d,2> L_i){
            L = L_i;
            for(int i=0;i<d;++i){
                center(i)= (L(i,1)-L(i,0))/2;
            }
        }
        vec2    getExitTimes(const Eigen::Matrix<numtype,d,1> &r,
                             const Eigen::Matrix<numtype,d,1> &s) const; //Defined below
        vec2vec getExitTimes(const std::vector<Eigen::Matrix<numtype,d,1>> &r,
                             const std::vector<Eigen::Matrix<numtype,d,1>> &s) const; //Defined below
        bool isInside(Eigen::Matrix<numtype,d,1> r){
            bool x = 1;
            for(int i=0;i<d;++i){
                x*= (r(i)>=L(i,0))&&(r(i)<=L(i,1));
            }
            return x;
        }
    private:
        Eigen::Matrix<numtype,d,2> L;       // Rectangle dimensions
        Eigen::Matrix<numtype,d,2> center;  // Rectangle center
};

// Exit time implementation (2D.  Combine w/ 3D?)
template<>
vec2 RectDomain<2>::getExitTimes(const vec2 &r,
                                 const vec2 &s) const {
    vec2 tau;
    vec4 t;
    double tmin,tmax;
    // Check that the ray is at least pointing away from the boundary
    if( ((r(0)>L(0,1))&&(s(0)<=0))||
        ((r(0)<L(0,0))&&(s(0)>=0))||
        ((r(1)>L(1,1))&&(s(1)<=0))||
        ((r(1)<L(1,0))&&(s(1)>=0))){
        std::cout<<"Ray pointing in"<<std::endl;
        std::cout<<"r = "<<r<<"s = "<<s<<"xbounds = "<< L(0,0)<<","<<L(0,1)<<"ybounds = "<<L(1,0)<<","<<L(1,1)<<std::endl;
        tau << 0.0,0.0;
        return tau;
    }
    
    if(s(0)>=0){
        t(0) = (r(0)-L(0,1))/s(0);
        t(1) = (r(0)-L(0,0))/s(0);
    }else{
        t(0) = (r(0)-L(0,0))/s(0);
        t(1) = (r(0)-L(0,1))/s(0);
    }
    if(s(1)>=0){
        t(2) = (r(1)-L(1,1))/s(1);
        t(3) = (r(1)-L(1,0))/s(1);
    }else{
        t(2) = (r(1)-L(1,0))/s(1);
        t(3) = (r(1)-L(1,1))/s(1);
    }

    if(t(0)>t(2)){
        tmin = t(0);
    }else{
        tmin = t(2);
    }

    if(t(1)<t(3)){
        tmax = t(1);
    }else{
        tmax = t(3);
    }
    
    if(tmin>tmax){
        printf("No intersection!\n");
        tau << 0.0,0.0;
        return tau;
    }
    tau << tmin,tmax;
    return tau;
}

// Exit time implementation (2D, vector input  COMBINE W/ 3D??)
template<>
vec2vec RectDomain<2>::getExitTimes(const vec2vec &r,
                                    const vec2vec &s) const {
    vec2vec tau;
    vec4 t;
    for(int i=0;i<r.size();++i){
        double tmin,tmax;
        // Check that the ray is at least pointing away from the boundary
        if( ((r[i](0)>L(0,1))&&(s[i](0)<=0))||
            ((r[i](0)<L(0,0))&&(s[i](0)>=0))||
            ((r[i](1)>L(1,1))&&(s[i](1)<=0))||
            ((r[i](1)<L(1,0))&&(s[i](1)>=0))){
            std::cout<<"Ray pointing in"<<std::endl;
            std::cout<<"r = "<<r[i]<<"s = "<<s[i]<<"xbounds = "<< L(0,0)<<","<<L(0,1)<<"ybounds = "<<L(1,0)<<","<<L(1,1)<<std::endl;
            tau.push_back(vec2(0.0,0.0));
        }

        if(s[i](0)>=0){
            t(0) = (r[i](0)-L(0,1))/s[i](0);
            t(1) = (r[i](0)-L(0,0))/s[i](0);
        }else{
            t(0) = (r[i](0)-L(0,0))/s[i](0);
            t(1) = (r[i](0)-L(0,1))/s[i](0);
        }
        if(s[i](1)>=0){
            t(2) = (r[i](1)-L(1,1))/s[i](1);
            t(3) = (r[i](1)-L(1,0))/s[i](1);
        }else{
            t(2) = (r[i](1)-L(1,0))/s[i](1);
            t(3) = (r[i](1)-L(1,1))/s[i](1);
        }

        if(t(0)>t(2)){
            tmin = t(0);
        }else{
            tmin = t(2);
        }

        if(t(1)<t(3)){
            tmax = t(1);
        }else{
            tmax = t(3);
        }

        if(tmin>tmax){
            printf("No intersection!\n");
            tau.push_back(vec2(0.0,0.0));
        }
        tau.push_back(vec2(tmin,tmax));
    }
    return tau;
}

// Exit time implementation (3D)
template<>
vec2 RectDomain<3>::getExitTimes(const vec3 &r,
                                 const vec3 &s) const {
    vec2 tau(0.0,0.0);
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    vec3 inv_dir(-1/s(0),-1/s(1),-1/s(2));
    
    
    int sign[3];
    sign[0] = (inv_dir(0) < 0);
    sign[1] = (inv_dir(1) < 0);
    sign[2] = (inv_dir(2) < 0);
    //printf("r = (%f,%f,%f), inv_dir = (%f,%f,%f), sign = (%i,%i,%i)\n",r[0],r[1],r[2],inv_dir(0),inv_dir(1),inv_dir(2),sign[0],sign[1],sign[2]);
    tmin = (L(0,sign[0]) - r(0)) * inv_dir(0);
    tmax = (L(0,1-sign[0]) - r(0)) * inv_dir(0);
    if(isnan(tmin))
        tmin = -1.0/0.0;  // -inf
    if(isnan(tmax))
        tmax = 1.0/0.0;   // +inf
    //printf("(tmin,tmax) = (%f,%f)\n",tmin,tmax);
    tymin = (L(1,sign[1]) - r(1)) * inv_dir(1);
    tymax = (L(1,1-sign[1]) - r(1)) * inv_dir(1);
    if(isnan(tymin))
        tymin = -1.0/0.0; // -inf
    if(isnan(tymax))
        tymax = 1.0/0.0;  // +inf
    //printf("(tymin,tymax) = (%f,%f)\n",tymin,tymax);
    if ( (tmin > tymax) || (tymin > tmax) ) 
        return tau;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    tzmin = (L(2,sign[2]) - r(2)) * inv_dir(2);
    tzmax = (L(2,1-sign[2]) - r(2)) * inv_dir(2);
    if(isnan(tzmin))
        tzmin = -1.0/0.0; // -inf
    if(isnan(tzmax))
        tzmax = 1.0/0.0;  // +inf
    
    //printf("(tzmin,tzmax) = (%f,%f)\n",tzmin,tzmax);
    
    if ( (tmin > tzmax) || (tzmin > tmax) ) 
        return tau;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    
    //printf("(tmin,tmax) = (%f,%f)\n",tmin,tmax);
    tau(0) = tmin; 
    tau(1) = tmax; 
    return tau; 
}

// Exit time implementation (3D, vector input)
template<>
vec2vec RectDomain<3>::getExitTimes(const vec3vec &r,
                                    const vec3vec &s) const {
    vec2vec tau;
    for(int i=0;i<r.size();++i){
        vec2 taui = getExitTimes(r[i],s[i]);
        tau.push_back(taui);
    }
    return tau; 
}


}// namespace CPU

#endif