// -----------------------------------------------------------------------
// ************************ OAIQSIM TOOLBOX ******************************
// file: linear_algebra_cpu.h
// purpose: Contains various typedefs, includes and utility functions for 
//          linear algebra operations within the toolbox.
// Author:  Nick Henscheid
// Date:    5-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_COMPUTE_LIN_ALG_CPU_H__
#define __OAIQSIM_COMPUTE_LIN_ALG_CPU_H__

// Linear algebra includes
#include "compute/util_cpu.h"
#include <compute/extern/eigen3/Eigen/Dense>
//#include <vector>
//#include <complex>

// Vector typedefs

namespace CPU{
    typedef Eigen::Matrix<numtype,2,1> vec2;
    typedef Eigen::Matrix<numtype,3,1> vec3;
    typedef Eigen::Matrix<numtype,4,1> vec4;
    typedef Eigen::Matrix<complexnumtype,2,1> complexvec2;
    typedef Eigen::Matrix<complexnumtype,3,1> complexvec3;
    typedef Eigen::Matrix<complexnumtype,4,1> complexvec4;
    typedef std::vector<vec2> vec2vec;
    typedef std::vector<vec3> vec3vec;
    typedef std::vector<vec4> vec4vec;
    typedef std::vector<vec2> complexvec2vec;
    typedef std::vector<vec3> complexvec3vec;
    typedef std::vector<vec4> complexvec4vec;
    typedef std::vector<numtype> numvec; 
    typedef std::vector<complexnumtype> complexnumvec; 
    typedef Eigen::Matrix<numtype,2,2> mat2x2;
    typedef Eigen::Matrix<numtype,3,2> mat3x2;

// Utility functions 

//Function to convert a C-style array to an std::vector of Eigen::Vectorxd.
template<int d> //Dimension of vector
std::vector<Eigen::Matrix<numtype,d,1> > ArrayToVec(const numtype* x_in,int n)
{
    Eigen::Matrix<numtype,d,1> v_temp;
    std::vector<Eigen::Matrix<numtype,d,1> > x_out;
    for(int i=0;i<n;++i){
        for(int j=0;j<d;++j){
            v_temp(j) = x_in[d*i+j];
        }
        x_out.push_back(v_temp);
    }
    return x_out;
}
// Function to convert std::vector of Eigen::Vectorxd to C-style array
template<int d> //Dimension of vector
CPU::numtype* VecToArray(const std::vector<Eigen::Matrix<numtype,d,1> >& x_in)
{
    int n = x_in.size();
    numtype* x_out = (numtype*) malloc(n*d*sizeof(numtype));
    
    for(int i=0;i<n;++i){
        for(int j=0;j<d;++j){
            x_out[d*i+j] = x_in[i](j);
        }
    }
    return x_out;
}
// Function to convert a C-style array to an std::vector of numtypes
std::vector<numtype> ArrayToNumVec(const numtype* x_in,int n)
{
    std::vector<numtype> x_out;
    for(int i=0;i<n;++i){
        x_out.push_back(x_in[i]);
    }
    return x_out;
}
// Function to convert std::vector of numtypes to a C-style array
void NumVecToArray(const std::vector<numtype>& x_in,numtype* x_out)
{
    // Warning! Assumes that x_out is already allocated! bad design. fix.
    int n = x_in.size();
    for(int i=0;i<n;++i){
        x_out[i] = x_in[i];
    }
}


}
#endif