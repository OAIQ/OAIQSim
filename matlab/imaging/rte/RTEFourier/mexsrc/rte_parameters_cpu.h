// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_parameters_cpu.h
// purpose: This file contains defintions of helper functions used in the 
// RTE solver
// 
// Author:  Nick Henscheid
// Date:    5-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// INSTRUCTIONS: 
//     To compile a new RTE mex file, you must write your own helper functions.
//     Use these ones as templates. 
//     Each helper function must have at the very least a parenthesis operator defined 
//     i.e. double operator()(vec2 r, vec2 s, double e)
// ---------------------------------------------------------------------------
#ifndef __CGRI_IMAGING_RTE_PARAMETERS_CPU_H__
#define __CGRI_IMAGING_RTE_PARAMETERS_CPU_H__

#include "compute/linalg/linear_algebra_cpu.h"   //Some linear algebra things

struct XiFun 
{
    vec2 k;
    double operator()(vec2 r,vec2 s,double e){
    return cos(2*M_PI*k.dot(r))/(4.0*M_PI);
    }
};

// 2-dimensional complex Fourier mode with wavevector k
struct FourierMode2d
{
    vec2 k;  // Mode
    complexnumtype operator()(vec2 r,vec2 s,double e){
    return exp(2*M_PI*I*k.dot(r))/(4.0*M_PI);
    }
};

struct FourierMode3d
{
    vec3 k;  // Mode
    complexnumtype operator()(vec3 r,vec3 s,double e){
    return exp(2*M_PI*I*k.dot(r))/(4.0*M_PI);
    }
};

struct MuFun
{
    double var = 0.05;
    double mx  = 0.0;
    double my  = 0.0;
    double A   = 1.0;
    double operator()(vec2 r,vec2 s,double e){
        return A*exp(-(pow(r(0)-mx,2.0)+pow(r(1)-my,2.0))/(2*var));
    }
};

#endif