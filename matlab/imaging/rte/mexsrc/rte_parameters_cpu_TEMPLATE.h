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
// INSTRUCTIONS 
//   To compile a new RTE solver, 
// ---------------------------------------------------------------------------
#ifndef __CGRI_IMAGING_RTE_PARAMETERS_H__
#define __CGRI_IMAGING_RTE_PARAMETERS_H__

struct XiFun //: public XiFun
{
    vec2 k;
    XiFun(double k1,double k2){k = vec2(k1,k2);}
    XiFun(){};
    double operator()(vec2 r,vec2 s,double e){
        return cos(2*M_PI*k.dot(r))/(4.0*M_PI);
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