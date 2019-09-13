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

#include "compute/linalg/linear_algebra_cpu.h"
#include "objects/lumpy/lumpy_bgnd_cpu.h"

namespace CPU{
struct LumpyXi 
{
    LumpyBgnd<double>* XiFun;
    LumpyXi(){ 
        XiFun = new LumpyBgnd<double>;
    }
    template<int dim>
    complexnumtype operator()(Eigen::Matrix<numtype,dim,1> r,
                              Eigen::Matrix<numtype,dim,1> s,
                              double e){
        return complexnumtype((*XiFun)(r.data()));
    }
};

struct LumpyMu
{
    LumpyBgnd<double>* MuFun;
    LumpyMu(){ 
        MuFun = new LumpyBgnd<double>;
    }
    template<int dim>
    complexnumtype operator()(Eigen::Matrix<numtype,dim,1> r,
                              Eigen::Matrix<numtype,dim,1> s,
                              double e){
        return complexnumtype((*MuFun)(r.data()));
    }
};

}//namespace CPU

#endif