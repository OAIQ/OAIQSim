// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_fourier_cpu_mex.h
// purpose: Contains the mex interface to compute ballistic RTE solutions  
// Author:  Nick Henscheid
// Date:    4-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include <mex.h>
#include <iostream>
#include "imaging/rte/cpu/rte_cpu.h"  // RTE Class file 
#include "imaging/rte/cpu/domain_cpu.h"  // Domain definitions
#include "rte_parameters_cpu.h"  // Custom parameters

// This is a specialized mex interface to compute attenuated ray transforms
// of Fourier basis elements Phi_k(r) = exp(2*pi*i*k\cdot r)chi(r) where 
// chi(r) is the indicator function of a disc.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    // Gather inputs
    double* R  = (double*)mxGetData(prhs[0]);  // positions  (2xneval, col major)
    double* S  = (double*)mxGetData(prhs[1]);  // directions (2xneval, col major)
    double* k  = (double*)mxGetData(prhs[2]);  // Wavenumber 
    int neval  = (int)mxGetScalar(prhs[3]);    // Number of evaluation points
    double cm  = (double)mxGetScalar(prhs[4]); // Speed of light in the medium
    int nstep  = (int)mxGetScalar(prhs[5]);    
    // Prepare output
    plhs[0]    = mxCreateDoubleMatrix(neval, 1, mxCOMPLEX);
    double* Wr  = (double*)mxGetPr(plhs[0]);
    double* Wi  = (double*)mxGetPi(plhs[0]);

    // Initialize RTE helper functions (from domain_cpu.h and rte_parameters.h)
    DiscDomain D; //Standard unit disc domain
    MuFun mu; mu.A = 0.0; //No attenuation
    // Initialize the RTE integrator
    RK4Integrator Int;
    Int.SetStep(nstep);
    // Initialize the ballistic RTE solver
    RTEBallistic<DiscDomain,MuFun,RK4Integrator> w(mu,D,Int);
    w.SetLightSpeed(cm);
   
    // Set up input vectors 
    vec2vec r = ArrayToVec<2>(R,neval);
    vec2vec s = ArrayToVec<2>(S,neval);
    numvec  e;

    for(int i=0;i<neval;++i){
        e.push_back(0.0);
    }
    
    // Set up evaluation function (defined in rte_parameters.h)
    FourierMode2d xi;
    xi.k << k[0],k[1];
    //std::cout<<"xi.k = "<<xi.k<<std::endl;
    //std::cout<<"xi(r0) = "<<xi(r[0],s[0],e[0])<<std::endl;
    complexnumvec result = w(xi,r,s,e);
    
    for(int i=0;i<result.size();++i){
        Wr[i] = result[i].real();
        Wi[i] = result[i].imag();
        //printf("(w_res[i].real(),w_res[i].imag()) = (%2.3e,%2.3e)\n",Wr[i],Wi[i]);
    }
}

