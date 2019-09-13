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
#include "compute/ode/ode_solvers_cpu.h"  // RTE Class file 

// This is a specialized mex interface to compute attenuated ray transforms
// of Fourier basis elements Phi_k(r) = exp(2*pi*i*k\cdot r)chi(r) where 
// chi(r) is the indicator function of a disc.

struct rhs
{
    complexnumtype lambda;
    complexnumtype operator()(double t, complexnumtype y){
        return lambda*y;
    }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    // This function solves y' = lambda*y using RK4 with nstep steps on 
    // the time interval (t[0],t[1])

    // Inputs: lambdar,lambdai: real and complex parts of lambda
    //         nstep: number of steps.  t: time interval.

    double lambdar = (double)mxGetScalar(prhs[0]);  
    double lambdai = (double)mxGetScalar(prhs[1]);
    int    nstep   = (int)mxGetScalar(prhs[2]);
    double* t      = (double*)mxGetData(prhs[3]);
    
    plhs[0]    = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    double* Wr  = (double*)mxGetPr(plhs[0]);
    double* Wi  = (double*)mxGetPi(plhs[0]);
  
    mexPrintf("(lambdar,lambdai,nstep,t[0],t[1]) = (%f,%f,%d,%i,%i)\n",lambdar,lambdai,nstep,t[0],t[1]);
    // Initialize the RK4 integrator
    RK4Integrator Integrate;
    Integrate.SetStep(nstep);
    
    std::vector<rhs> rhs_vec;
    rhs F; F.lambda = complexnumtype(lambdar,lambdai);
    rhs_vec.push_back(F);

    vec2 T(t[0],t[1]);
    vec2vec Tinterval; Tinterval.push_back(T);
    
    complexnumtype y0(1.0,0.0);
    complexnumvec  Y0; Y0.push_back(y0);

    complexnumvec sol = Integrate(rhs_vec,Tinterval,Y0);

    std::cout<<" F.lambda, T = "<<F.lambda<<","<<T<<std::endl;
    std::cout<<" sol[0] = "<<sol[0]<<std::endl;

    
    for(int i=0;i<sol.size();++i){
        Wr[i] = sol[i].real();
        Wi[i] = sol[i].imag();
        printf("(res[i].real(),res[i].imag()) = (%2.8e,%2.8e)\n",Wr[i],Wi[i]);
    }
}
