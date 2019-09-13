// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_fourier_cpu_mex.h
// purpose: Contains the mex interface to compute ballistic RTE solutions  
// Author:  Nick Henscheid
// Date:    4-2017, 3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:   
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include <mex.h>
#include <iostream>
#include "imaging/rte/gpu/rte_gpu.cuh"  // RTE Class file 
#include "imaging/rte/gpu/domain_gpu.cuh"  // Domain definitions
#include "rte_parameters_gpu.cuh"  // Custom parameters
#include "compute/util_gpu.cuh"

// This is a specialized mex interface to compute attenuated ray transforms of Fourier basis elements 
// Phi_k(r) = exp(2*pi*i*k\cdot r)chi(r) 

using namespace GPU;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    int dim    = (int)mxGetScalar(prhs[0]);    // Dimension (2 or 3)
    if((dim!=2)&&(dim!=3)){
        mexErrMsgTxt("Error! Only 2D and 3D!\n");
    }
    
    numtype* R  = (numtype*)mxGetData(prhs[1]);  // positions  (2xneval, col major)
    numtype* S  = (numtype*)mxGetData(prhs[2]);  // directions (2xneval, col major)
    numtype* k  = (numtype*)mxGetData(prhs[3]);  // Wavenumber 
    int neval  = (int)mxGetScalar(prhs[4]);    // Number of evaluation points
    numtype cm  = (numtype)mxGetScalar(prhs[5]); // Speed of light in the medium
    numtype hmin  = (numtype)mxGetScalar(prhs[6]);

    plhs[0]    = mxCreateNumericMatrix(neval, 1, mxDOUBLE_CLASS,mxCOMPLEX);
    numtype* Wr  = (numtype*)mxGetPr(plhs[0]);
    numtype* Wi  = (numtype*)mxGetPi(plhs[0]);
    
    numvec e(neval,0.0);  // Create a vector of energies (all set to zero!)
    
    complexnumvec result;
    if(dim==2){
        vec2vec r; // Convert to host vector-of-2D-vectors
        vec2vec s; 
        for(int i=0;i<neval;++i){
            r.push_back(vec2(R[2*i],R[2*i+1]));
            s.push_back(vec2(S[2*i],S[2*i+1]));
        }
        SphereDomain<2>* D = new SphereDomain<2>();
        MuFun<2>* mu = new MuFun<2>();
        XiFun<2>* xi = new XiFun<2>(k);
        RTEBallistic<2,SphereDomain<2>,MuFun<2>,XiFun<2> > W(mu,xi,D);
        W.setLightSpeed(cm);
        W.setStep(hmin);
        result = W(r,s,e);
        for(int i=0;i<neval;++i){
            Wr[i] = result[i].real();
            Wi[i] = result[i].imag();
        }
    }else{ //Already checked that dim==2 or dim==3
        vec3vec r,s; // Convert to host vector-of-3D-vectors
        for(int i=0;i<neval;++i){
            r.push_back(vec3(R[3*i],R[3*i+1],R[3*i+2]));
            s.push_back(vec3(S[3*i],S[3*i+1],S[3*i+2]));
        }
        SphereDomain<3>* D = new SphereDomain<3>();
        MuFun<3>* mu = new MuFun<3>();
        XiFun<3>* xi = new XiFun<3>(k);
        RTEBallistic<3,SphereDomain<3>,MuFun<3>,XiFun<3> > W(mu,xi,D);
        W.setLightSpeed(cm);
        W.setStep(hmin);
        result = W(r,s,e);
        for(int i=0;i<neval;++i){
            Wr[i] = result[i].real();
            Wi[i] = result[i].imag();
        }
    } 
}

