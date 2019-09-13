// ---------------------------------------------------------------------------
// **********************CGRI SIMULATION TOOLBOX***************************
// file: rte_lumpy_mex_cpu.cc
// purpose: Contains the mex interface to compute ballistic RTE solutions
//          for lumpy-type functions.
// Author:  Nick Henscheid
// Date:    3-2019
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without technical
// support, and with no warranty, express or implied, as to its usefulness for
// any purpose.
// ---------------------------------------------------------------------------
#include <mex.h>
#include <iostream>
#include "imaging/rte/cpu/rte_cpu.h"     // RTE Class file 
#include "imaging/rte/cpu/domain_cpu.h"  // Domain definitions
#include "rte_parameters_cpu.h"          // Custom parameters

using namespace CPU;

template<typename T>
void getLumpyBgnd(LumpyBgnd<T>* Lbg, int dim,const mxArray *input)
{   
    // need to validate inputs!
    Lbg->set_dim(dim);
    Lbg->set_K((int)mxGetScalar(mxGetField(input,0, "K")));
    Lbg->set_B0((double)mxGetScalar(mxGetField(input,0, "B0")));
    Lbg->set_b0((double*)mxGetData(mxGetField(input,0, "b0")));
    Lbg->set_precision((double*)mxGetData(mxGetField(input,0, "cov")));
    Lbg->set_centers((double*)mxGetData(mxGetField(input,0, "centers")));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    // Gather inputs
    int dim = (int)mxGetScalar(prhs[0]);
    if((dim!=2)&&(dim!=3)){
        mexErrMsgTxt("Error: Invalid Dimension!");
    }
    // Source function
    // Set up evaluation function (defined in rte_parameters.h)
    // Xi is in slot 1, Mu is in slot 2.
    LumpyXi xi;
    getLumpyBgnd(xi.XiFun,dim,prhs[1]);
    
    // Attenuation function
    LumpyMu mu;
    getLumpyBgnd(mu.MuFun,dim,prhs[2]);
    
    // Evaluation points
    double* R  = (double*)mxGetData(prhs[3]);  // positions  (2xneval, col major)
    double* S  = (double*)mxGetData(prhs[4]);  // directions (2xneval, col major)
    int neval  = (int)mxGetScalar(prhs[5]);    // Number of evaluation points
    double cm  = (double)mxGetScalar(prhs[6]); // Speed of light in the medium
    int nstep  = (int)mxGetScalar(prhs[7]);    
  
    // Prepare output
    plhs[0]    = mxCreateDoubleMatrix(neval, 1, mxCOMPLEX);
    double* Wr  = (double*)mxGetPr(plhs[0]);
    double* Wi  = (double*)mxGetPi(plhs[0]);

    // Initialize the RTE integrator
    RK4Integrator Int;
    Int.SetStep(nstep);  // Number of steps for integrator.  Should make adaptive...
    // Set up input vectors and domain
    numvec  e;  //Not used ATM
    switch(dim){ // Doing it this way because dimension is a template parameter for the domain,
                 // but a runtime variable for RTE. 
        case 2:{
            mexPrintf("dim = 2\n");
            const vec2vec r = ArrayToVec<2>(R,neval);
            const vec2vec s = ArrayToVec<2>(S,neval);
            double* L2 = (double*)mxGetData(prhs[8]);
            mat2x2  L; L << L2[0],L2[1],L2[2],L2[3];
            RectDomain<2> D(L); // In theory, should make this more flexible...
            RTEBallistic<RectDomain<2>,LumpyMu,RK4Integrator> w(mu,D,Int);
            w.SetLightSpeed(cm);
            complexnumvec result = w(xi,r,s,e);
            for(int i=0;i<neval;++i){
                Wr[i] = result[i].real();
                Wi[i] = result[i].imag();
            }
            break;
        }
        case 3:{
            mexPrintf("dim = 3\n");
            const vec3vec r = ArrayToVec<3>(R,neval);
            const vec3vec s = ArrayToVec<3>(S,neval);
            double* L3 = (double*)mxGetData(prhs[8]);
            mat3x2  L; L << L3[0],L3[1],L3[2],L3[3],L3[4],L3[5];
            RectDomain<3> D(L);
            RTEBallistic<RectDomain<3>,LumpyMu,RK4Integrator> w(mu,D,Int);
            w.SetLightSpeed(cm);
            complexnumvec result = w(xi,r,s,e);
            for(int i=0;i<neval;++i){
                Wr[i] = result[i].real();
                Wi[i] = result[i].imag();
            }
            break;
        }
    }
    
    
}

