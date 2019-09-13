// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// File:    lumpy_kernel_cpu.h
// Purpose: This header contains the CPU compute kernel for the 
//          standard lumpy background 
// Author:  Nick Henscheid
// Date:    9-2016
// Contact: nph@email.arizona.edu
// References: "Factors Influencing Lesion Detection in Medical Imaging" 
//              J.P. Rolland, 1990 and 
//             "Foundations of Image Science", H.H. Barrett and K.J. Myers 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_OBJECTS_LUMPY_KERNEL_CPU_H__
#define __OAIQSIM_OBJECTS_LUMPY_KERNEL_CPU_H__
#include <math.h>

template<typename T,int d> 
double LumpyKernel(T* center,T* eval,T B0, T b0,T rb2)
{
    
}


template<typename T,int d>
void LumpyBgnd(T* centers,T* evalpts, T* result,int K,int nEval,T B0,T* b0,T* rb2)
{
    T exponent = 0;
    T power2 = 2.0; //Must convert to T else strange bugs appear
   
    T xEval,yEval,A;

    for(int ix=0;ix<nEval;++ix){
        xEval = evalpts[2*ix];
        yEval = evalpts[2*ix+1];
        result[ix] = B0;
        for(int iLump=0;iLump<K;++iLump){
            if((rb2[iLump]>0)&&(b0[iLump]>0)){
                A = (b0[iLump]/(M_PI*rb2[iLump])); 
                exponent = pow(xEval-centers[2*iLump],power2)+pow(yEval-centers[2*iLump+1],power2); 
                result[ix] += A*exp(-exponent/rb2[iLump]); 
            } else{
                //printf("Invalid input\n");
                result[ix] += 0.0;
            }
            
        }
    }
}

#endif