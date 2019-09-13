// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// File:    lumpy_kernel_cpu.h
// Purpose: This header contains the CPU compute kernel for the 
//          standard lumpy background 
// Author:  Nick Henscheid
// Date:    9-2016, 2-2019
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

namespace CPU{

template<typename T>
bool CheckInput(T*,int,int);

template<typename T> 
T LumpyKernel(T*,T*,int,int,T,T*,T*);

template<typename T>
class LumpyBgnd{
    public:
        // Default constructor
        LumpyBgnd(){
            dim = 2; K = 1;
            B0 = 0; 
            b0 = (T*)malloc(sizeof(T)); b0[0] = 1.0;
            rb2 = (T*)malloc(2*sizeof(T)); rb2[0] = 0.1;rb2[1] = 0.1;
            centers = (T*)malloc(2*sizeof(T)); centers[0] = 0.5;centers[1] = 0.5;
        };
        // Construct default object of dimension dim
        LumpyBgnd(int dim_i){
            if(dim_i>0){
                dim = dim_i; K = 1;
                B0 = 0; 
                b0 = (T*)malloc(sizeof(T)); b0[0] = 1.0;
                rb2 = (T*)malloc(dim*sizeof(T)); 
                centers = (T*)malloc(dim*sizeof(T));
                for(int i=0;i<dim;++i){
                    rb2[i] = 0.1;
                    centers[i] = 0.5;
                }                
            }
        };
        LumpyBgnd(int dim_i,int K_i,T* centers_i,T B0_i,T* b0_i,T* rb2_i){
            if((dim_i>0)&&(K_i>0)&&(CheckInput(rb2_i,K_i,dim_i))){
                dim = dim_i; K = K_i;
                centers = centers_i;
                B0 = B0_i; b0 = b0_i; rb2 = rb2_i;
            }
        }
        
        void set_dim(int d){dim = d;}
        int  dimension(void){return dim;}
        void set_K(int k){K=k;}
        int  num_lumps(void){return K;}
        void set_B0(numtype B){B0 = B;}
        numtype DC_val(void){return B0;}
        void set_centers(numtype* p){if(p){centers = p;}}// if(p){} checks for null ptr
        void set_b0(numtype* p){if(p){b0 = p;}}
        void set_rb2(numtype* p){if(p){rb2 = p;}}
        
        void eval(int nEval,T* evalpts,T* result)
        {
            for(int ix=0;ix<nEval;++ix){
                result[ix] = LumpyKernel<T>(centers,evalpts+dim*ix,dim,K,B0,b0,rb2);
            }
        }
        T operator()(T* evalpt){
            // Evaluate the texture for a single point
            return LumpyKernel(centers,evalpt,dim,K,B0,b0,rb2);
        }
        
    private:
        int dim;
        int K;
        T B0;
        T* centers;
        T* b0;
        T* rb2;
};// LumpyBgnd

template<typename T> 
T LumpyKernel(T* centers,T* eval,int dim,int K,T B0,T* b0,T* rb2)
{
    // Should check if eval satisfies a domain condition?
    T exponent = 0;
    T power2   = 2.0;  // Must convert to T.
    T result   = B0;
    T A;
    for(int iLump=0;iLump<K;++iLump){
        // Should figure out how to limit the loop to only those lumps that
        // are "close" to the evaluation point.  Pseudogrid?
        if(b0[iLump]!=0){
            A = b0[iLump]/pow(2.0*M_PI,((T)dim)/2.0); 
            exponent = 0.0;
            for(int i=0;i<dim;++i){
                exponent += pow(eval[i]-centers[dim*iLump+i],power2)/rb2[dim*iLump+i]; 
                A*=1/sqrt(rb2[dim*iLump+i]);
            }
            result += A*exp(-0.5*exponent); 
        } 
    }
    return result;
}

// Doing this check in the constructor prevents divide-by zero 
template<typename T>
bool CheckInput(T* rb2,int K,int dim)
{
    bool x = 1;
    for(int i=0;i<dim*K;++i){
        if(rb2[i]<=0){
            x = 0;
        }
    }
    return x;
   
}

}//namespace CPU
#endif