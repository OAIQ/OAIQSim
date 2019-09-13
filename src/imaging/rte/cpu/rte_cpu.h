// -----------------------------------------------------------------------
// ************************ OAIQSIM TOOLBOX ******************************
// File:    rte_cpu.h
// Purpose: Contains the class definitions for the CPU version of the  
//          ballistic RTE solver 
// Author:  Nick Henscheid
// Date:    4-2017
// Contact: nhenscheid@math.arizona.edu
// References: 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_IMAGING_RTE_CPU_H__
#define __OAIQSIM_IMAGING_RTE_CPU_H__

#include "compute/ode/ode_solvers_cpu.h"

namespace CPU{
template<class Domain,class MuFun,class Integrator>
class RTEBallistic
{
    public:
        // This is the ode right hand side that will get passed to the ode integrator
        template<typename XiFun,int dim>
        struct RTERHS
        {
            numtype tau;  // Backward exit time
            Eigen::Matrix<numtype,dim,1>   r0;   // Initial point
            Eigen::Matrix<numtype,dim,1>   s0;   // Initial direction
            numtype e0;   // Energy
            numtype cm;   // Speed of light
            MuFun *mu;
            XiFun *xi;
            RTERHS(MuFun *mu_in,XiFun *xi_in,numtype cm_in){
                mu = mu_in;xi = xi_in;cm = cm_in;
            }
            RTERHS(){};
            complexnumtype operator()(numtype t,complexnumtype y){
                Eigen::Matrix<numtype,dim,1> temp = r0 + (t-tau)*s0;
                return (-(*mu)(temp,s0,e0))*y + ((*xi)(temp,s0,e0))/cm;
            }
        };
        // Constructor 
        RTEBallistic(MuFun mu_in, Domain D_in, Integrator I_in){
            mu = mu_in;D = D_in;Int = I_in;
        }
        //void SetStep(int n){nstep = n;}
        void SetLightSpeed(numtype c){cm = c;}
        // Evaluation function.  Allows for simple w(r,s,e) evaluation. 
        template<class argtype,int dim>
        complexnumvec operator()(argtype xi,const std::vector<Eigen::Matrix<numtype,dim,1>>& R,
                                            const std::vector<Eigen::Matrix<numtype,dim,1>>& S,
                                            const numvec& E){
            vec2vec tau = D.GetExitTimes(R,S); // Get exit times (need a version of getexittimes that accepts vectors...)
            vec2vec times;
            complexnumvec  initvals;
            std::vector<RTERHS<argtype,dim> > rhs_vec;
            RTERHS<argtype,dim> F(&mu,&xi,cm);
            for(int i=0;i<R.size();++i){
                F.r0  = R[i];
                F.s0  = S[i];
                F.tau = tau[i](1);
                rhs_vec.push_back(F);
                if(isnan(tau[i](1)))
                    printf("NaN time detected!\n");
                times.push_back(vec2(0.0,tau[i](1)));
                initvals.push_back(0.0+0.0*I);
            }
            complexnumvec w = Int(rhs_vec,times,initvals);  // Solve the odes
            return w;
        }
        
    private: 
        MuFun mu;     
        Domain D; 
        Integrator Int;
        numtype cm = 299792458; // Unless otherwise specified     
};

} //namespace CPU
#endif