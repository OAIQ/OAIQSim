// ------------------------------------------------------------------------
// *************************** OAIQSIM TOOLBOX ****************************
// file: ode_solvers_cpu.h
// purpose: Contains various ODE solvers.  Used mainly by the RTE 
//          solver in the imaging module.
// Author:  Nick Henscheid
// Date:    4-2017, 10-2018
// Contact: nhenscheid@math.arizona.edu
// References: 
// Notes:   - For RK4Kernel, only the result at the final time is returned.
//          - 
// This software is in the public domain, furnished "as is", without 
// technical support, and with no warranty, express or implied, as to its 
// usefulness for any purpose.
// ------------------------------------------------------------------------
#ifndef __OAIQSIM_COMPUTE_ODE_SOLVERS_CPU_H__
#define __OAIQSIM_COMPUTE_ODE_SOLVERS_CPU_H__

#include "compute/linalg/linear_algebra_cpu.h"


namespace CPU{
template<typename rhs>   
void RK4Kernel(std::vector<rhs>& F,CPU::complexnumvec& result,vec2vec& Tinterval,
               CPU::complexnumvec& y0,int n_step)
{
    // Solves an ensemble of scalar ODEs: one for each rhs in the vector F,  
    // init cond in y0, and time interval [t_i,t_f] in Tinterval.  
    // !!NOTE!! Only result at t_f is returned !!  
    CPU::complexnumtype y;
    numtype t,h;
    CPU::complexnumtype k1,k2,k3,k4;
    for(int j=0;j<result.size();++j){
        rhs f = F[j];
        t = Tinterval[j](0);
        y = y0[j];
        h = (Tinterval[j](1)-Tinterval[j](0))/((numtype)n_step);
        for(int i=0;i<n_step;++i) 
        {
            k1 = f(t,y);
            k2 = f(t+ h/2.0,y + h*k1/2.0);
            k3 = f(t+ h/2.0,y + h*k2/2.0);
            k4 = f(t+ h,y + h*k3);
            y = y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
            t = t + h; 
        }
        result[j] = y;
    }
}

    class RK4Integrator{
        public:
            RK4Integrator(void){n_step = 1000;}
            void SetStep(int n){n_step = n;}

            template<typename rhs>
            complexnumvec operator()(std::vector<rhs>& rhs_vec,vec2vec& Tinterval,complexnumvec& y0){
                int N = rhs_vec.size();
                complexnumvec result(N);
                RK4Kernel(rhs_vec,result,Tinterval,y0,n_step);
                return result;
            }
        private:
            int n_step;
    };
}//namespace CPU

#endif 