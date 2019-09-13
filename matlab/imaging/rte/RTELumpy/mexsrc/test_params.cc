#include "compute/util_cpu.h"
#include "compute/linalg/linear_algebra_cpu.h"

#include "rte_parameters_cpu.h"


int main()
{
    


    LumpyXi Xi;
    LumpyMu Mu;
    
    CPU::vec2 r(0.25,0.5);
    CPU::vec2 s(0.0,0.0);
    std::cout<<r<<std::endl;
    std::cout<<s<<std::endl;
    
    
    double* r_a = r.data();
    
    std::cout<<r_a[0]<<"\n"<<r_a[1]<<std::endl;
    
    printf("Xi.L.centers = (%f,%f)\n",Xi.Xi.centers[0],Xi.Xi.centers[1]);
    
    printf("Xi(0,0,0) = %f\n",Xi(r,s,0.0));
    printf("Mu(0,0,0) = %f\n",Mu(r,s,0.0));



    return 0;
}