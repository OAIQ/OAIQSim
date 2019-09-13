#ifndef __UTIL_GPU_MATLAB_H__
#define __UTIL_GPU_MATLAB_H__

#include <compute/util_gpu.cuh>
//#include <mat.h>
//#include <mex.h>

#define MEXCHECK(call)						   							    \
{																			\
	const cudaError_t error = call;											\
	if(error != cudaSuccess)												\
	{																		\
		printf("Error: %s:%d, ",__FILE__,__LINE__);							\
		printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));	\
		mexErrMsgTxt("Returning to Matlab\n"); 								\
	}																		\
}	

template<typename T>
inline int MatReadNumElem(const char* file,const char* var)
{
	mwSize num = 0;
    MATFile* pmat = matOpen(file, "r");
    if (pmat == NULL) {
        std::cout<<file<<" not found"<<std::endl;
        return 0;
    }
    // extract the specified variable
    //mxArray* struc = matGetVariable(pmat, "u");
    mxArray* arr = matGetVariable(pmat,var);
    if(arr == NULL){
        std::cout<<"NULL VARIABLE"<<std::endl;
        return 0;
    }

    if(mxIsEmpty(arr)){
        std::cout<<"Empty!"<<std::endl;
		return 0;
	}
    if (arr != NULL && mxIsSingle(arr) && !mxIsEmpty(arr)) {
        // copy data
        num = mxGetNumberOfElements(arr);
		
	}
    mxDestroyArray(arr);
    matClose(pmat);
	return (int)num;   
}


// Read in a matlab .mat format file
template<typename T>
inline void MatRead(T* x,int N,const char* file,const char* var)
{
    std::cout<<"Reading in "<<file<<std::endl;
    // open MAT-file
    MATFile* pmat = matOpen(file, "r");
    if (pmat == NULL) {
        std::cout<<file<< " not found"<<std::endl;
        return;
    }
    // extract the specified variable
    //mxArray* struc = matGetVariable(pmat, "u");
    mxArray* arr = matGetVariable(pmat,var);
    if(arr == NULL){
        std::cout<<"NULL VARIABLE"<<std::endl;
        return;
    }

    if(mxIsEmpty(arr)){
        std::cout<<"Empty!"<<std::endl;
		return;
	}
    if (arr != NULL && mxIsSingle(arr) && !mxIsEmpty(arr)) {
        // copy data
        mwSize num = mxGetNumberOfElements(arr);
		if(num != N){
			std::cout<<"Input data and array sizes do not match!"<<std::endl;
			return;
		}
		
        float *pr = (float*)mxGetData(arr);
        if (pr != NULL) {
            for(int i=0;i<N;++i){
            	x[i] = pr[i];
            }
        }
    }

    // cleanup
	//mxDestroyArray(struc);
    mxDestroyArray(arr);
    matClose(pmat);
}



#endif /* !UTIL_GPU_MATLAB_CUH */

