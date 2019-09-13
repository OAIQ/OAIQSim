#ifndef __UTIL_H__
#define __UTIL_H__

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>    // For std::string
#include<string.h>  // For strcmp, etc.
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include "mat.h"
#include<complex>  // std::complex (superseded by thrust/complex on GPU)
#include<vector>   // std::vector<T>

// Namespace is to avoid typedef conflicts with GPU code.
namespace CPU{
    typedef double numtype;  // Should probably put this in a more central location
    typedef std::complex<numtype> complexnumtype;
    const   complexnumtype I(0.0,1.0);   // imaginary unit.  Allows for "1.0+1.0*I" notation.
}

//Macro for "verbose mode"  //WHERE IS THIS USED...?
#define VERB(call) if(VMODE){call;}

//Macro for cublas style array indexing (column major)
#define IDX2C(i,j,ld)((j)*(ld)+(i))
                                   
// ***************************************
// *************  File I/O ***************
// ***************************************
// Get number of lines in a txt file
inline int GetNumLines(const char* fname){
	int N = 0;
	std::string line;
	std::ifstream file(fname);
	if(!file){ 
		std::cout<<"Cannot open file.\n";
		return 0;
	}
	while(std::getline(file,line)){
		++N;
	}
	return N;
}

// Load file into array of T's
template<typename T>
inline void LoadFile(T* x, int N, const char* fname) {
	int i;
	std::ifstream in(fname);

	if (!in) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

	for(i=0;i<N;++i){
		in >> x[i];
	}
	
  	in.close();
}

// Save array of T's to a file
template<typename T>
inline void SaveFile(T* x, int N, const char* fname) {
	int i;
    int precision = 16; // Double precision saving
	std::ofstream out;
	out.open(fname);
	out.precision(precision);
	if (!out) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

	for(i=0;i<N;++i){
		out << x[i];
		out << "\n";
	}
	
  	out.close();
}

// Save N arrays of T's to a csv
// The array is assumed to be unrolled i.e. x = [x1;x2;...;xM]
// Column headings are contained in the std::vector<string> cnames
// ASSUMING THAT EACH COLUMN HAS SAME NUMBER OF ROWS

template<typename T>
inline void SaveCSV(T* x,int M, int N, const char* fname, std::vector<std::string> cnames) {
	int i,j;
    int precision = 16; // Double precision saving
	std::ofstream out;
	out.open(fname);
	out.precision(precision);
	if (!out) {
    	std::cout << "Cannot open file.\n";
    	return;
	}

    for(int icol = 0;icol<N-1;++icol){
        out << cnames[icol] << ",";
    }
    out<< cnames[N-1] <<"\n";
	for(i=0;i<M;++i){ //rows
        for(j=0;j<N-1;++j){ //cols
    		out << x[j*M + i]<<",";
        }
        out << x[(N-1)*M + i] << "\n";
	}
	
  	out.close();
}

// ***************************************
// **********  Matlab File I/O ***********
// ***************************************

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

// ***************************************
// ***********  Misc. Math Ops ***********
// ***************************************
//Compute the p norm of an array. T should be a floating point type, e.g. float or double.
//This should get replaced with something better.
template<typename T>
inline T PNorm(T* x,int N,T p){
	T result = 0;
	for(int i=0;i<N;++i){
		result += pow(abs(x[i]),p);
	}
	result = pow(result,1/p);
	return result;
}
//Compute p distance of two arrays 
template<typename T>
inline T PDist(T* x,T* y,int N,T p){
	T result = 0;
	for(int i=0;i<N;++i){
		result += pow(abs(x[i]-y[i]),p);
	}
	result = pow(result,1/p);
	return result;
}


#endif /* !UTIL_CUH */

