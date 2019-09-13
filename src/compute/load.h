#ifndef __LOADSAMPLES_H__
#define __LOADSAMPLES_H__

#include"util/util.h"

template<typename T>
void loadSamples(T*** A, char* SDIR, char* FPRE, int s0, int sN, int* nPerGPU, int* sDim, int NGPU)
{
	printf("Loading samples from directory %s\n",SDIR);
	*A = (T**) malloc(NGPU*sizeof(T*));
	//This will create a sample matrix on NGPU devices
	int nSamp = (sN-s0)+1;
	char str[80]; //File name buffer
	if(nSamp % NGPU != 0){
		fprintf(stderr,"ERROR: number of columns must be divisible by NGPU\n");
		exit(1);
	}
	*nPerGPU = nSamp/NGPU; //Assumes NSamp is divisible by NGPU.
	sprintf(str,"%s%s%i.csv",SDIR,FPRE,s0);
	*sDim = GetNumLines(str);

	//Malloc memory 
	for(int igpu=0;igpu<NGPU;++igpu){
		CHECK(cudaMallocManaged(&((*A)[igpu]),(*nPerGPU)*(*sDim)*sizeof(T)));
	}

	for(int igpu=0;igpu<NGPU;++igpu){
		for(int j=0;j<*nPerGPU;++j){
			sprintf(str,"%s%s%i.csv",SDIR,FPRE,igpu*(*nPerGPU) + j);
			LoadFile<T>((*A)[igpu]+j*(*sDim),*sDim,str);
		}
	}
}


template<typename T>
void loadSamplesMat(T*** A, char* SDIR, char* FPRE, char* VARNAME, int s0, int sN, int* nPerGPU, int* sDim, int NGPU)
{
	printf("Loading samples from directory %s\n",SDIR);
	*A = (T**) malloc(NGPU*sizeof(T*));
	//This will create a sample matrix on NGPU devices
	int nSamp = (sN-s0)+1;
	char str[80]; //File name buffer
	if(nSamp % NGPU != 0){
		fprintf(stderr,"ERROR: number of columns must be divisible by NGPU\n");
		exit(1);
	}
	*nPerGPU = nSamp/NGPU; //Assumes NSamp is divisible by NGPU.
	sprintf(str,"%s%s%i.mat",SDIR,FPRE,s0);
	printf(str);printf("\n");
	*sDim = MatReadNumElem<T>(str,VARNAME);
	if(*sDim==0){
		std::cout<<"ERROR in loadSamplesMat"<<std::endl;
		return;
	}

	//Malloc memory 
	for(int igpu=0;igpu<NGPU;++igpu){
		CHECK(cudaMallocManaged(&((*A)[igpu]),(*nPerGPU)*(*sDim)*sizeof(T)));
	}

	for(int igpu=0;igpu<NGPU;++igpu){
		for(int j=0;j<*nPerGPU;++j){
			sprintf(str,"%s%s%i.mat",SDIR,FPRE,igpu*(*nPerGPU) + j);
			MatRead<T>((*A)[igpu]+j*(*sDim),*sDim,str,VARNAME);
		}
	}
}




// ****************SINGLE GPU IMPLEMENTATION NOT NECESSARY ****************//
/*
template<typename T>
void loadSamplesSingleGPU(T** A, char* SDIR, char* FPRE, int s0, int sN, int* sDim)
{
	int nSamp = (sN-s0)+1;
	char str[80]; //File name buffer
	// Get sample dimension from the first file, assuming all files same dim
	sprintf(str,"%s%s%i.csv",SDIR,FPRE,s0);
	*sDim = GetNumLines(str);
	CHECK(cudaMallocManaged(A,nSamp*(*sDim)*sizeof(T)));

	for(int j=0;j<nSamp;++j){
			sprintf(str,"%s%s%i.csv",SDIR,FPRE,j);
			LoadFile<T>(*A+j*(*sDim),*sDim,str);
	}
}
*/

#endif /* __LOADSAMPLES_H__ */
