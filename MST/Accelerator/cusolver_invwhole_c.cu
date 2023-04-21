#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <complex.h>
//#include "cudaDoubleComplex.hpp"
//#include "DeviceStorage.hpp"

extern "C"
void cusolver_invwhole_c_(int *m, double _Complex *a, int *block_size, double _Complex *b)
{
static bool initialized = false;
static cuDoubleComplex  *aDev;
static cuDoubleComplex  *workArray;
static cuDoubleComplex  *tau00Dev;
static cudaError_t error;
static int *pivotArray;
static int *infoArray;
static int Lwork;
float time_copyin=0;
float time_copyout=0;
float time_compute=0;

// LFu: fix segmentation fault when allocating the array
double _Complex *tau00;
tau00 = (double _Complex *) malloc(sizeof(double _Complex) * *m * *m);

for (int i=0;i<*m;i++){
for (int j=0;j<*m;j++){
	if (i==j){tau00[i * *m + j]={1,0};}
        else{tau00[i * *m + j]={0,0};}
}
}

// print input variables
//printf("In cuda: %f\t%f\n",creal(*a[0][0]),cimag(*a[0][0]));
//printf("In cuda: %f\t%f\n",creal(*tau00[0][0]),cimag(*tau00[0][0]));
//printf("In cuda: %d\t%d\t%d\t%d\n",*m,*n,*lda,*kkr_size);

//for (int i=0;i< *m * *n;i++)
//{
//	printf("In cuda: %f\t%f\n",creal(a[i]),cimag(a[i]));
//}

//printf("initialized or not? %d\n",initialized);
if (!initialized)
{
printf("CUDA memory assigned \n");
error=cudaMalloc((void**)&aDev,  sizeof(cuDoubleComplex)* *m * *m);
if (error != cudaSuccess) fprintf(stderr,"\nError1: %s\n",cudaGetErrorString(error));
error=cudaMalloc((void**)&pivotArray,  sizeof(int) * *m);
if (error != cudaSuccess) fprintf(stderr,"\nError2: %s\n",cudaGetErrorString(error));
error=cudaMalloc((void**)&infoArray,  sizeof(int));
if (error != cudaSuccess) fprintf(stderr,"\nError3: %s\n",cudaGetErrorString(error));
error=cudaMalloc(&tau00Dev, sizeof(cuDoubleComplex)* *m * *m);
if (error != cudaSuccess) fprintf(stderr,"\nError4: %s\n",cudaGetErrorString(error));
initialized = true;
}

cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
cusolverDnHandle_t cusolverHandle;
cusolverStatus_t cusolverStatus;
cusolverDnCreate(&cusolverHandle);
cusolverDnZgetrf_bufferSize(cusolverHandle, *m, *m, aDev, *m, &Lwork);
//printf("Lwork is %d\n",Lwork);

error=cudaMalloc((void**)&workArray, Lwork*sizeof(cuDoubleComplex));
if (error != cudaSuccess) fprintf(stderr,"\nError5: %s\n",cudaGetErrorString(error));

cudaEventRecord(start); 
error = cudaMemcpy(aDev, a, sizeof(cuDoubleComplex)* *m * *m, cudaMemcpyHostToDevice);
if (error != cudaSuccess) fprintf(stderr,"\nError6: %s\n",cudaGetErrorString(error));
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_copyin, start, stop);

error = cudaMemcpy(tau00Dev, tau00, sizeof(cuDoubleComplex)* *m * *m, cudaMemcpyHostToDevice);
if (error != cudaSuccess) fprintf(stderr,"\nError7: %s\n",cudaGetErrorString(error));



cudaEventRecord(start);
cusolverStatus = cusolverDnZgetrf(cusolverHandle, *m, *m, aDev, *m, workArray, pivotArray, infoArray);
//if (cusolverStatus == CUSOLVER_STATUS_SUCCESS)
//  printf("cuSOLVER ZGETRF SUCCESSFUL! \n");
//else
//  printf("cuSOLVER ZGETRF UNSUCCESSFUL! \n");

cusolverStatus = cusolverDnZgetrs(cusolverHandle,CUBLAS_OP_N,*m,*m,aDev,*m, pivotArray,tau00Dev,*m,infoArray); 
//if (cusolverStatus == CUSOLVER_STATUS_SUCCESS)
//  printf("cuSOLVER ZGETRS SUCCESSFUL! \n");
//else
//  printf("cuSOLVER ZGETRS UNSUCCESSFUL! \n");
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_compute, start, stop);     

//printf("4\n");
//cudaMemcpy(a, tau00Dev, *m * *m *sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);//copy the whole inverse back to a
//printf("In cuda:tau00(10)(5) %f\t%f\n",creal(tau00[4][9]),cimag(tau00[4][9]));
//printf("In cuda:tau00(1)(1) %f\t%f\n",creal(tau00[0][0]),cimag(tau00[0][0]));

cudaEventRecord(start);
for (int i=0;i<*block_size;i++)
{
	for (int j=0;j<*block_size;j++)
	{
		cudaMemcpy(&b[i+*block_size * j], &tau00Dev[i+*m * j], sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
	}
}
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&time_copyout, start, stop);

//Print the time (in ms) for GPU data transfer and GPU compute
//printf("Time for copyin: %f\tfor copyout: %f\tfor compute inverse: %f\n",time_copyin*0.001,time_copyout*0.001,time_compute*0.001);

//clean up
cusolverDnDestroy(cusolverHandle);
cudaFree(workArray);

// LFu: clean up tau00
free(tau00);
}
