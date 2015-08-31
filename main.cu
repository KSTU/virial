#include <stdlib.h>
#include <time.h>
#define N 1024

// declare the kernel
__global__ void daxpy(int n, double a, double *x, double *y){
int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < N){
		y[i] += a*x[i];
	}
}

int main(void){
	double *x, *y, a, *dx, *dy;
	x = (double *)malloc(sizeof(double)*N);
	y = (double *)malloc(sizeof(double)*N);
	// initialize x and y
	srand(time(NULL));

	// allocate device memory for x and y
	cudaMalloc(dx, N*sizeof(double));
	cudaMalloc(dy, N*sizeof(double));
	// copy host memory to device memory
	cudaMemcpy(dx, x, N*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dy, y, N*sizeof(double), cudaMemcpyHostToDevice);
	// launch the kernel function
	a=0.1;
	daxpy<<<N/64,64>>>(N, a, dx, dy);
	// copy device memory to host memory
	cudaMemcpy(y, dy, N*sizeof(double), cudaMemcpyDeviceToHost);
	// deallocate device memory
	cudaMemFree(dx);
	cudaMemFree(dy);
	free(x);
	free(y);
}
