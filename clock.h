#pragma once
#include<sys/time.h>

//gives time in sec
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
if (result != cudaSuccess) {
fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
assert(result == cudaSuccess);
}
#endif
return result;
}

class StopWatch_GPU
{
  cudaEvent_t startTime, stopTime;
  public:
	StopWatch_GPU()  { cudaEventCreate(&startTime); cudaEventCreate(&stopTime); }
	  ~StopWatch_GPU() { cudaEventDestroy(startTime); cudaEventDestroy(stopTime); }
	  inline void start() { cudaEventRecord(startTime,0); }
	  inline double stop() { return elapsed(); }
	  inline double elapsed() 
  	{
    	cudaEventRecord(stopTime,0);
    	cudaEventSynchronize(stopTime);
    	float result;
    	cudaEventElapsedTime(&result, startTime, stopTime);
    	return result/1000.0;    // 1000 mSec per Sec
  	}
};


class StopWatch_CPU
{
 	 timeval startTime, stopTime, diffTime;
 	public:
  	
	StopWatch_CPU() { start(); }
  	inline void start() { gettimeofday(&startTime,NULL); }
  	inline double stop() { return elapsed(); }
  	inline double elapsed() {
    	gettimeofday(&stopTime,NULL);
    	timersub(&stopTime, &startTime, &diffTime);
    	return diffTime.tv_sec + diffTime.tv_usec/1000000.0; // 10^6 uSec per Sec
  	}
};
