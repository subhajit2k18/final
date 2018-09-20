#pragma once
#include<stdio.h>

template<typename T>
__device__ void print(T *array,int nrows,int size)
{
	for(int i=1;i<=size;i++)
	{
		printf("%f ",(float)array[i-1]);
		if((i)%nrows==0)
		printf("\n");
	}

}
