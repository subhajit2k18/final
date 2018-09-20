#include"print_fn.h"
#include<stdio.h>
#include"header.h"
//extern __device__ int *color_data;
//template<typename T>
//__device__ void print(T *array,int nrows,int size);
template<typename T>
__device__ void printc(T *array,int nrows,int size)
{
	for(int i=1;i<=size;i++)
	{
		printf("%f ",array[i-1]);
		if((i)%nrows==0)
		printf("\n");
	}

}


__global__ void Assembly(int *element_color_data,int mxelement,float *d_material_data,int connect_index, float *d_element_matrix, double* d_rho, double penal )
{
//	if(threadIdx.x==0) printf("max element %d\n",mxelement);
	//int local_elementno=threadIdx.x;
	int local_elementno=(gridDim.x*blockIdx.y*blockDim.x*blockDim.y)+(blockIdx.x*blockDim.x*blockDim.y)+ blockDim.x*threadIdx.y+threadIdx.x;
	int local_threadid=blockDim.x*threadIdx.y+threadIdx.x;
//	printf("local element %d\n",local_elementno);
	if(local_elementno < mxelement)
	{
		int global_elementno=element_color_data[local_elementno]+1;
		int node_connect[8];	//for 8 noded brick element 
		//int csr_index[64];
//	printf("global no %d\n",global_elementno);	
		for(int i=0;i<8;i++)
		{node_connect[i]=dg_connect_matrix[connect_index*8+i*mxelement+local_elementno];
		//printf("%d\n",node_connect[i]);
		}
//printf("thread %d ,index %d\n",local_elementno,dg_connect_matrix[connect_index*8+0*mxelement+local_elementno]);
		//if(global_elementno==15)
		//printc(node_connect,1,8);
	 float elemental_data[576];
	for(int i=0;i<576;i++) elemental_data[i]=0;
//		for(int i=local_threadid;i<576;i=i+blockDim.x*blockDim.y) elemental_data[i]=d_element_matrix[i];
		//if(local_threadid==0) 
		//for(int i=0;i<576;i++) elemental_data[i]=d_element_matrix[i];
//__syncthreads();
	element_stiffness(local_elementno,d_material_data,elemental_data,node_connect,connect_index,mxelement);
//__syncthreads();
		//if(global_elementno==33767)
		//printc(elemental_data,24,576);
		//#pragma unroll	
		for(int i=0;i<8;i++)
		{
			
			for(int j=0;j<24;j++)
			{
			dg_global_matrix_data[dg_global_matrix_ptr[3*node_connect[i]]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]  += elemental_data[3*i*24+j]*powf(d_rho[global_elementno-1],penal);
			dg_global_matrix_data[dg_global_matrix_ptr[3*node_connect[i]+1]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]+=elemental_data[(3*i+1)*24+j]*powf(d_rho[global_elementno-1],penal); 
			dg_global_matrix_data[dg_global_matrix_ptr[3*node_connect[i]+2]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]+=elemental_data[(3*i+2)*24+j]*powf(d_rho[global_elementno-1],penal); 
			//printf("%d ",3*dg_csr_index_data[(global_elementno-1)*64+i*8+j/3]+j%3);
			dg_global_matrix_col[dg_global_matrix_ptr[3*node_connect[i]]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]=3*node_connect[j/3]+j%3;
			dg_global_matrix_col[dg_global_matrix_ptr[3*node_connect[i]+1]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]=3*node_connect[j/3]+j%3;
			dg_global_matrix_col[dg_global_matrix_ptr[3*node_connect[i]+2]+3*dg_csr_index_data[(connect_index)*64+i*8*mxelement+(j/3)*mxelement+local_elementno]+j%3]=3*node_connect[j/3]+j%3;
	//	__syncthreads();
			}
		//__syncthreads();
//printf("\n");
		}

//__syncthreads();

//printf("hello\n");
//printf("%d\n",element_color_data[threadIdx.x]);
//printf("color data %d\n",color_data[threadIdx.x]);
//printf("%p\n",color_data);
//printf("%p\n",element_color_data);
//if(threadIdx.x==0)
//printf("global element no %d\n",global_elementno);
//printf("value %d\n",2/3);
//if(global_elementno==3)
//printc(dg_global_matrix_col,24,576);

//printf("%f %f\n",elemental_data[0],elemental_data[1]);
	}
//printf("hello\n");
//__syncthreads();






}
