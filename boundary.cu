#include<stdio.h>
#include<math.h>
#include"header.h"
__global__ void boundary_condtn(int *d_boundary_nodes,int maxnodes,int n_node,float end_load,int noof_endloaded_nodes,float *d_global_vector)
{
        int threadid=blockIdx.x*blockDim.x+threadIdx.x;
        if(threadid<maxnodes)
        {
        int k1=dg_global_matrix_ptr[d_boundary_nodes[threadid]];
        int k2=dg_global_matrix_ptr[d_boundary_nodes[threadid]+1];

                for(int i=0;i<(k2-k1);i++)
                {
                if(dg_global_matrix_col[k1+i]==d_boundary_nodes[threadid])
                dg_global_matrix_data[k1+i]=1;
                else
                dg_global_matrix_data[k1+i]=0;

                }

	}

        if(threadid < noof_endloaded_nodes)
        {
                end_load=end_load/noof_endloaded_nodes;
                 int k1=d_boundary_nodes[maxnodes+threadid];
                d_global_vector[k1]+=-end_load;
        }
//__syncthreads();
}
