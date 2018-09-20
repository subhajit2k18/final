
#ifndef HEADER_H
#define HEADER_H
#define nint 2

	#include<cassert>

	extern float *material_data;
	extern int *connect_matrix;
	extern float *nodal_data;
	//extern float * elemental_data;
	extern float * elemental_vector;
	//extern float **global_matrix;
	extern float *global_vector;
	extern int *Nodes_noofelements;//keeps the count of  elements containing a particular node....data arranged  nodewise
	extern int *Nodesharingelements;//keeps the list of elements number sharing a node....data arranged nodewise
	extern int *Nodesharingelements_ptr;// points to first location in Nodesharingelements for each node
	
	extern float *global_matrix_data;
	extern int *global_matrix_col;
	extern int *global_rowptr;
	extern int *element_color;//stores color of elements...data stored elementwise
	
	extern int n_node;
	extern int nnz;
	extern int *boundary_nodes;	
		
	 extern __device__ __constant__ int *dg_connect_matrix;
	 extern __device__ __constant__ float *dg_nodal_data;
	 extern __device__ __constant__ float *dg_global_matrix_data;
	 extern __device__ __constant__ int *dg_global_matrix_col;
	 extern __device__ __constant__ int *dg_global_matrix_ptr;
	 extern __device__ __constant__ int *dg_csr_index_data;


void meshing(float d_length,float d_breadth,float d_height,int x_nelement,int y_nelement,int z_nelement);
void color_mesh(int nelement);

void Calc_boundary_nodes(int *,int x_nelement,int y_nelement ,int z_nelement);

__device__ void getGausspoints(double *gausspoints,double *weights,int );
void solve(float **k_m,float *f_m,float *sol,int size);

__device__ double d_phi(int i,double zeta,double eta,double zi,char type);
__device__ double phi(int i,double zeta,double eta,double zi);
__device__ void element_stiffness(int elementno,float *d_material_data,float *elemental_data,int *node_connect,int node_index,int mxelement);

__global__ void Assembly(int *,int,float *,int connect_index,float *, double *, double);
void sort(float *data,int *col,int *ptr);

__global__ void boundary_condtn(int *boundary_nodes,int maxnodes,int n_node,float end_load,int noof_endloaded_nodes,float *d_global_vector);

#endif
