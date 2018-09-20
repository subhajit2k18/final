
#ifndef HEADER_CPP_H
#define HEADER_CPP_H
	
#define nint 2	
	
	#include<cassert>
	extern float *material_data;
	extern int *connect_matrix;
	extern float *nodal_data;
	extern double * elemental_data;
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
	extern float *sol;
	extern int n_node;
	extern int nnz;
	extern int *boundary_nodes;
	extern double *u;
void meshing(float d_length,float d_breadth,float d_height,int x_nelement,int y_nelement,int z_nelement);
void color_mesh(int nelement);

void boundary(float end_load,int x_nelement,int y_nelement ,int z_nelement);

void getGausspoints(double *gausspoints,double *weights,int );
void solve(float **k_m,float *f_m,float *sol,int size);

double d_phi1(int i,double zeta,double eta,double zi,char type);
double phi1(int i,double zeta,double eta,double zi);
void element_stiffness1(int elementno ,double *elemental_data, float * material_data);
void Reorder_connectivity(int *element_color,int *color_size,int noofcolors,int *Reorder_connect);
void Reorder_nodal_data(int *Reorder_connect,int *color_size,int noofcolors,float *Reorder_nodal_data);
void Reorder_csr_index(int *element_color,int *csr_index,int *re_csr_index,int noofcolors,int *color_size);
void data_read(int n_node,int nelement);
void element_stiffness2(int elementno,float *d_material_data,float *elemental_data,int node_index,int mxelement,float *);
int viz(double *sol, double *rho, int *connect_matrix,int loop);








#endif
