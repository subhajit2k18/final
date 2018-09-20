#include<iostream>
#include<sstream> //for stringstream 
#include<cstdlib> //for exit
#include"header.h"
#include<fstream>
#include"header_cpp.h"
#include"clock.h"
#include<cusp/print.h>
#include<cusp/csr_matrix.h>
#include<cusp/krylov/cg.h>
#include<cusp/copy.h>
#include<cmath>
#include<cusp/monitor.h>
using namespace std;
//#define nint 2
	//	int *boundary_nodes;
__device__ __constant__ int *dg_connect_matrix;
__device__ __constant__ float *dg_nodal_data;
__device__ __constant__ float *dg_global_matrix_data;
__device__ __constant__ int *dg_global_matrix_col;
__device__ __constant__ int *dg_global_matrix_ptr;
__device__ __constant__ int *dg_csr_index_data;
void print(int,int,int);
template<typename T>void print(T *,int size,int rowlength);
void print(float **,int);
void Calc_global_rowptr();
void Cal_csrindex(int elementno,int *csr_index);
//float GetGPUmeminfo();
void FE(int x_nelement, int y_nelement, int z_nelement, int loop, double * rho, double penal)
{
	//cudaSetDevice(0);
	StopWatch_GPU clock;
	StopWatch_CPU watch;
	
	int *boundary_nodes;
	int *re_connect_matrix;
	float *re_nodal_data;
	

	float end_load;
	float y_modulus;
	float p_ratio;

	end_load=51;
	y_modulus=1;
	//p_ratio=(float)1/3;	
	p_ratio=0.3;


	//nnz=(x_nelement*3+1)*(y_nelement*3+1)*18*2+(x_nelement*3+1)*(y_nelement*3+1)*27*(z_nelement-1);//No of nonzeros in global assembled matrix
	//cout<<endl<<nnz<<endl;
	int nelement=x_nelement*y_nelement*z_nelement;
	n_node=3*(x_nelement+1)*(y_nelement+1)*(z_nelement+1);//total no of node for 3 DOF per node
        nnz=(x_nelement*3+1)*(y_nelement*3+1)*18*2+(x_nelement*3+1)*(y_nelement*3+1)*27*(z_nelement-1);//No of nonzeros in global assembled matrix for structured mess
	cout<<"Elements "<<x_nelement<<"\t"<<y_nelement<<"\t"<<z_nelement<<endl;
////////////////////////**************************************material data initialising
//cout<<"taking isotropic material"<<"\n";

	material_data=new float[36];
	for(int i=0;i<36;i++)
	material_data[i]=0;
	{
	float lambda=p_ratio*y_modulus/((1+p_ratio)*(1-2*p_ratio));
	float shear_m=y_modulus/(2*(1+p_ratio));
	material_data[0]=material_data[7]=material_data[14]=lambda+2*shear_m;
	material_data[1]=material_data[2]=material_data[6]=material_data[8]=material_data[12]=material_data[13]=lambda;
	material_data[21]=material_data[28]=material_data[35]=shear_m;
	}
	//material_data[0]=material_data[1]=material_data[3]=material_data[4]=material_data[8]=1;
	//print(material_data,36,6);
/////************************************************************************
	int *csr_index_data=	new int[64*nelement];
	int *re_csr_index_data= new int[64*nelement];

	re_connect_matrix=	new int[x_nelement*y_nelement*z_nelement*8];
	boundary_nodes=		new int[(z_nelement+1)*(y_nelement+1)*3+y_nelement+1];//keeps fixed nodes of cantilever as well as end nodes for loading
	re_nodal_data=		new float[24*nelement];

	for(int i=0;i<576;i++)  elemental_data[i]=0;
//********************************************************************************mesh generation
	cout<<"after initialisation"<<endl;
	if(loop==1)
	{	
	data_read(n_node,nelement);
	element_stiffness1(1,elemental_data,material_data);
	}
	//print(elemental_data,576,24);

	cout<<"In middle"<<endl;		
	//meshing(d_length,d_breadth,d_height, x_nelement,y_nelement,z_nelement);
	//color_mesh(x_nelement*y_nelement*z_nelement);
	//print(x_nelement,y_nelement,z_nelement);
	int *element_color_data;//keeps the colourwise element no list
	
		int colors_size[8]={0,0,0,0,0,0,0,0};//Since at max 8 colors are used for structured mesh
	
		element_color_data=new int[x_nelement*y_nelement*z_nelement];	
		for(int i=0;i<x_nelement*y_nelement*z_nelement;i++)
		colors_size[element_color[i]]++;
		int a[8]={0};
		int k=0;
		for(int i=0;i<8;i++) a[i]=colors_size[i];
		for(int i=0;i<x_nelement*y_nelement*z_nelement;i++)
		{
			k=element_color[i];
			element_color_data[k*colors_size[k]+colors_size[k]-a[k]]=i;
			a[k]-=1;	

		}	




	for(int i=0;i<8;i++)
	cout<<colors_size[i]<<endl;

	

//for(int i=0;i<x_nelement*y_nelement*z_nelement;i++)
//cout<<element_color_data[i]<<endl;
Calc_global_rowptr();//calculate global_matrix_ptr 
Calc_boundary_nodes(boundary_nodes,x_nelement,y_nelement,z_nelement);
	for(int i=1;i<=nelement;i++)
	{
		//cout<<"yes"<<endl;
		Cal_csrindex(i,csr_index_data+(64*(i-1)));		//calculate csr_matrix data index for global assembly
		//csr_index_data+=64;	
	}
//cout<<csr_index_data[0]<<endl;
//print(csr_index_data,64*nelement,8);

watch.start();
Reorder_connectivity(element_color_data,colors_size,8,re_connect_matrix);

Reorder_nodal_data(re_connect_matrix,colors_size,8,re_nodal_data);
Reorder_csr_index(element_color_data,csr_index_data,re_csr_index_data,8,colors_size);
cout<<"reordered connectivity time :"<<watch.elapsed()<<endl;

//for(int i=0;i<576;i++) elemental_data[i]=0;
//for(int i=0; i<24*nelement;i++) nodal_data[i]=re_nodal_data[i];
//element_stiffness2(0,material_data,elemental_data,0,32768,re_nodal_data);
//print(elemental_data,576,24);




//print(csr_index_data,4*64,8);
//cout<<"Reordered index"<<endl;
//print(re_csr_index_data,4*64,8);
//print(re_nodal_data,96,1);
//print(re_connect_matrix,128,1);
//print(connect_matrix,80,1);
global_matrix_col=new int[nnz];
global_matrix_data=new float[nnz];
//print(boundary_nodes,(z_nelement+1)*(y_nelement+1)*3+y_nelement+1,1);
//cout<<boundary_nodes[0]<<endl;
//*********************************************************************************Device variables
	int *d_element_color_data;
	float *d_nodal_data;
	int *d_connect_matrix; 
	float *d_material_data;
	float *d_global_vector;
	float *d_global_matrix_data;	//Array for CSR matrix format
	float *d_sol;
	int *d_global_matrix_col;
	int *d_global_matrix_ptr;
	int *d_csr_index_data;
	int *d_boundary_nodes;
	float *d_element_matrix;
	double *d_rho;
	//int *d_re_connect_matrix;
		cudaMalloc(&d_rho,              sizeof(double)*nelement);
		cudaMalloc(&d_element_matrix,   sizeof(float)*576);
		cudaMalloc(&d_sol,		sizeof(float)*n_node);
		cudaMalloc(&d_element_color_data,sizeof(int)*nelement);
		cudaMalloc(&d_nodal_data,	sizeof(float)*24*nelement);
		cudaMalloc(&d_connect_matrix,	sizeof(int)*nelement*8);
		cudaMalloc(&d_material_data,	sizeof(float)*36);
		cudaMalloc(&d_global_matrix_data,sizeof(float)*nnz);		
		cudaMalloc(&d_global_matrix_col,sizeof(int)*nnz);
		cudaMalloc(&d_global_matrix_ptr,sizeof(int)*(n_node+1));
		cudaMalloc(&d_csr_index_data,	sizeof(int)*64*nelement);
		cudaMalloc(&d_global_vector,	sizeof(float)*n_node);
		cudaMalloc(&d_boundary_nodes,	sizeof(int )*((z_nelement+1)*(y_nelement+1)*3+y_nelement+1));
	//	cudaMalloc(&d_re_connect_matrix,sizeof(int)*nelement*8);	
			cudaMemset(d_global_matrix_data,0,sizeof(float)*nnz);
			cudaMemset(d_global_matrix_col,0,sizeof(int)*nnz);
			cudaMemset(d_global_vector,0,	sizeof(float)*n_node);
			cudaMemset(d_sol,0,sizeof(float)*n_node);
clock.start();
		//cudaMemcpy(d_re_connect_matrix,re_connect_matrix,sizeof(int)*nelement*8,cudaMemcpyHostToDevice);
		cudaMemcpy(d_rho, rho,  sizeof(double)*nelement,cudaMemcpyHostToDevice);
		cudaMemcpy(d_element_matrix, elemental_data, sizeof(float)*576, cudaMemcpyHostToDevice);
		cudaMemcpy(d_element_color_data,element_color_data,sizeof(int)*nelement,cudaMemcpyHostToDevice);
		cudaMemcpy(d_nodal_data,re_nodal_data,sizeof(float)*nelement*24,cudaMemcpyHostToDevice);
		cudaMemcpy(d_connect_matrix,re_connect_matrix,sizeof(int)*nelement*8,cudaMemcpyHostToDevice);
		cudaMemcpy(d_material_data,material_data,sizeof(float)*36,cudaMemcpyHostToDevice);
		cudaMemcpy(d_csr_index_data,re_csr_index_data,sizeof(int)*64*nelement,cudaMemcpyHostToDevice);
		cudaMemcpy(d_global_matrix_ptr,global_rowptr,sizeof(int)*(n_node+1),cudaMemcpyHostToDevice);
		cudaMemcpy(d_boundary_nodes,boundary_nodes,sizeof(int)*((z_nelement+1)*(y_nelement+1)*3+y_nelement+1),cudaMemcpyHostToDevice);
				cudaError_t result=cudaMemcpyToSymbol(dg_nodal_data,&d_nodal_data,sizeof(float*));
			   cudaMemcpyToSymbol(dg_connect_matrix,&d_connect_matrix,sizeof(int *));
			   cudaMemcpyToSymbol(dg_global_matrix_data,&d_global_matrix_data,sizeof(float *));
			   cudaMemcpyToSymbol(dg_global_matrix_col,&d_global_matrix_col,sizeof(int *));
			   cudaMemcpyToSymbol(dg_csr_index_data,&d_csr_index_data,sizeof(int *));
			   cudaMemcpyToSymbol(dg_global_matrix_ptr,&d_global_matrix_ptr,sizeof(int *));
cout<<"Elapsed time in cudaMemcpy "<<clock.elapsed()<<endl;	
//assert(result==cudaSuccess);
//cout<<d_element_color_data[0]<<endl;
//dim3 threadperblock(2,3);
/////////////////////////////////////////////******************Kernel launch for Assembly and boundary condition.
//dim3 blocks(1,1);
 k=0;
int noofblocks=ceil((float)colors_size[0]/64);
cout<<"No of blocks :"<<noofblocks<<endl;
dim3 blocks(noofblocks);
dim3 threadperblock(64);


//clock.start();
for(int i=0;i<8;i++)
{

Assembly<<<blocks,threadperblock>>>(d_element_color_data+k,colors_size[i],d_material_data,k,d_element_matrix,d_rho,penal);
k+=colors_size[i];
//cout<<k<<endl;
//cudaThreadSynchronize();
//cudaError_t error=cudaDeviceSynchronize();
//if(error!=cudaSuccess)
//cerr<<"error::"<<cudaGetErrorString(error)<<endl;
}

//cout<<"elapsed time in kernel :"<<clock.elapsed()<<endl;

noofblocks=((z_nelement+1)*(y_nelement+1)*3)/257+1;
threadperblock.x=256;
cout<<"No of boundary nodes "<<(z_nelement+1)*(y_nelement+1)*3<<endl;
cout<<"No of blocks boundary "<<noofblocks<<endl;
blocks.x=noofblocks;
boundary_condtn<<<noofblocks,threadperblock>>>(d_boundary_nodes,(z_nelement+1)*(y_nelement+1)*3,n_node,end_load,y_nelement+1,d_global_vector);

//cudaThreadSynchronize();
//cudaDeviceSynchronize();
//cout<<"GPU memory usage :"<<GetGPUmeminfo()<<endl;
			//cudaMemcpy(global_matrix_data,d_global_matrix_data,sizeof(float)*nnz,cudaMemcpyDeviceToHost);
//			cudaMemcpy(global_matrix_col,d_global_matrix_col,sizeof(int)*nnz,cudaMemcpyDeviceToHost);
//***************************************************************************************************************

		cudaFree(d_element_color_data);
		cudaFree(d_nodal_data);
		cudaFree(d_connect_matrix);
		cudaFree(d_material_data);
		//cudaFree(d_global_matrix_data);
		//cudaFree(d_global_matrix_col);
		//cudaFree(d_global_matrix_ptr);
		cudaFree(d_csr_index_data);
		cudaFree(d_boundary_nodes);
		//cudaFree(d_global_vector);
//print(global_matrix_data,576,24);
//element_stiffness(4);
//print(elemental_data,576,24);

	

/*	for(int i=0;i<n_node/2;i++)
		{for(int j=0;j<(global_rowptr[i+1]-global_rowptr[i]);j++)
		cout<<global_matrix_data[global_rowptr[i]+j]<<"\t";
		cout<<endl;
		}
*/	

//print(sol,n_node,1);
//******************************************************************************************CUSP
	//thrust::device_ptr<float> wrapped_b_vect(d_global_vector);
	typedef typename cusp::array1d_view<thrust::device_ptr<int> > DeviceIndexArrayView;
    	typedef typename cusp::array1d_view<thrust::device_ptr<float> > DeviceValueArrayView;
	
		thrust::device_ptr<float> wrapped_b_vect(d_global_vector);	
		thrust::device_ptr<float> wrapped_d_global_matrix_data(d_global_matrix_data);	
		thrust::device_ptr<int>   wrapped_d_global_matrix_col(d_global_matrix_col);
		thrust::device_ptr<int>	  wrapped_d_global_matrix_ptr(d_global_matrix_ptr);
		thrust::device_ptr<float> wrapped_d_sol(d_sol);
	DeviceIndexArrayView row_offsets (wrapped_d_global_matrix_ptr,wrapped_d_global_matrix_ptr+n_node+1);
	DeviceIndexArrayView col_indices (wrapped_d_global_matrix_col,wrapped_d_global_matrix_col+nnz);
	DeviceValueArrayView values (wrapped_d_global_matrix_data,wrapped_d_global_matrix_data+nnz);
	DeviceValueArrayView sol_d  (wrapped_d_sol,wrapped_d_sol+n_node);
	//HostValueArrayView x (sol,sol+n_node);
	//DeviceValueArrayView b_vect (d_global_vector,d_global_vector+n_node);
	DeviceValueArrayView b_vect(wrapped_b_vect,wrapped_b_vect+n_node);
	typedef cusp::csr_matrix_view<DeviceIndexArrayView,DeviceIndexArrayView,DeviceValueArrayView> DeviceView;
	DeviceView A(n_node,n_node,nnz,row_offsets,col_indices,values);	

	//cusp::array1d<float,cusp::device_memory>x(n_node,0);
//cudaDeviceSynchronize();	
	//cusp::print(b_vect);
	//cusp::print(A.values);
//	cusp::default_monitor<float> monitor(b_vect, 10000, 1e-5);
	cusp::monitor<double> monitor(b_vect,2000,1e-5);
	cusp::krylov::cg(A,sol_d,b_vect,monitor);
//	print(sol_d);
	//cout<<x(n_node);
	float *sol=new float[n_node];
	cudaMemcpy(sol, d_sol,sizeof(float)*n_node,cudaMemcpyDeviceToHost);
	for(int i=0;i<n_node;i++) u[i]=sol[i];
	delete [] sol;
//cudaDeviceSynchronize();
//	cusp::print(b_vect);
//	print(sol_d);



	cudaFree(d_global_vector);
	cudaFree(d_global_matrix_data);
	cudaFree(d_global_matrix_col);
	cudaFree(d_global_matrix_ptr);
	cudaFree(d_sol);


//*********************************************************************************************************************
//cout<<"last node "<<2*(h_nelement*v_nelement+1)+1<< "y displacement "<<sol[2*(h_nelement*v_nelement+1)+1]<<"\n";
/*	ostringstream outfile("");
	outfile<<"file_"<<d_length<<".dat";	
	string filename=outfile.str();
	cout<<filename<<endl;
	ofstream outf(filename.c_str(),ios::app);
	//outf<<h_nelement*v_nelement<<"\t"<<sol[2*(h_nelement*v_nelement+1)+1]<<"\n";
	int testnode=h_nelement*((v_nelement/2)+1)+(v_nelement/2);
	outf<<h_nelement*v_nelement<<"\t"<<sol[2*testnode+1]<<"\n";
*/	//cout<<testnode<<endl;
//delete []re_nodal_data;
//delete []re_connect_matrix;
//delete []boundary_nodes;
//delete []csr_index_data;
//delete []re_csr_index_data;
//delete []global_matrix_data;
//delete []global_matrix_col;
//delete []global_rowptr;
//delete []global_vector;
//delete []material_data;
//delete []elemental_data;
//delete []elemental_vector;
//delete []element_color;
//delete []element_color_data;
//delete []Nodes_noofelements;
//delete []Nodesharingelements;
//delete []Nodesharingelements_ptr;
//delete []nodal_data;
//delete []connect_matrix;
//delete []sol;
//return 0;
}
template<typename T>
void print(T *data,int size,int rowlength)
{
for(int i=0;i<size/rowlength;i++)
{
	for(int j=0;j<rowlength;j++)
		cout<<data[i*rowlength+j]<<"\t";
	cout<<"\n";
}
}

void print(float** data,int size)
{
for(int i=0;i<size;i++)
{
        for(int j=0;j<size;j++)
                cout<<data[i][j]<<"\t";
        cout<<"\n";
}


}



void print(int h_nelement,int v_nelement,int z_nelement)
{
	for(int i=0;i<h_nelement*v_nelement*z_nelement;i=i+1)
	{
	cout<<"values of i"<<i<<endl;
	for(int j=0;j<8;j++)
	cout<<connect_matrix[i*8+j];
	cout<<"\n";
	}
	for(int i=0;i<3*(h_nelement+1)*(v_nelement+1)*(z_nelement+1);i=i+3)
	cout<<nodal_data[i]<<"\t"<<nodal_data[i+1]<<"\t"<<nodal_data[i+2]<<"\n";

}

void Calc_global_rowptr()
{
 	global_rowptr=new int [n_node+1];
		global_rowptr[0]=0;
                global_rowptr[1]=24;
                global_rowptr[2]=48;
                global_rowptr[n_node]=nnz;
                for(int j=1;j<=n_node/3;j++)
                {
                int k=0;
                if(Nodes_noofelements[j-1]==1)
                {k=24;}
                else if(Nodes_noofelements[j-1]==2)
                k=36;
                else if(Nodes_noofelements[j-1]==4)
                k=54;
                else if(Nodes_noofelements[j-1]==8)
                k=81;
                //global_rowptr[j]= Nodes_noofelement[j-1]*
                global_rowptr[3*j-2]=global_rowptr[3*j-3]+k;
                global_rowptr[3*j-1]=global_rowptr[3*j-2]+k;
                global_rowptr[3*j]=global_rowptr[3*j-1]+k;
                //global_rowptr[3*j+1]=global_rowptr[3*j]+k;
                //global_rowptr[3*j+2]=global_rowptr[3*j+1]+k;
                }




}

float GetGPUmeminfo()
{

	size_t free_byte ;
	size_t total_byte ;

	cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte );

	if ( cudaSuccess != cuda_status ){
	printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
	exit(1);
	}

float free_db = (float)free_byte ;
float total_db = (float)total_byte ;
float used_db = total_db - free_db ;

//	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
	return used_db/1024/1024;
}



