#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"header_cpp.h"
#include"clock.h"
#include<fstream>
#include<iostream>
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
void FE(int x_nelement, int y_nelement, int z_nelement, int loop, double * rho, double penal);
void OCupdate(int nx,int ny,int nz,double *rho,double volfrac,double MinDens,double *d_comp);
void check(int nx,int ny,int nz,double rmin,double *rho, double *d_comp);
      
        float *nodal_data;
        int *connect_matrix;
        int *element_color;
        double *elemental_data;  //stores elemental stiffness matrix
        float *elemental_vector;  //stores elemental source vector
        float *material_data;   //stores material properties
        //float **global_matrix;        //stores global stiffness matrix
        int *Nodes_noofelements;//keeps the no of elements containing a particular node
        int *Nodesharingelements;// keeps the list of element number sharing a node
        int *Nodesharingelements_ptr;
        float *global_matrix_data;
        int *global_matrix_col;
        int *global_rowptr;
        float * global_vector;  //stores global source vector
        //float * sol;
        int n_node;
        int nnz;
        int nelement;
	double *u;
//      int *boundary_nodes;

int main()
{
StopWatch_CPU cpu_clock;
StopWatch_GPU feTime;
StopWatch_GPU OCTime;
StopWatch_GPU CheckTime;
int nx = 200;
int ny = 100;
int nz = 100;
int n_node= (nx+1)*(ny+1)*(nz+1);
int n_elem= nx*ny*nz;
	//elemental_data=		new float[576];////////for linear shape function  on rectangular element having 3 degree of freedom per node..24*24
	//elemental_vector=	new float[24];
			
/*	sol= new float[3*n_node];
	for(int i=0;i<3*n_node;i++)
	{ 
	sol[i]=0;
	}
*/

	elemental_data     =    new double[576];
	nodal_data = new float[3*n_node];
        connect_matrix = new int[8*n_elem];
        element_color = new int[n_elem];
	data_read(n_node, n_elem);


double volfrac = 0.3;
double penal = 3;
double rmin = 3.0;
double MinDens = 0.1;

u = (double *)malloc(3*n_node * sizeof(double));
for(int i = 0; i < 3*n_node; i++)
{
	u[i] = 1.0;
}
double ue[24];
double Ke[24][24];
double uKe[24];
double uKeu = 0.0;

for(int i = 0; i < 24; i++)
{
	ue[i] = 0.0;
	uKe[i] = 0.0;
	for(int j = 0; j < 24; j++)
	{
		Ke[i][j] = 0.0;
		
	}
}
double comp = 0.0; //Compliance
double *d_comp = (double *)malloc(nx*ny*nz*sizeof(double));
for(int i = 0; i < nx*ny*nz; i++)
{
	d_comp[i] = 0.0;
}

int *connectMatrix = (int *)malloc(nx*ny*nz*8*sizeof(int));
for(int i = 0; i < nx*ny*nz*8; i++)
{
	connectMatrix[i] = connect_matrix[i];
}

double *rho =(double *) malloc(nx*ny*nz*sizeof(double));
double *rho_old = (double *)malloc(nx*ny*nz*sizeof(double));

for(int i = 0; i < nx*ny*nz; i++)
{
	rho[i] = volfrac;
	rho_old[i] = volfrac;
}

int loop = 0;
double change = 1;

double checkMax = 0.0;
double checkSum = 0.0;
//Main loop

cpu_clock.start();
while( change > 0.01 && loop <5)
{
	//viz(u,rho,connectMatrix,loop);
	loop++;
	printf("\n Loop no :%d\nCompliance = %lf\n",loop,comp);
	for(int i = 0; i < nx*ny*nz ; i++)
	{
		rho_old[i] = rho[i];
	}
	
	///////////// u = FEA();
	feTime.start();
	FE(nx,ny,nz,loop,rho,penal);
	feTime.stop();
	//for(int i=0;i<3*n_node;i++) printf("%lf \n",u[i]);
	if(loop==1)
	{
	for(int i=0;i<24;i++)
	for(int j=0;j<24;j++){
      	Ke[i][j]=elemental_data[i*24+j];
		}
	}
	///////////// Ke = FEA();
	
	comp = 0.0;

	for(int i = 0; i < nx * ny * nz; i++)
	{
		for(int ii = 0; ii < 8; ii++)
		{
			ue[ii * 3] =  u[connectMatrix[i * 8 + ii] * 3];
			ue[ii * 3 + 1] = u[connectMatrix[i * 8 + ii] * 3 + 1];
			ue[ii * 3 + 2] = u[connectMatrix[i * 8 + ii] * 3 + 2];
			
		}
//		if(i == 167){for(int j = 0; j < 24; j++) printf("\nUe = %lf", ue[j]);}
		for(int j = 0; j <24; j++)
		{
			uKe[j] = 0.0;
		}
		for(int j = 0; j < 24; j++)
		{
			for(int k = 0; k < 24; k++)
			{
				uKe[j] += ue[k] * Ke[k][j]; 
			}
		}
	/*	if(i == 167)
		{
			printf("\n");	
			for(int l = 0; l < 24; l++)
			{
				for(int m = 0; m < 24; m++)
				{
					printf("%lf ",Ke[l][m]);
				}
				printf("\n");
			}
		}
		*/
 //		if(i == 167){for(int j = 0; j < 24; j++) printf("\nUKe = %lf", uKe[j]);}
		uKeu = 0.0;
		for(int j = 0; j < 24; j++)
		{
			uKeu += uKe[j] * ue[j];
		}
//		if(i == 167){for(int j = 0; j < 24; j++) printf("\nUe = %lf", ue[j]);}
//		{printf("\ni = %d UKeu = %lf",i, uKeu);}
		comp = comp + pow(rho[i], penal) * uKeu;
		d_comp[i] = - penal * pow(rho[i], (penal - 1)) * uKeu;
	}
	
	printf("\n Second Compliance = %lf\n",comp);	
	//Sensititivity filter
//	for(int i=0;i<nx * ny * nz;i++) printf("\n d_comp = %lf\n",d_comp[i]);
	/////////////////////////////////////////////////////////////////////////////////////////////////////




	CheckTime.start();
	check(nx,ny,nz,rmin,rho,d_comp);
	CheckTime.stop();
//	for(int i=0;i<nx * ny * nz;i++) printf("\n updated_d_comp = %lf\n",d_comp[i]);


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Design update
	OCTime.start();
	OCupdate(nx,ny,nz,rho,volfrac,MinDens,d_comp);
	OCTime.stop();
//	for(int i=0;i<nx * ny * nz;i++)	printf("\n Rho = %lf\n",rho[i]);	
	
	for(int i = 0; i < nx * ny * nz; i++)
	{
		if(abs(rho[i] - rho_old[i]) > checkMax)
		{
			checkMax = abs(rho[i] - rho_old[i]);
		}
	}
	
	change = checkMax;


	printf("\n				FE :%lf\n                              Check :%lf\n                              OC :%lf",feTime.elapsed(),CheckTime.elapsed(),OCTime.elapsed());		

	
//	change = 0;
}
	printf("Time taken in 3 iteration: %lf \n",cpu_clock.elapsed());

	std::ofstream myfile;
        myfile.open("rho2.txt");
        if(myfile.is_open())
        {
                for(int i=0;i< nx * ny * nz;i++){
                myfile<<rho[i]<<"\n";
                }
        myfile.close();
        }
	
	myfile.open("disp2.txt");
        if(myfile.is_open())
        {
                for(int i=0;i< 3*n_node ;i++){
                myfile<<u[i]<<"\n";
                }
        myfile.close();
        }




return 0;
}

void check(int nx,int ny,int nz,double rmin,double *rho, double *d_comp)
{
   double *New_d_comp = (double*)malloc(sizeof(double)*nx*ny*nz);
  
  	for(int i=0; i<nx*ny*nz; i++)
  	{  New_d_comp[i] = 0.0;     }
  	
  	for(int k=0; k<nz; k++)
  	{
  	   for(int j=0; j<ny; j++)
  	   {
  	      for(int i=0; i<nx; i++)
  	      {
  	      	   double sum = 0.0;
  	       	   // Global no. of this element
  	      	   int id1 = i+(j*nx)+(k*nx*ny);
  	      	   
  	      	   for(int ek = max(k-(int)rmin,0); ek<= min(k+(int)rmin,nz); ek++)
  	      	   {
  	      	      for(int ej = max(j-(int)rmin,0); ej<= min(j+(int)rmin,ny); ej++)
  	      	      {
  	      	          for(int ei = max(i-(int)rmin,0); ei<= min(i+(int)rmin,nx); ei++)
  	      	          {	
  	      	              // Global no. of this element
  	       	              int id2 = ei+(ej*nx)+(ek*nx*ny);
  	      	              
  	      	              double fac = rmin-sqrt((k-ek)*(k-ek)+(j-ej)*(j-ej)+(i-ei)*(i-ei));
              		      sum += max(0.0, fac);
//              		      printf("sum %lf\n",sum); 
              		      New_d_comp[id1] += d_comp[id2]*max(0.0, fac)*rho[id2];
              		      
  	      	          }
  	      	      }
  	      	   }
  	         
  	         New_d_comp[id1] =  New_d_comp[id1]/(rho[id1]*sum);
  	      }
  	   }
  	}
  	
  	for(int i=0; i<nx*ny*nz; i++)
  	{  
  	   d_comp[i] = New_d_comp[i];     
    	}
  printf("Check working");
}

void OCupdate(int nx,int ny,int nz,double *rho,double volfrac,double MinDens,double *d_comp)
{
  double l1 = 0.0;
  double l2 = 1e5;
  int counter = 0;	// Bi-section iteration counter
  double move = 0.2;
  
  double *New_rho = (double*)malloc(sizeof(double)*nx*ny*nz);
  
  	for(int i=0; i<nx*ny*nz; i++)
  	{  New_rho[i] = 0.0;     }
    	double derivative;
	while(l2-l1 > 1e-4 && counter < 1e3)
	{
		   counter++;
		   
	    	  // Bi-section 
	    	  double lmid = 0.5*(l2+l1);
	    	  
	    	  //loop over all elements
		  for(int i=0; i<nx*ny*nz; i++)
	    	  {
	    	       //the gradient must be less than zero everywhere, so
          	       //everything bigger than 0.0 will be truncated
          	       derivative = max(-1.0*d_comp[i], 0.0);
          	       
          	       // Updated d_comp
          	       New_rho[i] = max(MinDens, max(rho[i]-move, min(1.0, min(rho[i]+move, rho[i]*sqrt(derivative/lmid)))));
          	  }
	    	   
	    	   // Summing the gradients
	    	   double sum = 0.0;
	    	   for(int i=0; i<nx*ny*nz; i++)
	    	   {
	    	   	sum += New_rho[i];
	    	   }
	    	 	    	 
	    	 if(sum - volfrac*nx*ny*nz > 0)
      		   {	l1 = lmid;	}
    		 else
      		   {	l2 = lmid;	}
		  
		 // Updating
		 
	    	 
   //free(New_d_comp);
	}
	for(int i=0; i<nx*ny*nz; i++)
                 {
                    rho[i] = New_rho[i];
                 }
	printf("OC working");
}
