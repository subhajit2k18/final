#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<cstdlib>
#include<fstream>
using namespace std;
//#define SinglePrecision
//void changeUform(double *sol, REAL4 *u, int **connect_matrix, double *rho);
//Enter node number
int NumX = 101;
int NumY = 51;
int NumZ = 51;


#include "DataTypes.h"
#include "visualization.h"
/*struct REAL4
{
	double x; 
	double  y;
	double  z;
	double  w;
};*/

int viz(double *sol, double *rho, int *connect_matrix,int loop)
{
	
	int size = NumX*NumY*NumZ;
	int Nele = (NumX-1)*(NumY-1)*(NumZ-1);
//	sol = (double*)malloc(3*size*sizeof(double));
//	rho = (double*)malloc(Nele*sizeof(double));
//	connect_matrix = (int*) malloc(sizeof(int )* Nele*8);		// Rows
	   	
	struct REAL4 *u =  new REAL4[size];
	
	/*..........................................................................*/
	/*for(int i=0; i< size; i++)
	{	sol[i] = 1.0;	}
	
	for(int i=0; i< Nele; i++)
	{	rho[i] = 1.0;	}*/
	
/*	FILE *f1 = fopen("disp.txt" , "r");          
         if (f1 == NULL)
         {
        	printf("Error Reading File 1\n");
        	exit (0);
    	 }
      		for(int i=0; i<size*3; i++)
      		{
           		
             		fscanf(f1, "%lf", &sol[i]);
             	         		
        	}
       	fclose(f1);
	
	FILE *f2 = fopen("rho.txt" , "r");          
         if (f2 == NULL)
         {
        	printf("Error Reading File 1\n");
        	exit (0);
    	 }
      		for(int i=0; i<Nele; i++)
      		{
           		
             		fscanf(f2, "%lf", &rho[i]);
             	         		
        	}
       	fclose(f2);
*/	
	/*..........................................................................*/
/*	 int a;
	 FILE *f3 = fopen("connect_matrix.txt" , "r");          
         if (f3 == NULL)
         {
        	printf("Error Reading File 1\n");
        	exit (0);
    	 }
      		for(int i=0; i<Nele*8; i++)
      		{
           		
             		fscanf(f3, "%d", &a);
             		connect_matrix[i] = a;
           		
        	}
       	fclose(f3);
*/       	
	/*..........................................................................*/
	
	changeUform(sol, u, connect_matrix, rho,size, Nele);
	char sBuffer[512];
	sprintf(sBuffer, "./output/Data.%04d.bin", loop);
	WriteBinary(sBuffer, NumX, NumY, NumZ, u);
	return 0;
}
