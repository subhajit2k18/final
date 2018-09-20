
#include<iostream>
#include"header_cpp.h"
using namespace std;
//boundary_nodes=new int[(z_nelement+1)*(y_nelement+1)*3];

void Calc_boundary_nodes(int *boundary_nodes,int x_nelement,int y_nelement,int z_nelement)
{
//**********************************boundary condition for fixed end of cantilever beam for structured mess
//	boundary_nodes=new int[(z_nelement+1)*(y_nelement+1)*3];
        //find the node number and store in array
        int size=(z_nelement+1)*(y_nelement+1)*3;
	for(int i=0;i<y_nelement+1;i++)
        {
        int a=i*x_nelement+i;
        boundary_nodes[3*(z_nelement+1)*i]=3*a;
        boundary_nodes[3*(z_nelement+1)*i+1]=3*a+1;
        boundary_nodes[3*(z_nelement+1)*i+2]=3*a+2;
                for(int j=0;j<z_nelement;j++)
                {
                boundary_nodes[3*(z_nelement+1)*i+3*j+3]=3*a+(j+1)*3*(x_nelement+1)*(y_nelement+1);
                boundary_nodes[3*(z_nelement+1)*i+3*j+4]=boundary_nodes[3*(z_nelement+1)*i+3*j+3]+1;
                boundary_nodes[3*(z_nelement+1)*i+3*j+5]=boundary_nodes[3*(z_nelement+1)*i+3*j+3]+2;
                }
        }
	
	for(int i=0;i<y_nelement+1;i++)
	{
	 int a=i*(x_nelement+1)+x_nelement;//node number in lowest layer
	//boundary_nodes[size+i]=3*(a+(x_nelement+1)*(y_nelement+1)*z_nelement)+2;     // Force on top layer
	boundary_nodes[size+i]=3*(a)+2;		// Force on bottom layer
	}


}
/*
	for(int i=0;i<((z_nelement+1)*(y_nelement+1)*3);i++)
	{
//	int k=global_rowptr[nodenumber[i]];
		for(int k=0;k<n_node;k++)	//can be further optimised to check whether the column index exist in particular node nnz
		{
		       
		 int ptr=global_rowptr[k];
		
			if(k!=nodenumber[i])
			{
				for(int j=0;j<global_rowptr[k+1]-global_rowptr[k];j++)
					if(global_matrix_col[ptr+j]==nodenumber[i])
					{global_matrix_data[ptr+j]=0;break;}		
			}
			else if(k==nodenumber[i])
			{
			for(int j=0;j<global_rowptr[k+1]-global_rowptr[k];j++)
				if(global_matrix_col[ptr+j]!=nodenumber[i])
					global_matrix_data[ptr+j]=0;
				else
					global_matrix_data[ptr+j]=1;		
			}
		}
	
		global_vector[nodenumber[i]]=0;
	}

//***********************************************************************boundary condition for end point load
	end_load=end_load/(y_nelement+1);
	for(int i=0;i<y_nelement+1;i++)
        {
        	int a=i*(x_nelement+1)+x_nelement;//node number in lowest layer
        	global_vector[3*(a+(x_nelement+1)*(y_nelement+1)*z_nelement)+2]+=-end_load;
        //cout<<3*(a+(x_nelement+1)*(y_nelement+1)*z_nelement)+1<<endl;
        }





}*/
