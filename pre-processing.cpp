#include<iostream>
#include"header_cpp.h"
using namespace std;

void Reorder_connectivity(int *element_color,int *color_size,int noofcolors,int *Reorder_connect)
{
	int index=0;
	for(int i=0;i<noofcolors;i++)
	{
		if(i>0) index+=color_size[i-1];
	//	cout<<"index "<<index<<endl;
		for(int j=0;j<color_size[i];j++)
		{
			for(int k=0;k<8;k++)
			{Reorder_connect[(index*8)+color_size[i]*k+j]=connect_matrix[element_color[index+j]*8+k];}
		}

	}

}

void Reorder_nodal_data(int *Reorder_connect,int *color_size,int noofcolors,float *Reorder_nodal_data)
{
	int index=0;int k=0;
	for(int i=0;i<noofcolors;i++)
	{
		if(i>0) index+=color_size[i-1];
		cout<<"index "<<index<<endl;
		for(int j=0;j<color_size[i]*8;j++)
		{
			k=index*24+(j/color_size[i])*color_size[i]*3+(j%color_size[i]);
			Reorder_nodal_data[k]=nodal_data[3*Reorder_connect[index*8+j]];
			Reorder_nodal_data[k+color_size[i]]=nodal_data[3*Reorder_connect[index*8+j]+1];
			Reorder_nodal_data[k+2*color_size[i]]=nodal_data[3*Reorder_connect[index*8+j]+2];
			
			//cout<<index*24+(j/color_size[i])*color_size[i]*3+j%color_size[i]<<endl;
		}
	}



}

void Reorder_csr_index(int *element_color,int *csr_index,int *re_csr_index,int noofcolors,int *color_size)
{
	int index=0;
	for(int i=0;i<noofcolors;i++)
	{
		if(i>0) index+=color_size[i-1];
		for(int j=0;j<color_size[i];j++)
		{
			for(int k=0;k<64;k++)
			{
			
			re_csr_index[index*64+k*color_size[i]+j]=csr_index[64*element_color[index+j]+k];



			}
		

		}
	}


}


