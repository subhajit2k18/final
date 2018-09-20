#include<iostream>
#include"header_cpp.h"
using namespace std;

//connect_matrix stores values elementwise sequence like  n1n2n3n4n5n6n7n8.....n1n2n3n4....

void meshing(float d_length,float d_breadth,float d_height,int x_nelement,int y_nelement, int z_nelement)
{
int y,x,z,p;
	Nodes_noofelements=new int[(x_nelement+1)*(y_nelement+1)*(z_nelement+1)];
	for(int i=0;i<n_node/3;i++)
	Nodes_noofelements[i]=0;
//	Nodesharingelements=new int[];
	nodal_data =new float[3*(x_nelement+1)*(y_nelement+1)*(z_nelement+1)];
	connect_matrix=new int [x_nelement*y_nelement*z_nelement*8];
	for(int i=1;i<=x_nelement*y_nelement*z_nelement;i++)
	{
	p=(i-1)/(x_nelement*y_nelement);
	y=i/x_nelement;
	x=i%x_nelement;
		if(x==0)
			{y=y*(x_nelement+1)-2;
			connect_matrix[(i-1)*8]=y+p*(x_nelement+1);
			connect_matrix[(i-1)*8+1]=y+1+p*(x_nelement+1);
			connect_matrix[(i-1)*8+2]=y+x_nelement+2+p*(x_nelement+1);
			connect_matrix[(i-1)*8+3]=y+x_nelement+1+p*(x_nelement+1);
			connect_matrix[(i-1)*8+4]=connect_matrix[(i-1)*8]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+5]=connect_matrix[(i-1)*8+1]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+6]=connect_matrix[(i-1)*8+2]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+7]=connect_matrix[(i-1)*8+3]+(x_nelement+1)*(y_nelement+1);
			 }
		else
		{
			y=y*(x_nelement+1)+x-1;
			connect_matrix[(i-1)*8]=y+p*(x_nelement+1);
			connect_matrix[(i-1)*8+1]=y+1+p*(x_nelement+1);
			connect_matrix[(i-1)*8+2]=y+x_nelement+2+p*(x_nelement+1);
			connect_matrix[(i-1)*8+3]=y+x_nelement+1+p*(x_nelement+1);			
			connect_matrix[(i-1)*8+4]=connect_matrix[(i-1)*8]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+5]=connect_matrix[(i-1)*8+1]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+6]=connect_matrix[(i-1)*8+2]+(x_nelement+1)*(y_nelement+1);
			connect_matrix[(i-1)*8+7]=connect_matrix[(i-1)*8+3]+(x_nelement+1)*(y_nelement+1);
		}
	for(int j=0;j<8;j++)
	Nodes_noofelements[connect_matrix[(i-1)*8+j]]++;
			
	
	}	
//**************************Find the element no sharing a particular node
{
	int k=0;
	for(int i=0;i<n_node/3;i++)
	k+=Nodes_noofelements[i];
	
	Nodesharingelements=new int [k];
	Nodesharingelements_ptr=new int[n_node/3];
	Nodesharingelements_ptr[0]=0;
	for(int i=1;i<n_node/3;i++)	
	{
	Nodesharingelements_ptr[i]=Nodesharingelements_ptr[i-1]+Nodes_noofelements[i-1];
	}

	for(int i=0;i<(x_nelement*y_nelement*z_nelement);i++)
	{
		for(int j=0;j<8;j++)
		{
		Nodesharingelements[Nodesharingelements_ptr[connect_matrix[i*8+j]]]=i;
		Nodesharingelements_ptr[connect_matrix[i*8+j]]+=1;	
			
		}






	}

	for(int i=0;i<n_node/3;i++)
		Nodesharingelements_ptr[i]=Nodesharingelements_ptr[i]-Nodes_noofelements[i];


//*************************printing the matrix
/*	for(int i=0;i<n_node/3;i++)
	{
		//node_ptr[i]=node_ptr[i]-Nodes_noofelements[i];
		for(int j=0;j<Nodes_noofelements[i];j++)
		cout<<Nodesharingelements[Nodesharingelements_ptr[i]+j]<<"\t";
		cout<<"\n";
	}*/
	//for(int m=0;m<(n_node/3);m++)
	//std::cout<<node_ptr[m]<<endl;
	
	
	

}

	
std::cout<<"conenct matrix generated"<<"\n";
float e_length=d_length/x_nelement;
float e_breadth=d_breadth/y_nelement;
float e_height=d_height/z_nelement;
/////////////////////////////////////////////////////stores nodal data in (x,y,z),(x,y,z).....nodewise sequence
//node numbering starts from 0.	
	for(int i=0;i<(x_nelement+1)*(y_nelement+1)*(z_nelement+1);i++)
	{
	p=i/((x_nelement+1)*(y_nelement+1));
	y=i/(x_nelement+1);
	y=y-p*(y_nelement+1);
	x=i%(x_nelement+1);
	z=i/((x_nelement+1)*(y_nelement+1));
	nodal_data[i*3]=e_length*x;
	nodal_data[i*3+1]=e_breadth*y;
	nodal_data[i*3+2]=e_height*z;
	}
}

void Cal_csrindex(int elementno,int *csr_index)
{
 	//cout<<"inside"<<endl;   
    int element_location=0;;
        //for(int i=0;i<n_node/3;i++)
        //Nodesharing_size+=Nodes_noofelements[i];









        int elementno_node[8];        //For quadrilateral mesh elements having 8 nodes
        for(int i=0;i<8;i++)
        elementno_node[i]=connect_matrix[(elementno-1)*8+i];
//	cout<<elementno_node[0]<<endl;
        for(int i=0;i<8;i++)            //Loop over all the node of elementno
        {
                int k=Nodes_noofelements[elementno_node[i]];
                int *a=new int[k*8];
                int *b=new int[k*8];
                int *c=new int[k*8];
                for(int j=0;j<k;j++)
                {
                        for(int l=0;l<8;l++)
                        a[j*8+l]=connect_matrix[Nodesharingelements[Nodesharingelements_ptr[elementno_node[i]]+j]*8+l];
                        if(Nodesharingelements[Nodesharingelements_ptr[elementno_node[i]]+j]==elementno-1)
                        element_location=j;
                }
                for(int l=0;l<k*8;l++)
                {b[l]=l;c[l]=l;}

//cout<<"Before insertion"<<endl;       
//      for(int m=0;m<k*8;m++)
//      cout<<a[m]<<"\t";
//      cout<<endl;     

////////////////////////////insertion sort algorithm for increasing order
        int temp=0;
        for(int m=1;m<(k*8);m++)
	{
		for(int j=m;j>=1;j--)
                {

                if (a[j] < a[j-1])
                    {

                        temp = a[j];

                        a[j] = a[j-1];

                        a[j-1] = temp;
                        temp=b[j];
                        b[j]=b[j-1];
                        b[j-1]=temp;
                        c[b[j]]=j;
                        c[b[j-1]]=j-1;
                    }
                          else
                        break;


                }
        }
///////////////////////////////*******************************Removing the duplicate entries and correcting the C array
                int p=0,q=0;
                for(int m=1;m<(k*8);m++)
                {
                        if(a[m]!=a[p])
                        {
                        p++;
                        a[p]=a[m];
                        c[b[m]]-=q;
                        }
                        else
                        {       q++;
                        c[b[m]]=c[b[m-1]];
                        }
                }
//**********************************************

//cout<<"till assignment"<<endl;
for(int m=0;m<(8);m++)
        csr_index[i*8+m]=c[element_location*8+m];

//      for(int m=0;m<k*8;m++)
//      cout<<a[m]<<"\t"<<b[m]<<"\t"<<c[m]<<"\n";
        delete []b;
        delete []a;
        delete []c;
        }


//delete []Nodesharingelements_ptr;
}

