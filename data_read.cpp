#include<iostream>
#include<fstream>
#include"header_cpp.h"
using namespace std;
void data_read(int n_node,int nelement)
{

	int size=0;
//	nodal_data = new float[3*n_node];
  //      connect_matrix = new int[8*nelement];
	Nodes_noofelements = new int[n_node/3];
	Nodesharingelements_ptr=new int[n_node/3];
	//Nodesharingelements = new int[];

	ifstream myfile;
        myfile.open("./Mesh/nodal_data.txt");
        if(myfile.is_open())
        {
                for(int i=0;i<3*n_node;i++){
                myfile>>nodal_data[i];
                }
	myfile.close();
        } else
        cerr<<"File read unsuccessful"<<endl;
     
	 myfile.open("./Mesh/connect_matrix.txt");
        if(myfile.is_open())
        {
                for(int i=0;i<8*nelement;i++){
                myfile>>connect_matrix[i];
                }
	myfile.close();
        } else cerr<<"File read unsuccessful"<<endl;

	 myfile.open("./Mesh/Nodesharingelements_ptr.txt");
        if(myfile.is_open())
        {
                for(int i=0;i<n_node/3;i++){
                myfile>>Nodesharingelements_ptr[i];
                }
	myfile.close();
        } else cerr<<"File read unsuccessful"<<endl;	
   
        myfile.open("./Mesh/Nodesharingelements.txt");
        if(myfile.is_open())
        {
		myfile>>size;
		Nodesharingelements = new int[size];
                for(int i=0;i<size;i++){
                myfile>>Nodesharingelements[i];
                }

	myfile.close();	
        } else cerr<<"File read unsuccessful"<<endl;
	myfile.open("./Mesh/element_color.txt");
        if(myfile.is_open())
        {
             	for(int i=0;i<nelement;i++){
                myfile>>element_color[i];
                }
	myfile.close();
        } else cerr<<"File read unsuccessful"<<endl;

	 myfile.open("./Mesh/Nodes_noofelements.txt");
        if(myfile.is_open())
        {
                for(int i=0;i<n_node/3;i++){
                myfile>>Nodes_noofelements[i];
                }
        myfile.close();
        } else cerr<<"File read unsuccessful"<<endl;

}
