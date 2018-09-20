#include<stdio.h>
#include<cassert>
#include<cmath>
#include"header.h"
using namespace std;
//__device__ double jacobian[9];
//__device__ double matrix_b[144];  ///(6x24 matrix contains derivative of shape function in physical coordinates)
///elemental computation uses isoparametric shape function 

__device__ void element_stiffness(int elementno,float *d_material_data,float *elemental_data,int *node_connect,int node_index,int mxelement)
{
	float jacobian[9];
	float matrix_b[144];
	double *gausspoints=new double [nint];
	double *weights=new double [nint];
	getGausspoints(gausspoints,weights,nint);
//printf("gausspoints : %f %f\n",gausspoints[0],gausspoints[1]);	
	float d_shapefn[24];
//	printf("elementno :%d\n",elementno);	
for(int m=0;m<nint;m++)
{

for(int n=0;n<nint;n++)
{
	
for(int p=0;p<nint;p++)
{
//	cout<<"LOOP"<<"\n";
	//cout<<"Points "<<gausspoints[m]<<"\t"<<gausspoints[n]<<"\t"<<gausspoints[p]<<"\n";
	for(int i=0;i<144;i++) {matrix_b[i]=0;}// initialising

	float determinant=0;
///////////////////////////keeping values of derivative of shape function in an array
	for(int i=0;i<24;i++)
{	if(i<8)
	d_shapefn[i]=d_phi(i+1,gausspoints[m],gausspoints[n],gausspoints[p],'z');
	else if(i>=8 && i<16)
	d_shapefn[i]=d_phi(i-7,gausspoints[m],gausspoints[n],gausspoints[p],'e');
	else 
	d_shapefn[i]=d_phi(i-15,gausspoints[m],gausspoints[n],gausspoints[p],'j');
}

//__syncthreads();
/*
if(elementno==1)
{	for(int i=0;i<3;i++)
	{
	for(int j=0;j<8;j++)
	{printf("%f  ",d_shapefn[i*8+j]);}
	printf("hello\n");
	}}*/
//__syncthreads();
////////////////calculating the jacobian value
	jacobian[0]=jacobian[1]=jacobian[2]=jacobian[3]=0;		//initialising 
	jacobian[4]=jacobian[5]=jacobian[6]=jacobian[7]=jacobian[8]=0;
//__syncthreads();	
	for(int j=0;j<8;j++)
		{
		//cout<<"node :"<<nodal_data[3*(connect_matrix[(elementno-1)*8+j])]<<nodal_data[3*(connect_matrix[(elementno-1)*8+j])+2];
		float x,y,z;
		x=dg_nodal_data[node_index*24+j*mxelement*3+elementno];
		y=dg_nodal_data[node_index*24+j*mxelement*3+elementno+mxelement];
		z=dg_nodal_data[node_index*24+j*mxelement*3+elementno+2*mxelement];
		//x=dg_nodal_data[3*node_connect[(j)]];
		//y=dg_nodal_data[3*(node_connect[j])+1];
		//z=dg_nodal_data[3*(node_connect[j])+2];
		jacobian[0]=jacobian[0]+d_shapefn[j]*x;
		jacobian[3]=jacobian[3]+d_shapefn[j+8]*x;
		jacobian[6]=jacobian[6]+d_shapefn[j+16]*x;
		jacobian[1]=jacobian[1]+d_shapefn[j]*y;
		jacobian[4]=jacobian[4]+d_shapefn[j+8]*y;
		jacobian[7]=jacobian[7]+d_shapefn[j+16]*y;
		jacobian[2]=jacobian[2]+d_shapefn[j]*z;
		jacobian[5]=jacobian[5]+d_shapefn[j+8]*z;
		jacobian[8]=jacobian[8]+d_shapefn[j+16]*z;
		}
//__syncthreads();	
////////////////////////////***********Determinant
	determinant=jacobian[0]*(jacobian[4]*jacobian[8]-jacobian[5]*jacobian[7])-jacobian[1]*(jacobian[3]*jacobian[8]-jacobian[5]*jacobian[6])+jacobian[2]*(jacobian[3]*jacobian[7]-jacobian[4]*jacobian[6]);

//	printf("determinant :%f\n",determinant);

/*	cout<<"Jacobian "<<endl<<endl;
	for(int i=0;i<3;i++)
        {
        for(int j=0;j<3;j++)
	cout<<jacobian[i*3+j]<<" ";
	cout<<endl;
	}
*/



////////**************************************inverse of jacobian
{
	float k[9];
	k[0]=(jacobian[4]*jacobian[8]-jacobian[5]*jacobian[7]);
	k[1]=-1*(jacobian[3]*jacobian[8]-jacobian[5]*jacobian[6]);	
	k[2]=(jacobian[3]*jacobian[7]-jacobian[4]*jacobian[6]);
	k[3]=-1*(jacobian[1]*jacobian[8]-jacobian[2]*jacobian[7]);
	k[4]=(jacobian[0]*jacobian[8]-jacobian[2]*jacobian[6]);
	k[5]=-1*(jacobian[0]*jacobian[7]-jacobian[1]*jacobian[6]);
	k[6]=(jacobian[1]*jacobian[5]-jacobian[2]*jacobian[4]);
	k[7]=-1*(jacobian[0]*jacobian[5]-jacobian[2]*jacobian[3]);
	k[8]=(jacobian[0]*jacobian[4]-jacobian[1]*jacobian[3]);
	jacobian[0]=k[0]/determinant;
	jacobian[3]=k[1]/determinant;
	jacobian[6]=k[2]/determinant;
	jacobian[1]=k[3]/determinant;
	jacobian[4]=k[4]/determinant;
	jacobian[7]=k[5]/determinant;
	jacobian[2]=k[6]/determinant;
	jacobian[5]=k[7]/determinant;
	jacobian[8]=k[8]/determinant;

}
/*	cout<<"Inverse of jacobian "<<endl<<endl;
	for(int i=0;i<3;i++)
        {
        for(int j=0;j<3;j++)
	cout<<jacobian[i*3+j]<<" ";
	cout<<endl;
	}

*/

////////////////////////////////////////////////calculating matrix_b
for(int i=0;i<8;i++)
{

	matrix_b[i*3]=jacobian[0]*d_shapefn[i]+jacobian[1]*d_shapefn[i+8]+jacobian[2]*d_shapefn[i+16];
	//cout<<jacobian[3]<<" "<<jacobian[4]<<" "<<jacobian[5]<<""<<endl;
	//cout<<d_shapefn[i]<<" "<<d_shapefn[i+8]<<" "<<d_shapefn[i+16]<<endl;
	matrix_b[i*3+25]=jacobian[3]*d_shapefn[i]+jacobian[4]*d_shapefn[i+8]+jacobian[5]*d_shapefn[i+16];
	matrix_b[i*3+50]=jacobian[6]*d_shapefn[i]+jacobian[7]*d_shapefn[i+8]+jacobian[8]*d_shapefn[i+16];
	matrix_b[i*3+73]=jacobian[6]*d_shapefn[i]+jacobian[7]*d_shapefn[i+8]+jacobian[8]*d_shapefn[i+16];
	matrix_b[i*3+74]=jacobian[3]*d_shapefn[i]+jacobian[4]*d_shapefn[i+8]+jacobian[5]*d_shapefn[i+16];
	matrix_b[i*3+96]=jacobian[6]*d_shapefn[i]+jacobian[7]*d_shapefn[i+8]+jacobian[8]*d_shapefn[i+16];
	matrix_b[i*3+98]=jacobian[0]*d_shapefn[i]+jacobian[1]*d_shapefn[i+8]+jacobian[2]*d_shapefn[i+16];
	matrix_b[i*3+120]=jacobian[3]*d_shapefn[i]+jacobian[4]*d_shapefn[i+8]+jacobian[5]*d_shapefn[i+16];
	matrix_b[i*3+121]=jacobian[0]*d_shapefn[i]+jacobian[1]*d_shapefn[i+8]+jacobian[2]*d_shapefn[i+16];

}
/*	cout<<"\n\n";
	for(int i=0;i<6;i++)
		{for(int j=0;j<24;j++)
		cout<<matrix_b[i*24+j]<<"  ";
		cout<<"\n";}
*/
/////////////////////////////////multiplication B'CB(~~~~~~~another way is possible(suitable for GPU) in which each non zero element is calculated individually without using matrix_b)
for(int i=0;i<24;i++)
{
	float a=0,b=0,c=0,d=0,e=0,f=0;
	a=(matrix_b[i]*d_material_data[0]+matrix_b[i+24]*d_material_data[6]+matrix_b[i+48]*d_material_data[12]+matrix_b[i+72]*d_material_data[18]+matrix_b[i+96]*d_material_data[24]+matrix_b[i+120]*d_material_data[30]);
	b=(matrix_b[i]*d_material_data[1]+matrix_b[i+24]*d_material_data[7]+matrix_b[i+48]*d_material_data[13]+matrix_b[i+72]*d_material_data[19]+matrix_b[i+96]*d_material_data[25]+matrix_b[i+120]*d_material_data[31]);
	c=(matrix_b[i]*d_material_data[2]+matrix_b[i+24]*d_material_data[8]+matrix_b[i+48]*d_material_data[14]+matrix_b[i+72]*d_material_data[20]+matrix_b[i+96]*d_material_data[26]+matrix_b[i+120]*d_material_data[32]);
	d=(matrix_b[i]*d_material_data[3]+matrix_b[i+24]*d_material_data[9]+matrix_b[i+48]*d_material_data[15]+matrix_b[i+72]*d_material_data[21]+matrix_b[i+96]*d_material_data[27]+matrix_b[i+120]*d_material_data[33]);
	e=(matrix_b[i]*d_material_data[4]+matrix_b[i+24]*d_material_data[10]+matrix_b[i+48]*d_material_data[16]+matrix_b[i+72]*d_material_data[22]+matrix_b[i+96]*d_material_data[28]+matrix_b[i+120]*d_material_data[34]);
	f=(matrix_b[i]*d_material_data[5]+matrix_b[i+24]*d_material_data[11]+matrix_b[i+48]*d_material_data[17]+matrix_b[i+72]*d_material_data[23]+matrix_b[i+96]*d_material_data[29]+matrix_b[i+120]*d_material_data[35]);

//	cout<<"a "<<a<<" b "<<b<<" c "<<c<<" d "<<d<<" e "<<e<<" f "<<f<<"\n";
	for(int j=0;j<24;j++)
	{ 
		elemental_data[i*24+j]=elemental_data[i*24+j]+(a*matrix_b[j]+b*matrix_b[j+24]+c*matrix_b[j+48]+d*matrix_b[j+72]+e*matrix_b[j+96]+f*matrix_b[j+120])*determinant*weights[m]*weights[n]*weights[p];
	}

}
//__syncthreads();
/*	cout<<"\n\n";
	for(int i=0;i<24;i++)
		{for(int j=0;j<24;j++)
		cout<<elemental_data[i*24+j]<<"  ";
		cout<<"\n";}
*/
	
	//******************calculating source vector************
/*	for(int i=0;i<4;i++)
	{
		elemental_vector[2*i]+=phi(i+1,gausspoints[m],gausspoints[n],gausspoints[p])*determinant*weights[m]*weights[n]*0;
		elemental_vector[2*i+1]+=phi(i+1,gausspoints[m],gausspoints[n],gausspoints[p])*determinant*weights[m]*weights[n]*0;
	}*/
}
//__syncthreads();
}
//__syncthreads();
}


delete []gausspoints;
delete []weights;
}



/////////////////////////////////////////gives derivative of linear interpolation

/*double d_phi(int i,double point,char type)
{
assert(i>=1 && i<=4);
if(type=='z')
	{
	if(i==1 ||i==4)
		return (-1+pow(-1,i+1)*point)/4;
	return (1+point*pow(-1,i+1))/4;

	}
else
	{
	if(i==1 ||i==2)
		return (-1+point*pow(-1,i+1))/4;
	return (1+point*pow(-1,i+1))/4;	

	}

}*/

__device__ double d_phi(int i,double zeta,double eta,double zi,char type)
	{
	double k=0;
	double ret;
	assert(i>=1 && i<=8);
	switch(type)
	{
		case 'z':
			k=phi(i,zeta,eta,zi);
			if(i==1 || i==4 || i==5 || i==8)	
			ret=-1*k/(1-zeta);
			else
			ret=k/(1+zeta);
			break;
		case 'e':
			k=phi(i,zeta,eta,zi);
			if(i==1 || i==2 || i==5 || i==6)
			ret=-1*k/(1-eta);
			else
			ret=k/(1+eta);
			break;
		case 'j':
			k=phi(i,zeta,eta,zi);
			if(i==1 || i==2 || i==3 || i==4)
			ret=-1*k/(1-zi);
			else
			ret=k/(1+zi);
			break;
  		default :
			ret=0;//MA 	
			break;
	}
	return ret;



	}


///////////////////////////////////gives linear interpolation function
__device__ double phi(int i,double zeta,double eta,double zi)
{
assert(i>=1 && i<=8);
double ret=0;
switch(i)
{
	case 1:
		ret=0.125*(1-zeta)*(1-eta)*(1-zi);
		break;
	case 2:
		ret=0.125*(1+zeta)*(1-eta)*(1-zi);
		break;
	case 3:
		ret= 0.125*(1+zeta)*(1+eta)*(1-zi);
		break;
	case 4:
		ret= 0.125*(1-zeta)*(1+eta)*(1-zi);
		break;
	case 5:
		ret= 0.125*(1-zeta)*(1-eta)*(1+zi);
		break;
	case 6: 
		ret= 0.125*(1+zeta)*(1-eta)*(1+zi);
		break;
	case 7:
		ret= 0.125*(1+zeta)*(1+eta)*(1+zi);
		break;
	case 8:
		ret= 0.125*(1-zeta)*(1+eta)*(1+zi);
		break;
	default:
			
		break;

}
return ret;
}
