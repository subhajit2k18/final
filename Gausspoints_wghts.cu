#include<cstdlib>
#include<stdio.h>
using namespace std;

__device__ void getGausspoints(double *gausspoints,double *weights,int nint)
{
//cout<<"nint "<<nint<<endl;
switch (nint)
{
	case 1:
	{
	 double xk[1]={0.000000000000000000000000};
	 double wk[1]={2.000000000000000000000000};
	gausspoints[0]=xk[1];
	weights[0]=wk[0];
	
	}
	break;
	case 2:
	{ double x[1] = {0.5773502691896257645091488};
	 double w[1] = {1.0000000000000000000000000};
	for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}}
	break;
	case 4:
	{ double x[2] = {0.3399810435848562648026658,0.8611363115940525752239465};
	 double w[2] = {0.6521451548625461426269361,0.3478548451374538573730639};
	for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}
	}
	break;
	case 6:
	{ double x[3] = {0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016};
	 double w[3] = {0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961};
	for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}
	}
	break;
	case 8:
	{ double x[4] = {0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609};
	 double w[4] = {0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314};
		for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}
	}
	break;
	case 10:
	{ double x[5] = {0.1488743389816312108848260,0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640};
	 double w[5] = {0.2955242247147528701738930,0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688};
		for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}
	}
	break;
	case 12:
	{ double x[6] = {0.1252334085114689154724414,0.3678314989981801937526915,0.5873179542866174472967024,0.7699026741943046870368938,0.9041172563704748566784659,0.9815606342467192506905491};
	 double w[6] = {0.2491470458134027850005624,0.2334925365383548087608499,0.2031674267230659217490645,0.1600783285433462263346525,0.1069393259953184309602547,0.0471753363865118271946160};
		for(int i=0;i<nint/2;i++)
	{
		gausspoints[2*i]=x[i];
		gausspoints[2*i+1]=-x[i];
		weights[2*i]=w[i];
		weights[2*i+1]=w[i];
	}
	}
	break;
	case 3:
	{ double x[2] = {0.0000000000000000000000000,0.7745966692414833770358531};
	 double w[2] = {0.8888888888888888888888889,0.5555555555555555555555556};
	gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	case 5:
	{ double x[3] = {0.0000000000000000000000000,0.5384693101056830910363144,0.9061798459386639927976269};
	 double w[3] = {0.5688888888888888888888889,0.4786286704993664680412915,0.2369268850561890875142640};
	gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	case 7:
	{ double x[4] = {0.0000000000000000000000000,0.4058451513773971669066064,0.7415311855993944398638648,0.9491079123427585245261897};
	 double w[4] = {0.4179591836734693877551020,0.3818300505051189449503698,0.2797053914892766679014678,0.1294849661688696932706114};
	gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	case 9:
	{ double x[5] = {0.0000000000000000000000000,0.3242534234038089290385380,0.6133714327005903973087020,0.8360311073266357942994298,0.9681602395076260898355762};
	 double w[5] = {0.3302393550012597631645251,0.3123470770400028400686304,0.2606106964029354623187429,0.1806481606948574040584720,0.0812743883615744119718922};
gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	case 11:
	{ double x[6] = {0.0000000000000000000000000,0.2695431559523449723315320,0.5190961292068118159257257,0.7301520055740493240934163,0.8870625997680952990751578,0.9782286581460569928039380};
	 double w[6] = {0.2729250867779006307144835,0.2628045445102466621806889,0.2331937645919904799185237,0.1862902109277342514260976,0.1255803694649046246346943,0.0556685671161736664827537};
	gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	case 13:
	{ double x[7] = {0.0000000000000000000000000,0.2304583159551347940655281,0.4484927510364468528779129,0.6423493394403402206439846,0.8015780907333099127942065,0.9175983992229779652065478,0.9841830547185881494728294};
	 double w[7] = {0.2325515532308739101945895,0.2262831802628972384120902,0.2078160475368885023125232,0.1781459807619457382800467,0.1388735102197872384636018,0.0921214998377284479144218,0.0404840047653158795200216};
	gausspoints[0]=x[0];weights[0]=w[0];
	for(int i=0;i<(nint-1)/2;i++)
	{
	gausspoints[2*i+1]=x[i+1];
	gausspoints[2*i+2]=-x[i+1];
	weights[2*i+1]=w[i+1];
	weights[2*i+2]=w[i+1];
	}
	}
	break;
	default:
		{printf("NINT not found\n");	
	//	exit(0);
		asm("trap;");
		}
}







}

