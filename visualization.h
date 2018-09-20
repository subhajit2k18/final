#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "DataTypes.h"



void changeUform(double *sol, REAL4 *u, int *connect_matrix, double *rho, int size, int Nele)
{
	// Copying the nodal displacements
	for(int i=0; i< size; i++)
	{
		u[i].x = sol[i*3];
		u[i].y = sol[i*3+1];
		u[i].z = sol[i*3+2];
	}
	
	// Initializing elemental densities
	for(int i=0; i< Nele; i++)
	{
		u[i].w = 0.0;
	}

	// Copying elemental densities
	for(int i=0; i< Nele; i++)
	{
		u[connect_matrix[i*8]].w = rho[i];
	}
}

// Writing binary file
void WriteBinary(const char sFileName[], const int NX, const int NY, const int NZ, REAL4 *u)
{
	FILE *pFile;
	pFile = fopen (sFileName,"wb");
	if (pFile!=NULL)
	{
		int version = 1;
		int accuracy = sizeof(REAL4);
		fwrite(&version, sizeof(int), 1, pFile);
		fwrite(&NX, sizeof(int), 1, pFile);
		fwrite(&NY, sizeof(int), 1, pFile);
		fwrite(&NZ, sizeof(int), 1, pFile);
		fwrite(&accuracy, sizeof(int), 1, pFile);
		fwrite(u, sizeof(REAL4), NX*NY*NZ, pFile);
		fclose (pFile);
	}
	else
	{
		cerr << "ERROR WRITING OUTPUT FILE: "<<sFileName<<endl;
	}
}

// Reading binary file
bool ReadBinary(const char sFileName[], int *NX, int *NY, int *NZ, REAL4 **u)
{
	ifstream file(sFileName, ios::in|ios::binary);
	if(file.is_open())
	{
		cout << "Inputfile: " << sFileName << " opened sucessfully"<<endl;
		int version;
		file.read((char*)(&version), sizeof(int));
		cout << "Binary File Version: " << version << endl;
		if(version == 1)
		{
			int accuracy;
			file.read((char*)(NX), sizeof(int));
			file.read((char*)(NY), sizeof(int));
			file.read((char*)(NZ), sizeof(int));
			cout << "Data Dimensions: "<<*NX<< " " << *NY <<" "<< *NZ << endl;
			file.read((char*)(&accuracy), sizeof(int));
			cout << "Byte per Datapoint, i.g. REAL4: " << accuracy<<endl;
		
			if(accuracy != sizeof(REAL4))
			{
				cerr << "DATA SIZE MISMATCH! RECOMPILE WITH PROPER ACCURACY." << endl;
				cerr << "SIZE IN DATAFILE IS " << accuracy << " EXPECTED " << sizeof(REAL4)<<endl;
				exit(1);
			}

			REAL4 *u_1;
			u_1 = new REAL4[(*NX)*(*NY)*(*NZ)];
			file.read((char*)(u_1), (*NX)*(*NY)*(*NZ)*accuracy);
			file.close();
			*u = u_1;
			cout << "Read Binary complete"<<endl;
			return false;
		}
		else
		{
			cerr << "VERSION MISMATCH"<<endl;
			return true;
		}
	}
	else
	{
		cerr << "ERROR READING INPUT FILE: "<<sFileName<<endl;
		return true;
	}
}

// Writing .vti file

bool WriteVTI4DPadded(const char sFileName[], int NX, int NY, int NZ, const REAL4 u[])
{
	//Add ending ".vts"
	char sMyFileName[2048];
	sprintf(sMyFileName, "%s.vti", sFileName);
	printf("Writing output file %s\n", sMyFileName);
	FILE *fOutput = fopen(sMyFileName, "w");
	if(fOutput == NULL)
	{
		printf("ERROR: could not open file %s\n", sMyFileName);
		return true;
	}
	REAL4 *u_new = new REAL4[(NX+2)*(NY+2)*(NZ+2)];
	/*
	for(int k=0;k<NZ+2;k++)
	{
		for(int j=0;j<NY+2;j++)
		{
			for(int i=0;i<NX+2;i++)
			{
				u_new[i+(NX+2)*j+(NX+2)*(NY+2)*k].x = 0.0;
				u_new[i+(NX+2)*j+(NX+2)*(NY+2)*k].y = 0.0;
				u_new[i+(NX+2)*j+(NX+2)*(NY+2)*k].z = 0.0;
				u_new[i+(NX+2)*j+(NX+2)*(NY+2)*k].w = 0.0;
			}
		}
	}
	 */
	//copy boundary into halo
	for(int k=0;k<NZ;k++)
	{
		for(int j=0;j<NY;j++)
		{
			int i = 0;
			u_new[i+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
			i = NX-1;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
		}
	}
	for(int k=0;k<NZ;k++)
	{
		for(int i=0;i<NX;i++)
		{
			int j = 0;
			u_new[i+1+(NX+2)*(j)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+1+(NX+2)*(j)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+1+(NX+2)*(j)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+1+(NX+2)*(j)+(NX+2)*(NY+2)*(k+1)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
			j = NY-1;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
		}
	}
	for(int j=0;j<NY;j++)
	{
		for(int i=0;i<NX;i++)
		{
			int k=0;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
			k=NZ-1;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
			u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].w = 0.0;//u[i+NX*j+NX*NY*k].w;
		}
	}
	//copy values
	for(int k=0;k<NZ;k++)
	{
		for(int j=0;j<NY;j++)
		{
			for(int i=0;i<NX;i++)
			{
				u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].x = u[i+NX*j+NX*NY*k].x;
				u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].y = u[i+NX*j+NX*NY*k].y;
				u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].z = u[i+NX*j+NX*NY*k].z;
				u_new[i+1+(NX+2)*(j+1)+(NX+2)*(NY+2)*(k+1)].w = u[i+NX*j+NX*NY*k].w;
			}
		}
	}
	NX += 2;
	NY += 2;
	NZ += 2;
	//write header
	fprintf(fOutput, "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fOutput, "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", 0, NX-1, 0, NY-1, 0, NZ-1);
	fprintf(fOutput, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, NX-1, 0, NY-1, 0, NZ-1);
	fprintf(fOutput, "      <PointData Scalars=\"ScalarPointData\">\n");
	fprintf(fOutput, "        <DataArray type=\"Float32\" Name=\"Displacements\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	
	for(int k=0;k<NZ;k++)
	{
		for(int j=0;j<NY;j++)
		{
			for(int i=0;i<NX;i++)
			{
				const int Index = i + j*NX + k*NX*NY;
				REAL Data1 = u_new[Index].x;
				REAL Data2 = u_new[Index].y;
				REAL Data3 = u_new[Index].z;
				if(fabs(Data1) < 1e-12) Data1 = 0.0;
				if(fabs(Data2) < 1e-12) Data2 = 0.0;
				if(fabs(Data3) < 1e-12) Data3 = 0.0;
				fprintf(fOutput, "%e %e %e\n", Data1, Data2, Data3);
			}
		}
	}
	fprintf(fOutput, "        </DataArray>\n");
	/*
	fprintf(fOutput, "        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	
	for(int k=0;k<NZ;k++)
	{
		for(int j=0;j<NY;j++)
		{
			for(int i=0;i<NX;i++)
			{
				const int Index = i + j*NX + k*NX*NY;
				fprintf(fOutput, "%e\n", u_new[Index].w);
			}
		}
	}
	fprintf(fOutput, "        </DataArray>\n");
	*/
	fprintf(fOutput, "      </PointData>\n");
	fprintf(fOutput, "      <CellData>\n");
	fprintf(fOutput, "        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
	for(int k=0;k<NZ-1;k++)
	{
		for(int j=0;j<NY-1;j++)
		{
			for(int i=0;i<NX-1;i++)
			{
				const int Index = i + j*NX + k*NX*NY;
				REAL Data = u_new[Index].w;
				if(fabs(Data) < 1e-12) Data = 0.0;
				fprintf(fOutput, "%e\n", Data);
			}
		}
	}
	fprintf(fOutput, "        </DataArray>\n");
	fprintf(fOutput, "      </CellData>\n");
	fprintf(fOutput, "    </Piece>\n");
	fprintf(fOutput, "  </ImageData>\n");
	fprintf(fOutput, "</VTKFile>\n");
	fclose(fOutput);
	delete [] u_new;
	return false;
}



