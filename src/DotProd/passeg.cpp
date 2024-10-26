//=========================================================
// This code Takes two vectors and calculates their dot
// product.
//
// Author: Navid Afrasiabian, 2022
// ========================================================
#include <iostream>
#include <fstream>
#include <vector>
#include "../mpistuff.h"

using namespace std;
//=========================================================================
// The dotProduct function takes two vectors and their size and calculates
// the dot product of the vectors. Caveat: This program doesn't check the
// length of the two vectors to ensure they are the same size.
// 
// The readArray function takes a file handle, two vectors and the number
// of elements to read. It reads that number of elements from the file and
// saves them in the vectors.
//
// The functions are declared here and defined after the main body.
// ========================================================================
float dotProduct(std::vector<float> &, std::vector<float> &, int);
void readArray(fstream &, std::vector<float> &, std::vector<float> &, int);

// main body starts here
int main(int argc, char** argv)
{
	MPI_Stuff mympi(argc, argv);    //defining an instance of class MPI_Stuff
	MPI_Status status;		//MPI_Status variable for MPI_Recv
	float adotb,Gdot;		//local and global dot product variables
	int n,rsize, psize;		//n to store the total number of data points,
					//rsize for size of arrays on Boss core
					//psize for size of arrays on Worker cores
	
	fstream dotfh;
	dotfh.open("DotData.txt",ios::in | ios::binary);
	//==========================================================================
	// The Boss core checks if the file was opened successfully and otherwise
	// aborts the MPI processes. The number of data points is divided by the
	// number of Worker cores to determine the size of their portion of data
	// to handle. If not divisible, the remainder of the data will be handled by
	// the Boss core. The portion size is broadcasted to all cores.
	// The Boss core reads segments of the data file, sends it to corresponding
	// worker core until all Workers have received their segment. The Boss then
	// finds the dot product of the remainder if required.
	// =========================================================================
	if (mympi.nID == 0)
	{
		if (dotfh.is_open())
		{
			dotfh >> n;
			rsize = n % (mympi.nprocs-1);
			psize = n/(mympi.nprocs-1);
			
			MPI_Bcast(&psize,1,MPI_INT,0,MPI_COMM_WORLD);

			std::vector<float> a,b;
			a.reserve(psize);a.resize(psize);
			b.reserve(psize);b.resize(psize);
			
			for (int i =1;i <mympi.nprocs;i++)
			{
				readArray(dotfh,a,b,psize);
				MPI_Send(&a[0],psize,MPI_FLOAT,i,10,MPI_COMM_WORLD);
				MPI_Send(&b[0],psize,MPI_FLOAT,i,20,MPI_COMM_WORLD);
			}	
			if (rsize)
			{
				a.reserve(rsize);a.resize(rsize);
				b.reserve(rsize);b.resize(rsize);
			
				readArray(dotfh,a,b,rsize);
				adotb = dotProduct(a,b,rsize);
			}
			else
				adotb = 0;

		}
		else{
			cout << "Cannot open file!" << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else
	{
		std::vector<float> a,b;
		a.reserve(psize);a.resize(psize);
		b.reserve(psize);b.resize(psize);
		
		MPI_Recv(&a[0],psize,MPI_FLOAT,0,10,MPI_COMM_WORLD,&status);
		MPI_Recv(&b[0],psize,MPI_FLOAT,0,20,MPI_COMM_WORLD,&status);	
		
		adotb = dotProduct(a,b,psize);
	}
	//collecting all the local dot product and sum them up to find the total.
	MPI_Reduce(&adotb,&Gdot,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
	
	// The Boss prints out the result
	if (mympi.nID == 0)
		cout << "Dot product = " << Gdot << endl;
	dotfh.close();
	
	return 0;
}

//main body ends here. Function defining starts here.
//
float dotProduct(std::vector<float> &a,std::vector<float> &b,int size)
{	
	float sum=0;
	for (int i=0;i<size;i++)
		sum += a[i]*b[i];
	
	return sum;
}

void readArray(fstream & dotfile, std::vector<float> &a, std::vector<float> &b,int size)
{	
	for (int i = 0; i < size;i++)
		dotfile >> a[i] >> b[i];
}
