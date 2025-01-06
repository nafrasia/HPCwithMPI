//===============================================================
// Matrix-Vector multiplication is performed in parallel using
// MPI. The matrices are stored as boost arrays. This is the
// parallel counterpart of serialMV.cpp code.
//
// Compile: mpic++ -O3 parallelMV.cpp -o parallelMV.exe
//
// Run: mpirun -n <number of cores/cpus> parallelMV.exe
//
// Author: Navid Afrasiabian, 2022
//===============================================================

#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>
#include "../mpistuff.h"

using namespace std;

//============================================================
// Defining a (2d) matrix type and a vector type boost arrays.
// Boost arrays are defined globally since they are the inputs
// of the dotProd function.
// ===========================================================
typedef boost::multi_array<double, 2> mtype;
typedef boost::multi_array<double, 1> atype;

//===========================================================
// dotProd function receives two array types (Arow and b) and
// their size. Then, performs a dot product using for loop.
// ==========================================================
double dotProd(atype &, atype &, int);

int main(int argc,char** argv)
{
	MPI_Stuff mympi(argc,argv); //Initiating MPI

	int debugMode = 0; //This flag can be used to turn the debugging
			   //statements off and on

	//=======================================================
	// Declaring file handles for the matrix and vector files.
	// These files are generated using a serial cpp code.
	//=======================================================
	fstream matfh,vecfh;
	matfh.open("matrixA.txt", ios::in | ios::binary);
	vecfh.open("matrixB.txt", ios::in | ios::binary);

	//========================================================
	// Only 2 variables for matrix size is defined since
	// multiplication is only possible if brows==acols.
	// This condition is checked before calculations are done
	// and the program is aborted if it is not fulfilled
	//========================================================
	int arows, brows,bcols;
	//================================================================	
	// Reading the matrix and vector dimensions from file by the BOSS.
	// Then they are broadcasted to all the cores.
	//================================================================
	if(mympi.nID == 0)
	{
		int acols;
		if(matfh.is_open())
		{
			if(vecfh.is_open())
			{
				matfh >> arows >> acols;
				vecfh >> brows >> bcols;
				if(acols != brows)
				{
					cout << "Column-row mismatch!Cannot multiply" << endl;
					MPI_Abort(MPI_COMM_WORLD,-1);
				}
				if(arows<mympi.nprocs)
				{
					cout << "Too many processors!" << endl;
					MPI_Abort(MPI_COMM_WORLD,-1);
				}
			}
			else
			{
				cout << "Cannot open Vector file!" << endl;
				MPI_Abort(MPI_COMM_WORLD,-1);
			}
		}
		else
		{
			cout << "Cannot open Matrix file!" << endl;
			MPI_Abort(MPI_COMM_WORLD,-1);
			
		}
	}
	MPI_Bcast(&arows,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&brows,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	//=========================================================
	// Vectors b and Arow are defined. b is read from the file
	// by the Boss and broadcasted to everyone. Matrix A and
	// result vector c are defined on the Boss. A is read from
	// file. A row of A is sent to every processor. Then a new
	// row is sent to the next core that is available until
	// all the rows of matrix A are sent. The result is received
	// by the Boss and outputted. The total time spent (on 
	// calculations) is measured on the Boss and outputted.
	//=========================================================
	atype b(boost::extents[brows]);
	atype Arow(boost::extents[brows]);
	if(mympi.nID == 0)
	{
		for(int i = 0;i<brows;i++)
			vecfh >> b[i];
	}
	MPI_Bcast(&b[0],brows,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Status status;
	double Tcalc_all,Tcalc=0.,Ttot,Tcomm;
	if(mympi.nID == 0)
	{
		mtype A(boost::extents[arows][brows]);
		atype c(boost::extents[arows]);
		for(int i = 0;i<arows;i++)
			for(int j = 0;j<brows;j++)
				matfh >> A[i][j];
		double starttime=MPI_Wtime();
		int rowsent = 0;
		for(int i = 1;i<mympi.nprocs;i++){
			MPI_Send(&A[rowsent][0],brows,MPI_DOUBLE,i,rowsent+1,MPI_COMM_WORLD);
			rowsent++;
		}
		for(int i = 0;i<arows;i++)
		{
			double rowans;
			MPI_Recv(&rowans,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			int sender = status.MPI_SOURCE;
			int crow = status.MPI_TAG;
			c[crow-1] = rowans;
			if(rowsent < arows){
				MPI_Send(&A[rowsent][0],brows,MPI_DOUBLE,sender,rowsent+1,MPI_COMM_WORLD);
				rowsent++;
			}
			else
			{
				MPI_Send(MPI_BOTTOM,0,MPI_DOUBLE,sender,0,MPI_COMM_WORLD);
			}
		}
		Ttot = MPI_Wtime()-starttime;
		cout << "T= " << Ttot << endl;
		if (debugMode){
			for(int i = 0;i<arows;i++)
				cout << "c[" << i << "] = " << c[i] << endl;
		}
			
	}
	else
	{
		MPI_Recv(&Arow[0],brows,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		while(status.MPI_TAG != 0)
		{
			double sum;
			double tmp = MPI_Wtime();
			sum = dotProd(Arow,b,brows);
			Tcalc = MPI_Wtime() - tmp;
			
			MPI_Send(&sum,1,MPI_DOUBLE,0,status.MPI_TAG,MPI_COMM_WORLD);
			MPI_Recv(&Arow[0],brows,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}
	}
	MPI_Reduce(&Tcalc,&Tcalc_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(mympi.nID == 0)
	{
		cout << "Tcalc= " << Tcalc_all/mympi.nprocs << endl;
		cout << "Tcomm= " << Ttot - Tcalc_all << endl;
	}

	matfh.close();
	vecfh.close();
	return 0;
}

double dotProd(atype &a1,atype &a2,int bound)
{
	double sum=0.;
        for(int i = 0;i<bound;i++)
		sum += a1[i]*a2[i];

	return sum;
}
