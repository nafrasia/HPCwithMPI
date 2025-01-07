//=================================================================
// Matrix-Matrix multiplication is performed in parallel using
// MPI. The matrices are stored in boost arrays. To use boost arrays, 
// the boost libraries must be installed:
//
// sudo apt install libboost-all-dev (Linux inc. WSL)
//
// This program is the parallel counterpart of serialMM.cpp.
//
// Compile: mpic++ -O3 parallelMM.cpp -o parallelMM.exe
//
// Run: mpirun -n <number of cores> parallelMM.exe
//
// Author: Navid Afrasiabian, 2022
//=================================================================
#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>
#include "../mpistuff.h"

using namespace std;

//=======================================================
// A matrix type and an array type boost arrays are
// defined globally. These are inputs of the mrowProd
// function.
//=======================================================
typedef boost::multi_array<double, 2> mtype;
typedef boost::multi_array<double, 1> atype;

//==============================================================
// mrowProd receives a boost array and boost matrix and
// performs a matrix-vector multiplication using a nested 
// for loop (2 loops since 2d). The dimensions of the matrix
// and the vector are given to this function as input.
// A pointer to the result array is also passed to the function.
//==============================================================
void mrowProd(atype & ,mtype & ,double *,  int , int);


int main(int argc,char** argv)
{
	MPI_Stuff mympi(argc,argv);
	
	//=======================================================
        // Declaring file handles for the matrix files.
        // These files are generated using a serial cpp code.
        //=======================================================	
	fstream matAfh,matBfh;
	matAfh.open("matrixA.txt", ios::in | ios::binary);
	matBfh.open("matrixB.txt", ios::in | ios::binary);
	
	//========================================================
        // Only 3 variables for matrix size is defined since
        // multiplication is only possible if brows==acols.
        // This condition is checked before calculations are done
        // and the program is aborted if it is not fulfilled
        //========================================================	
	int arows, brows,bcols;
	
	//================================================================
        // Reading the dimensions of the matrices from file by the Boss.
        // Then they are broadcasted to all the cores.
        //================================================================	
	if(mympi.nID == 0)
	{
		int acols;
		if(matAfh.is_open())
		{
			if(matBfh.is_open())
			{
				matAfh >> arows >> acols;
				matBfh >> brows >> bcols;
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
	MPI_Bcast(&bcols,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//============================================================
        // Matrix B and vector Arow are defined. B is read from a file
        // by the Boss and broadcasted to everyone. Matrix A and
        // result Matrix C are defined on the Boss. A is read from
        // file. A row of A is sent to every processor. Then a new
        // row is sent to the next core that is available until
        // all the rows of matrix A are sent. The result is received
        // by the Boss and outputted. The total time spent (on
        // calculations) is measured on the Boss and outputted.
        //=========================================================
	mtype B(boost::extents[brows][bcols]);
	atype Arow(boost::extents[brows]);
	if(mympi.nID == 0)
	{
		for(int i = 0;i<brows;i++)
			for(int j = 0;j<bcols;j++)
				matBfh >> B[i][j];
	}
	MPI_Bcast(&B[0][0],brows*bcols,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Status status;
	double Tcalc_all,Tcalc=0.,Ttot,Tcomm;
	if(mympi.nID == 0)
	{
		mtype A(boost::extents[arows][brows]);
		mtype C(boost::extents[arows][bcols]);
		for(int i = 0;i<arows;i++)
			for(int j = 0;j<brows;j++)
				matAfh >> A[i][j];
		double starttime=MPI_Wtime();
		int rowsent = 0;
		for(int i = 1;i<mympi.nprocs;i++){
			MPI_Send(&A[rowsent][0],brows,MPI_DOUBLE,i,rowsent+1,MPI_COMM_WORLD);
			rowsent++;
		}
		for(int i = 0;i<arows;i++)
		{
			double * rowans;
			rowans = new double[bcols];
			MPI_Recv(&rowans[0],bcols,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			int sender = status.MPI_SOURCE;
			int crow = status.MPI_TAG;
			for(int j = 0;j<bcols;j++){
				C[crow-1][j] = rowans[j];
			}
			if(rowsent < arows){
				MPI_Send(&A[rowsent][0],brows,MPI_DOUBLE,sender,rowsent+1,MPI_COMM_WORLD);
				rowsent++;
			}
			else
			{
				MPI_Send(MPI_BOTTOM,0,MPI_DOUBLE,sender,0,MPI_COMM_WORLD);
			}
		
			delete []rowans;
		}
		Ttot = MPI_Wtime()-starttime;
		cout << "T= " << Ttot << endl;
		if (debugMode){
			for(int i = 0;i<arows;i++)
				for(int j=0;j<bcols;j++)
					cout << "c[" << i << "," <<j << "] = " << C[i][j] << endl;
		}
			
	}
	else
	{
		MPI_Recv(&Arow[0],brows,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		while(status.MPI_TAG != 0)
		{
			double * Crow;
			Crow = new double[bcols];
			double tmp = MPI_Wtime();
			mrowProd(Arow,B,&Crow[0],brows,bcols);
			Tcalc = MPI_Wtime() - tmp;
			
			MPI_Send(&Crow[0],bcols,MPI_DOUBLE,0,status.MPI_TAG,MPI_COMM_WORLD);
			MPI_Recv(&Arow[0],brows,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			delete []Crow;
		}
	}
	MPI_Reduce(&Tcalc,&Tcalc_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(mympi.nID == 0)
	{
		cout << "Tcalc= " << Tcalc_all/mympi.nprocs << endl;
		cout << "Tcomm= " << Ttot - Tcalc_all << endl;
	}

	matAfh.close();
	matBfh.close();
	return 0;
}

void mrowProd(atype &a,mtype &b,double *c,  int nrows, int ncols)
{

	for(int i = 0; i<ncols;i++){
		c[i]=0.0;
		for(int j = 0; j<nrows;j++){
			c[i]+=a[j]*b[j][i];
		}
	}
}
