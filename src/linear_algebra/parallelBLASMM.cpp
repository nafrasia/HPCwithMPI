
//=====================================================================
// Matrix-Matrix multiplication is performed in parallel using
// the BLAS library and MPI.
// The matrix-vector multiplication function of the BLAS library replaced
// the inner loops of serialMM.cpp program. The peformance of BLAS can
// be compared to basic C++ operators by comparing this program and
// serialMM.cpp. The rest of the program is similar.
//
// Keep in mind that BLAS and LAPACK are column major!!!
//
// Install LAPACK and BLAS using:
//
// sudo apt install liblapack-dev
//
// Compile: mpic++ serialBLAS.cpp -lblas -o parallelBLASMM.exe
//
// Run: mpirun -n <number of cores> ./parallelBLASMM.exe
//
// Author: Navid Afrasiabian, 2022
//=====================================================================
#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>
#include "mpistuff.h"

using namespace std;


typedef boost::multi_array<double, 2> mtype;
typedef boost::multi_array<double, 1> atype;

extern "C"{extern void dgemv_(char*, int*, int*,double*,double*,int*,double*,int*,double*,double*,int*);}

int main(int argc,char** argv)
{
	MPI_Stuff mympi(argc,argv);
	
	fstream matAfh,matBfh;
	matAfh.open("matrixA.txt", ios::in | ios::binary);
	matBfh.open("matrixB.txt", ios::in | ios::binary);
	int arows, brows,bcols;
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
		for(int i = 0;i<arows;i++)
			for(int j=0;j<bcols;j++)
				cout << "c[" << i << "," <<j << "] = " << C[i][j] << endl;
			
	}
	else
	{
		MPI_Recv(&Arow[0],brows,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		while(status.MPI_TAG != 0)
		{
			double * Crow;
			Crow = new double[bcols];
			char trans='T';
			double alpha= 1., beta=0.;
			int incA=1,incC=1;
			double tmp = MPI_Wtime();
			dgemv_(&trans,&brows,&bcols,&alpha,&B[0][0],&bcols,&Arow[0],&incA,&beta,&Crow[0],&incC);
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
