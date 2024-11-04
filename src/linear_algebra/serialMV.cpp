//===========================================================
// Matrix-Vector multiplication is performed on boost arrays.
// The matrices and vectors are generated using matrixmaker.cpp
// program.
//
// Author: Navid Afrasiabian, 2022 
//===========================================================
#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>

using namespace std;

int main()
{
	int debugMode = 0;//This flag can be used to
			  //turn debugging statements off and on

	//Importing data
	fstream matfh, vecfh;
	matfh.open("matrixA.txt", ios::in | ios::binary);
	vecfh.open("vectorb.txt", ios::in | ios::binary);
	
	//==================================================
	// For matrix-vector multiplication, we need
	// the dimensions of the matrix and vector to match.
	// Assuming right matrix multiplication, we require
	// the number of rows of the vector and columns
	// of the matrix to be the same.
	//==================================================
	if (matfh.is_open())
	{
		if(vecfh.is_open())
		{
			int arows,acols,brows;
			matfh >> arows >> acols;
			vecfh >> brows;
			if (acols != brows)
			{
				cout << "Matrix columns doesn't match vector rows!" << endl;
				matfh.close();
				vecfh.close();
				exit(-1);
			}

			//declaring boost arrays
			typedef boost::multi_array<double, 2> mtype;
			mtype A(boost::extents[arows][acols]);
			typedef boost::multi_array<double, 1> atype;
			atype b(boost::extents[brows]);
			typedef boost::multi_array<double, 1> atype;
			atype c(boost::extents[arows]);
			
			//Reading in the data and storing them in 
			//the boost arrays declared above
			for(int i = 0;i<arows;i++){
				for(int j = 0;j<acols;j++){
					matfh >> A[i][j];
				}
			}
			for(int i = 0;i<brows;i++){
				vecfh >> b[i];
			}
			//Performing matrix-vector multiplication
			for(int i = 0;i<arows;i++){
				c[i]=0;
				for(int j = 0;j<acols;j++)
					c[i] += A[i][j]*b[j];
			}
			if (debugMode){
				for(int i = 0;i<arows;i++)
					cout << "c" << i+1 << "= " << c[i] << endl;
			}
		}
		else
		{
			cout << "Cannot open vector file" << endl;
			matfh.close();
			vecfh.close();
			exit(-1);
		}
	}
	else
	{
		cout << "Cannot open matrix file" << endl;
		matfh.close();
		vecfh.close();
		exit(-1);
	}
	matfh.close();
	vecfh.close();
	return 0;	
}
