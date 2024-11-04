//=====================================================================
// Matrix-Matrix multiplication is performed using the BLAS library.
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
// Compile: g++ serialBLAS.cpp -lblas -o serialBLAS.exe
//
// Run: ./serialBLAS.exe
//
// Author: Navid Afrasiabian, 2022
//=====================================================================

#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>
#include <chrono>
using namespace std;

//============================================================
// Defining a (2d) matrix type and a vector type boost arrays.
// Boost arrays are defined globally since they are the inputs
// of the dotProd function.
// ===========================================================
typedef boost::multi_array<double, 2> mtype;
typedef boost::multi_array<double, 1> atype;

//=============================================================================
// BLAS functions need to be declared as an extern as they are written in C.
// The dgemv_ function is a "Double-precision GEneral Matrix Vector" multiplier.
// For more information, see "https://en.cppreference.com/w/cpp/numeric/linalg"
//=============================================================================
extern "C"{extern void dgemv_(char*, int*, int*,double*,double*,int*,double*,int*,double*,double*,int*);}

int main()
{
	int debugMode = 0;//debugMode = 1 activates print statements
			  //that can be used for debugging

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
        int arows, acols, brows, bcols;
	if(matfh.is_open())
	{
		if(vecfh.is_open())
		{
			matfh >> arows >> acols;
			vecfh >> brows >> bcols;
			if (brows!=acols)
			{
				cout << "Dimensions don't match!" << endl;
				exit(1);
			}
			mtype A(boost::extents[arows][acols]);
			mtype b(boost::extents[brows][bcols]);
			mtype c(boost::extents[arows][bcols]);

			//declaring dgemv_ argument variables
			double alpha= 1., beta=0.;
                        int incA=1,incC=1;
			char trans = 'T';
			for(int i=0;i<arows;i++)
				for(int j = 0;j<brows;j++)
					matfh >> A[i][j];
			for(int i = 0;i<brows;i++)
				for(int j=0;j<bcols;j++)
					vecfh >> b[i][j];
			auto start = std::chrono::steady_clock::now();
				for(int i = 0;i<arows;i++){
					dgemv_(&trans,&brows,&bcols,&alpha,&b[0][0],&bcols,&A[i][0],&incA,&beta,&c[i][0],&incC);//BLAS functions reference by address. Hence, & in the arguments
					}
			auto end = std::chrono::steady_clock::now();
			std::cout << "Total Wall time: " << std::chrono::duration<double> (end-start).count() << std::endl;
			if (debugMode){
			for(int i = 0;i<arows;i++)
				for(int j=0;j<bcols;j++)
					cout << "c[" << i << "," << j << "] " << c[i][j] << endl;
			}
			matfh.close();
			vecfh.close();

		}
		else
		{
			cout << "Cannot open vector file!" << endl;
			exit(1);
		}
	}
	else
	{
		cout << "Cannot open matrix file!" << endl;
		exit(1);
	}
	return 0;
}
