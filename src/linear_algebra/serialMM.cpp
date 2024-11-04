//======================================================================
// Matrix-Matrix multiplication is performed on boost arrays.
// boost arrays are very similar to static C arrays but allow
// for bound-checking and safer programming. The size of boost
// arrays can be adjusted after initial allocation using resize()
// function. To use boost array you might need to install
// the boost libraries using:
// 
// sudo apt install libboost-all-dev (Linux inc. WSL)
//
// Compile (with optimization): g++ -O3 serialMM.cpp -o serialMM.exe
//
// Run: ./serialMM.exe
//
// Author: Navid Afrasiabian, 2022
//======================================================================

#include <iostream>
#include "boost/multi_array.hpp"
#include <fstream>
#include <chrono>
using namespace std;


int main()
{
	int debugMode = 0; //This flag can be used to output
			   //debugging print statements
			   
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
			//Allocating boost arrays
			typedef boost::multi_array<double, 2> mtype;
			
			mtype A(boost::extents[arows][acols]);
			mtype b(boost::extents[brows][bcols]);
			mtype c(boost::extents[arows][bcols]);

			//Reading and storing the data from files
			for(int i=0;i<arows;i++)
				for(int j = 0;j<brows;j++)
					matfh >> A[i][j];
			for(int i = 0;i<brows;i++)
				for(int j=0;j<bcols;j++)
					vecfh >> b[i][j];
			//=============================================
			// Performing Matrix-Matrix multiplication.
			// Timer is started. The performance of the
			// program is tested by measure the time it
			// takes to perform the MM multiplication
			//=============================================
			auto start = std::chrono::steady_clock::now();
				for(int i = 0;i<arows;i++){
					for(int j=0;j<bcols;j++){
						c[i][j]=0.;
						for(int k=0;k<brows;k++)
							c[i][j]+=A[i][k]*b[k][j];
					}
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
