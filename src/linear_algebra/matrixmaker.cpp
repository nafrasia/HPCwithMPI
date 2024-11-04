//==============================================================
// This program creates matrix data files that can be read by
// the other programs and perform matrix-matrix operations on them
//
// Author: Navid Afrasiabian, 2022
//==============================================================

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	//variable declaration
	int arows,brows,bcols;
	double ** A;
	double ** B;

	//fstream declaration and out put opening files
	fstream matfh, vecfh;
	matfh.open("matrixA.txt", ios::out | ios::binary);
	vecfh.open("matrixB.txt", ios::out | ios::binary);
	
	//Receiving the dimensions of the matrices from USER
	cout << "Enter number of A rows: ";
	cin >> arows;
	cout << "Enter number of A columns(B rows): ";
	cin >> brows;
	cout << "Enter number of B columns: ";
	cin >> bcols;
	
	//Allocating A dynamically (non-contiguously)
	A = new double*[arows];
	for(int i = 0;i<arows;i++)
		A[i]= new double[brows];
	
	//Filling A (This current loop creates an Identity matrix)
	for(int i = 0;i<arows;i++)
		for(int j = 0;j<brows;j++){
			if(i==j)
				A[i][j]=1.0;
			else
				A[i][j]=(i+1)*(j+1);
		}

	//===============================================================
	//write the data into the output file. The dimensions are written
	//on top.
	//===============================================================
	matfh << arows;
	matfh << " ";
	matfh << brows;
	for(int i = 0;i<arows;i++){
		matfh << '\n';
		for(int j = 0;j<brows;j++){
			matfh << A[i][j] << '\t';
	}
	}

	delete []A;

	// Repeating a similar process to A for B.
	B = new double *[brows];
	for(int i = 0;i<brows;i++)
		B[i]=new double[bcols];
	
	for(int i = 0;i<brows;i++)
		for(int j = 0;j<bcols;j++){
			B[i][j]=(i+j)%10;
		}
		
	vecfh << brows << " " << bcols;
	for(int i = 0;i<brows;i++){
		vecfh << '\n';
		for(int j = 0;j<bcols;j++){
			vecfh << B[i][j] << '\t';
	}
	}
	delete []B;

	//Closing the files
	matfh.close();
	vecfh.close();
	return 0;
}
