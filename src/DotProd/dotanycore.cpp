//=========================================================
// This code Takes two vectors and calculates their dot
// product.
//
// Author: Navid Afrasiabian, 2022
// ========================================================
#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"

//Declaration and defining functions and classes starts here

class MPI_stuff 
{
public:
  int NProcs;
  int MyID;

  MPI_stuff(int &argc, char** &argv)
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
  }

  ~MPI_stuff()
  {
    MPI_Finalize();
  }
}; 


// Get the number of rows of data in the file from first row of file
int GetNumberElements(std::ifstream &fin)
{
  int n;

  if (fin.is_open())
    fin >> n;
  else { // no file to read from
    std::cout << "Input File not found\n";
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  return n;
}

// Read n elemnts of two column data from a file, store in arrays a,b
void ReadArrays(std::ifstream &fin, std::vector<float> &a, std::vector<float> &b, int n)
{
  if(fin.is_open()){
    for(int i = 0; i < n; i++)
      fin >> a[i] >> b[i];
  }
  else{ // fp null means no file to read from
    std::cout << "Input File not found\n";
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
}

//=======================================================================================
// The DotProducts takes vectors a and b, total number of elements, and class MPI_stuff as
// input. The size of the worker portion of the arrays is saved in pSize (for passing 
// size). The Boss node first sends the workers data to them, then calculates the sum of 
// the remainder of the data set.
// =======================================================================================
float DotProduct(std::vector<float> &a, std::vector<float> &b, 
		  int Ntot, const MPI_stuff &the_mpi)
{ 
  MPI_Status status;

  int pSize = Ntot/(the_mpi.NProcs-1);
  float sum=0;
 
  // Distribute or Scatter the array to the workers
  if(the_mpi.MyID == 0){
    int size = Ntot % (the_mpi.NProcs-1);
    for(int i = 1; i < the_mpi.NProcs; i++){
      MPI_Send(&a[pSize*(i-1)],pSize,MPI_FLOAT,i,10,MPI_COMM_WORLD);
      MPI_Send(&b[pSize*(i-1)],pSize,MPI_FLOAT,i,20,MPI_COMM_WORLD);
    }
    // Work out my part of the sum
    for(int i = Ntot-size; i < Ntot; i++)
	sum += a[i] * b[i];
  }
  else{
    MPI_Recv(&a[0],pSize,MPI_FLOAT,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(&b[0],pSize,MPI_FLOAT,0,20,MPI_COMM_WORLD,&status);
    // Work out my part of the sum
    for(int i = 0; i < pSize; i++)
    	sum += a[i] * b[i];
  }


  return sum;
}

//End of declaration and defining. Main body starts.

int main(int argc, char** argv)
{
  MPI_stuff the_mpi(argc, argv);
  
  std::ifstream fin;
  int n = 0;
  int size;
  float adotb,Gdot;

  if(the_mpi.MyID == 0) {
    fin.open("DotData.txt");
    n = GetNumberElements(fin);
  }
  //=========================================================================
  // An if condition checks if only one core has been requested. Since there
  // is no sending or receiving in this case, we just carry out the summation
  // in series. After defining the vectors, the ReadArrays function is used
  // to save the data from the file in to the vectors. A for loop goes over
  // all the data and sums the terms up into Gdot.
  // ========================================================================
  if (the_mpi.NProcs == 1)
  {
	std::vector<float> a,b;
	a.reserve(n); a.resize(n); // try to avoid std::vector reallocations
	b.reserve(n); b.resize(n);
    	
	ReadArrays(fin,a,b,n);
  	Gdot = 0;
	for (int i = 0; i < n ; i++)
		Gdot += a[i]*b[i]; 
  }
  else{
	//============================================================================
	// In the case of more than 1 core, the number of data points is broadcasted
	// to all the nodes. The result of summation of different portions of the
	// data is saved in adotb. Then, the total sum is calculated using a MPI_Reduce
	// command.
	// ==========================================================================
  	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);


  	std::vector<float> a,b;
  	if(the_mpi.MyID == 0){
    		a.reserve(n);  a.resize(n); 
    		b.reserve(n);  b.resize(n);
    
    		ReadArrays(fin,a,b,n);
  	}
  	else {
    		a.reserve(size);  a.resize(size);
    		b.reserve(size);  b.resize(size); 
  	}
	
  	adotb = DotProduct(a,b,n,the_mpi);

  	// Collect individual sums in Gsum
  	MPI_Reduce(&adotb,&Gdot,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
  }
  if(the_mpi.MyID == 0)
    std::cout << "The inner product is " << Gdot << std::endl;

  return 0;
}
