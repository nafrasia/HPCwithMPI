//===========================================================================
// Moments of a distribution are computed using Multiple-Instruction
// Single-Data (SIMD) parallelization. Each core has access to the
// entire dataset but they compute a different moment. For ease of 
// validity check, I created my dataset by sampling from a normal
// distribution, moments of which are known and analytically 
// computable.
//
// Compile by: mpic++ -O2 misd_moments.cpp -o <executable file name>
//
// Run the executable by: mpirun -n <number of cores> <executable file name>
//
// Author: Navid Afrasiabian, 2022
// ==========================================================================
#include <iostream>
#include "../mpistuff.h"
#include <random>
#include <cmath>
#include <chrono>
#include <ctime>

double momentSum(int, double *, int );

int main(int argc, char** argv)
{
	MPI_Stuff mympi(argc,argv); //Initializing MPI class
	
	double p[1000000];	 // Data Array 
	const int nrolls=1000000;  // number of experiments
	double mpiStartTime = MPI_Wtime();	
	double all_time;
	auto start_time = std::chrono::system_clock::now();
	//=============================================================
	// Using the random library, the boss core draws numbers from
	// a normal distribution of mean value 5 and standard deviation
	// of 2 and stores them in p. Then p is sent to all other cores
	// using broadcast command.
	// ============================================================
	if (mympi.nID == 0)
	{
		std::cout << "Number of data points: " << nrolls << std::endl;
  	  
		std::default_random_engine generator;
  	  	std::normal_distribution<double> distribution(5.0,2.0);

  	  	for (int i=0; i<nrolls; ++i) 
    			p[i] = distribution(generator);
	}
	
	MPI_Bcast(&p,nrolls,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//=========================================================================
	//momentSum function takes the rank/id of the core, the dataset, and the 
	//total number of data points and calculates the id-th moment of the data.
	//Therefore, we expect a result like this:
	//0-th moment = 1.0
	//1-th moment = mean = 5.0
	//2-th moment = std^2 + mean^2 = 29
	//3-th moment = mean^3 + 3*std^2*mean = 185
	//  .
	//  .
	//  .
	//  ======================================================================
	printf("%d-th moment is %f\n", mympi.nID, momentSum(mympi.nID, p, nrolls));
	MPI_Barrier(MPI_COMM_WORLD);	
	//=========================================================================
	// The timing and performance of the code is tested using both MPI and C++
	// functions. The system_clock routine from chrono library of the C++
	// standard libraries and the MPI_Wtime() from MPI.
	// An MPI_Reduce() is used to sum the time spend by each core. The boss core
	// prints out the total time spent.
	// ======================================================================== 	
	auto end_time = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_time = (end_time - start_time);
	std::cout << "Elapsed Time: " << elapsed_time.count() << " for core " << mympi.nID << std::endl;
	MPI_Reduce(&elapsed_time, &all_time, 1, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	double mpiElapsed = MPI_Wtime()-mpiStartTime;
	double mpiAllTime;	
	MPI_Reduce(&mpiElapsed, &mpiAllTime, 1, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if (mympi.nID == 0){
		std::cout << "total core time: " << all_time << std::endl;
		std::cout << "MPI time: " << mpiAllTime << std::endl;
	}
	
	return 0;

}

double momentSum(int nid, double * p, int expnum)
{
	double sum=0.;
	//std::cout << expnum << std::endl;	
	for (int i=0;i<expnum;i++)
		sum +=pow(p[i],nid);

	return sum/expnum;
}
