//=======================================================================
// The 4th moment of a distribution are computed using 
// Single-Instruction Multiple-Data (SIMD) parallelization.
// Each core has access to the entire dataset but they compute
// a different moment. For ease of validity check, I created
// my dataset by sampling from a normal distribution,
// moments of which are analytically known.
//
// Compile by: mpic++ -O2 simd_4thmoment.cpp -o <executable file name>
//
// Run the executable by: mpirun -n <number of cores> <executable file name>
//
//Author: Navid Afrasiabian, 2022
// ======================================================================
#include <iostream>
#include "mpistuff.h"
#include <random>

double myPow(double ,int);
double myMoment(double *,int,int);

int main(int argc, char** argv)
{
        double * p;
	double * pw;
	double sum, Gsum;
	double myMomentTime = 0., allMomentTime;
	MPI_Stuff mympi(argc,argv); //Initializing MPI class
	MPI_Status status;

        //double p[1000000];       // Data Array
        const int nrolls=1000000;  // number of experiments
	p = new double[nrolls];
        //=============================================================
        // Using the random library, the boss core draws numbers from
        // a normal distribution of mean value 5 and standard deviation
        // of 2 and stores them in p. Then different parts of p are
	// sent to the other cores using MPI_Send. Self-send was tested
	// and for the version of MPI used, it was blocking and did not
	// put the program in halt.
        // ============================================================
        if (mympi.nID == 0)
        {
                std::cout << "Number of data points: " << nrolls << std::endl;

                std::default_random_engine generator;
                std::normal_distribution<double> distribution(5.0,2.0);

                for (int i=0; i<nrolls; ++i)
                        p[i] = distribution(generator);
		
		//=====================================================
		//To account for indivisibility of the number of data
		//points by the number of processors, I divided the
		//entire dataset by the number of Worker cores. The Boss
		//would only participate in the summation process when
		//remainder is not zero. In that case, the Boss handles
		//the extra data points.
		//=====================================================
        	int sumSize = nrolls/(mympi.nprocs-1);
		
		for (int i = 1;i<mympi.nprocs;i++)
			MPI_Send(&p[sumSize*(i-1)],sumSize,MPI_DOUBLE,i,10,MPI_COMM_WORLD);
		int bosSize = nrolls%(mympi.nprocs-1);
		if (bosSize)
			sum = myMoment(&p[sumSize*mympi.nprocs],bosSize,4);
	}
	else{
		int sumSize = nrolls/(mympi.nprocs-1);
		pw = new double[sumSize];
		
		MPI_Recv(&pw[0], sumSize, MPI_DOUBLE,0, 10, MPI_COMM_WORLD, &status);
		double timeNow = MPI_Wtime();
		sum = myMoment(&pw[0],sumSize,4);
		myMomentTime += MPI_Wtime()-timeNow;
		//std::cout << "myMoment time: " << myMomentTime << std::endl;
	}

	MPI_Reduce(&sum,&Gsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&myMomentTime,&allMomentTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if (mympi.nID == 0)
	{
		Gsum /= nrolls;
		std::cout << "4th moment is " << Gsum << std::endl;
		std::cout << "Total myMoment time: " << allMomentTime << std::endl;
	}
	else{// pw is only allocated on the Worker cores so it is freed only on those.
		delete []pw;
	}
	delete []p;
	return 0;
}

double myPow(double base,int exponent)
{
	double bpow = 1.;
	for (int i=0;i<exponent;i++)
		bpow *= base;

	return bpow;
}

double myMoment(double *darr,int numofterms, int n)
{
	double dmoment = 0.;

	for (int i=0;i<numofterms;i++)
		dmoment += myPow(*(darr+i),n);

	return dmoment;
}
