//===================================================================
// Synchronization modes are tested in the series of codes named
// synctest*.cpp. In synctest2.cpp, blocking communication is used.
// A time consuming task is assigned to one of the cores to
// demonstrate the impact of blocking communications. The function of
// MPI_Barrier() is demonstrated.
//
// Compile: mpic++ -O2 synctest2.cpp -o synctest2.exe
//
// Run: mpirun -n <number of cores> synctest2.exe
//
// Author: Navid Afrasiabian, 2022
//===================================================================
#include <iostream>
#include "mpistuff.h"

int findPrime(int);

int main(int argc, char **argv)
{
	MPI_Stuff mympi(argc,argv); //Declaring an MPI_Stuff class.
                                    //This initializes MPI
	MPI_Status status;

	int numoftasks = 10, task,core,ctask;
	
	for(task=1;task<=numoftasks;task++){
		core = task % mympi.nprocs;
		if(core == mympi.nID && mympi.nID!=0)
		{
			//----------------------------------------
                        // Core 2 performs computations. All cores
                        // see the instructions. The if condition
                        // ensures that the core with coreID = 2
                        // runs findPrime() function
                        // ---------------------------------------
			if(core == 2)
				findPrime(10000);
			MPI_Send(&task,1,MPI_INT,0,10,MPI_COMM_WORLD);
		}
		//---------------------------------------------------------
                // The Boss core is responsible for printing results and
                // performs no computations. For this, the boss core
                // has to receive the data from the worker cores.
                // Since blocking communications are used, if a worker
                // sends data while the boss core is busy, the worker core
                // would stay on hold till the data is sent resulting in a
                // delay in the performance of that core.
                //---------------------------------------------------------
		if (mympi.nID == 0){
			if (core == 0)
				printf("Core %d will do %d\n", mympi.nID,task);
			else{
				MPI_Recv(&ctask,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
				printf("Core %d will do %d\n", status.MPI_SOURCE,ctask);
			}
		}
		
	}
	//-------------------------------------------------
	// MPI_Barrier() makes sure that all communications
	// are complete to this point before rest of 
	// the instructions are passed onto the cores.
	// This can be used when we want to synchronize all
	// cores and have the most up-to-date data. This is
	// useful/crucial when computation is performed on 
	// the updated data
	//-------------------------------------------------
	MPI_Barrier(MPI_COMM_WORLD);
	for(task=11;task<=10+numoftasks;task++){
		core = task % mympi.nprocs;
		if(core == mympi.nID && mympi.nID!=0)
		{
			if(core == 2)
				findPrime(10000);
			MPI_Send(&task,1,MPI_INT,0,10,MPI_COMM_WORLD);
		}

		if (mympi.nID == 0){
			if (core == 0)
				printf("Core %d will do %d\n", mympi.nID,task);
			else{
				MPI_Recv(&ctask,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&status);
				printf("Core %d will do %d\n", status.MPI_SOURCE,ctask);
			}
		}
		
	}
	return 0;
}

//=========================================================
// The findPrime function finds the number of prime numbers
// up to n. This function can be improved by ignoring even
// numbers and dividing only up to sqrt(n) but since the 
// purpose of this function is to waste time, those fixes
// are not implemented!
// ========================================================
int findPrime(int n)
{
	int num = 2, count=0, divisor,r=1;
	while(num<n)
	{
		divisor=2;
		while(divisor<num)
		{
			r = num % divisor;
			if (r==0)
				break;
			
			divisor++;

		}
		if (r!=0){
			//std::cout << num << std::endl;
			count++;
		}
		num++;
	}
	return count;
}
