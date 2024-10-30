//================================================================
// Synchronization modes are tested in the series of codes named
// synctest*.cpp. A time consuming task is assigned to one of the 
// cores to demonstrate the impact of blocking/non-blocking 
// communications.The impact of non-blocking communication is 
// tested in synctest5.cpp.
//
// Compile: mpic++ -O2 synctest5.cpp -o synctest5.exe
//
// Run: mpirun -n <number of cores> synctest5.exe
//
// Author: Navid Afrasiabian, 2022
//================================================================

#include <iostream>
#include "mpistuff.h"

int findPrime(int);

int main(int argc, char **argv)
{
	MPI_Stuff mympi(argc,argv);
	MPI_Status status;
	//------------------------------------------------------------------------
        // MPI_Request represents a handle on a non-blocking operation.
        // This is used by wait (MPI_Wait, MPI_Waitall, MPI_Waitany, MPI_Waitsome)
        // and test (MPI_Test, MPI_Testall, MPI_Testany, MPI_Testsome)
        // to know when the non-blocking operation handled completes.
        //------------------------------------------------------------------------
	MPI_Request reqr,reqs;
	reqs = MPI_REQUEST_NULL;
	reqr = MPI_REQUEST_NULL;

	int numoftasks = 10, task,core,ctask;
	
	for(task=1;task<=numoftasks;task++){
		core = task % mympi.nprocs;
		if(core == mympi.nID)
		{
			//----------------------------------------
                        // Core 2 performs computations. All cores
                        // see the instructions. The if condition
                        // ensures that the core with coreID = 2
                        // runs findPrime() function
                        // ---------------------------------------
			if(core == 2)
				findPrime(10000);
			MPI_Isend(&task,1,MPI_INT,0,10,MPI_COMM_WORLD,&reqs);
		}
	}
		//---------------------------------------------------
		// A separate loop is defined for the Boss core.
		// This combined with the non-blocking communications
		// allows the Boss core to carry on its tasks and
		// it does not have to wait for all the other cores
		// to finish their tasks.
		//---------------------------------------------------
		if (mympi.nID == 0){
			for(int i=0;i<numoftasks;i++){
				MPI_Irecv(&ctask,1,MPI_INT,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,&reqr);
				MPI_Wait(&reqr,&status);
				printf("Core %d will do %d\n", status.MPI_SOURCE,ctask);
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