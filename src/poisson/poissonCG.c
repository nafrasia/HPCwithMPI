//===========================================================
// Parallel version of the Conjugate Gradient solver for the
// Poisson problem. I parallelized this code using a 1D array
// of cores/processors.
//
// Compile:
// mpicc ./poissonCG.c -lm -o poissonCG.exe
//
// Run:
// mpirun -n <number of cores> ./poissonCG.exe
//
// Author: Navid Afrasiabian, 2022
//===========================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

//------------------------------
// Global parameters
//------------------------------
#define N  801
#define Tol  0.0001
#define maxitr 2*N
#define h 0.1

//-------------------------------------
// Function definition and declaration
//-------------------------------------
//MPI related data structures
struct mpi_vars {
  int nprocs;
  int nID;
}; 

struct mpi_vars mpi_start(int argc, char** argv)
{
  struct mpi_vars this_mpi;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &this_mpi.nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi.nID);

  return this_mpi;
}
//------------------------------------------------
// Initializing the matrix/geometry
//------------------------------------------------
double **matrix(int m, int n)
{   
  /* Note that you must allocate the array as one block in order to use */
  /* MPI derived data types on them.  Note also that calloc initializes */
  /* all entries to zero. */
  
  double **ptr = (double **)calloc(m, sizeof(double *));
  ptr[0]=(double *)calloc(m*n, sizeof(double));
  for(int i = 1; i < m ;i++)
    ptr[i]=ptr[i-1]+n;
  return (ptr);
}

//------------------------------------------------
// Assign the size of domain handled by each core
//------------------------------------------------
int compute_my_size(struct mpi_vars mympi)
{
  int remainder = (N - 1) % mympi.nprocs;
  int size = (N - 1 - remainder)/mympi.nprocs;
  if(mympi.nID < remainder)  // extra/remainder rows are managed by cores with coreID < remainder 
    size = size + 2;
  else
    size = size + 1;
  return size;
}

//-----------------------------------------------------------------------
// Initializing the solution and residual arrays
//-----------------------------------------------------------------------
void initialize(double ** r, double **x, int size, struct mpi_vars mympi)
{
  //Subtract off b from boundary conditions x = 1 on the boundary
  for(int i = 1; i < size; i++){
    r[i][1] -= 1;
    r[i][N-1] -= 1;
    x[i][0]=1;    // we won't use the x on the boundary but this will be useful for the output at the end
    x[i][N]=1;
  }
  if(mympi.nID == 0)
  {
  	for(int j = 1; j < N; j++){
    		r[1][j] -= 1;
   	 	x[0][j]=1;
  	}
	x[0][0]=x[0][N]=1;
  }
  
  if(mympi.nID == mympi.nprocs -1)
  {
  	for(int j = 1; j < N; j++){
    		r[size-1][j] -= 1;
    		x[size][j]=1;
  	}
  	x[size][0]=x[size][N]=1;    //fill in the corners
  }
}

//------------------------------------------------------
// Computes the matrix multiplication between matrix A in
// Ax = b equation and d. d is the direction vector for
// the conjugate gradient method. For Poisson problem, A
// is obtained from finite difference expansion of the
// derivative operators. The negative sign is to make
// A a positive definite matrix.
//------------------------------------------------------
void AD(double **Ad, double **d, int start, int finish)
{ // Compute Ad= A*d using the auxiliary layer around the outside that is all zero 
  // for d to account for boundary
  int i, j;
  for(i = start; i < finish; i++)
    for(j = 1; j < N; j++){
      Ad[i][j] = -(d[i+1][j] + d[i-1][j] + d[i][j+1] + d[i][j-1]-4.0*d[i][j]);
    }
  return;
}
//---------------------------------------------------------
// Outputs the solution on each core to a text file with
// the core id in the title
//---------------------------------------------------------
void output(double **new, int size, struct mpi_vars the_mpi)
{
  char str[20];
  FILE *fp;

  sprintf(str,"Solution%d.Txt",the_mpi.nID);
  fp = fopen(str,"w");
  if(the_mpi.nID == 0) {
    for(int j = 0; j < N + 1; j++)
      fprintf(fp,"%6.4f ",new[0][j]);
    fprintf(fp,"\n");
  }
  for(int i = 1; i < size; i++) {
    for(int j = 0; j < N + 1; j++)
      fprintf(fp,"%6.4f ",new[i][j]);
    fprintf(fp,"\n");
  }
  if(the_mpi.nID == the_mpi.nprocs - 1){
    for(int j = 0; j < N + 1; j++)
      fprintf(fp,"%6.4f ",new[size][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}

//---------------------------------------------------------------
// Generates a single solution file.
// A global view of the system is set up based on the parallelization
// geometry. We then iterate through the data on each core and
// place them in the right location based on their relative position
// to other core and system domain.
//---------------------------------------------------------------
void output_onefile(double **new, int ysize, struct mpi_vars mpi)
{
  int xsize = N;
  /* Create derived datatype for local interior grid (output grid) */
  MPI_Datatype grid;
  int start[2] = {1, 1};  // indices of interior "origin"
  int arrsize[2] = {xsize+1, ysize+1};  // full local array size
  int gridsize[2] = {xsize+1 - 2, ysize+1 - 2};  // size of interior

  MPI_Type_create_subarray(2, arrsize, gridsize,
			   start, MPI_ORDER_FORTRAN, MPI_DOUBLE, &grid);
  MPI_Type_commit(&grid);

  /* Create derived type for file view, how local interior fits into global system */
  MPI_Datatype view;
  int nnx = xsize+1-2, nny = ysize+1-2; 
  int startV[2] = { 0 , mpi.nID*nny };
  int arrsizeV[2] = { nnx, mpi.nprocs*nny };
  int gridsizeV[2] = { nnx, nny };
  
  MPI_Type_create_subarray(2, arrsizeV, gridsizeV,
			   startV, MPI_ORDER_FORTRAN, MPI_DOUBLE, &view);
  MPI_Type_commit(&view);

  /* MPI IO */
  MPI_File fp;

  MPI_File_open(MPI_COMM_WORLD, "1Doutput.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY,
		MPI_INFO_NULL, &fp);

  MPI_File_set_view(fp, 0, MPI_DOUBLE, view, "native", MPI_INFO_NULL);
  MPI_File_write_all(fp, &new[0][0], 1, grid, MPI_STATUS_IGNORE);
  MPI_File_close(&fp);
}

//-------------------------------
// Main code
//-------------------------------
int main(int argc, char** argv) 
{

  struct mpi_vars mympi = mpi_start(argc, argv);

  double delta, lambda,Glambda, olddelta, Gdelta, alpha;
  double Ttot, TtotG, TcalcG=0.,Tcalc=0.;
  int size;
  size = compute_my_size(mympi);
  
  /* create arrays and initialize to all zeros */
  /* Note: Although these are column vectors in the algorithm, they represent locations */
  /* in 2D space so it is easiest to store them as 2D arrays */
  double **x = matrix(size+1,N+1);
  double **r = matrix(size+1,N+1);
  double **d = matrix(size+1,N+1);
  double **u = matrix(size+1,N+1);
  
  MPI_Request req_send10, req_send20;
  MPI_Request req_recv10,req_recv20;

  //Work out A x(0)
  AD(r,x,1,size);
  
  initialize(r,x,size,mympi); // setups up boundary conditions
  
  //Fill in d
  for (int i=1; i<size; i++)
    for (int j=1; j<N; j++)
      d[i][j]=-r[i][j];
  
  /* Compute first delta */
  delta=0.0;
  for (int i=1; i<size; i++)
    for (int j=1; j<N; j++)
      delta+= r[i][j]*r[i][j];
  MPI_Allreduce(&delta, &Gdelta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  double starttime=MPI_Wtime();  
  
  /* Main loop */
  int itr;
  for (itr=0; itr< maxitr; itr++){

    /* send/receive at top of compute block */
    req_send10 = req_recv20 = MPI_REQUEST_NULL;
    if(mympi.nID < mympi.nprocs - 1){
      MPI_Isend(&d[size-1][1], N-1, MPI_DOUBLE, mympi.nID+1, 10,
		MPI_COMM_WORLD, &req_send10);
      MPI_Irecv(&d[size][1], N-1, MPI_DOUBLE, mympi.nID+1, 20,
		MPI_COMM_WORLD,&req_recv20);
    }
    /* send/receive at bottom of compute block */
    req_send20 = req_recv10 = MPI_REQUEST_NULL;
    if(mympi.nID > 0){
      MPI_Isend(&d[1][1], N-1, MPI_DOUBLE, mympi.nID-1, 20,
		MPI_COMM_WORLD, &req_send20);
      MPI_Irecv(&d[0][1], N-1, MPI_DOUBLE, mympi.nID-1, 10,
		MPI_COMM_WORLD, &req_recv10);
    }
    double tmptime=MPI_Wtime(); 
    //compute interior
    AD(u,d,2,size-1);
    Tcalc +=MPI_Wtime()-tmptime;
    //Ensure completion of all passing before ghost calculations
    if(mympi.nID<mympi.nprocs-1){
    	MPI_Wait(&req_recv20, MPI_STATUS_IGNORE);
    }
    tmptime=MPI_Wtime();
    //compute top
    AD(u,d,size-1,size);
    Tcalc +=MPI_Wtime()-tmptime;
    //Ensure completion of all passing
    if(mympi.nID > 0){
    	MPI_Wait(&req_recv10, MPI_STATUS_IGNORE);
    }
    tmptime=MPI_Wtime();
    //Compute bottom
    AD(u,d,1,2);
    Tcalc +=MPI_Wtime()-tmptime;
    
    tmptime=MPI_Wtime();

    lambda=0.0;
    for (int i=1; i<size; i++)
      for (int j=1;j<N; j++)
	lambda+= d[i][j]*u[i][j];
    
    Tcalc +=MPI_Wtime()-tmptime;
    
    MPI_Allreduce(&lambda,&Glambda,1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    
    tmptime=MPI_Wtime();
    
    Glambda=Gdelta/Glambda;
    //For debugging
    /* 
    if(mympi.nID == 0)
	printf("itr %d lambda= %f \n", itr,Glambda);    
    */
    for (int i=1; i<size; i++)
      for (int j=1;j<N; j++) {
	x[i][j] += Glambda*d[i][j];
	r[i][j] += Glambda*u[i][j];
      }
    
    olddelta=Gdelta;
    delta=0.0;
    for (int i=1; i<size; i++)
      for (int j=1; j<N; j++)
	delta+= r[i][j]*r[i][j];
    
    Tcalc +=MPI_Wtime()-tmptime;
    
    MPI_Allreduce(&delta, &Gdelta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    tmptime=MPI_Wtime();

    //Convergence check
    if (sqrt(Gdelta) < Tol)
      break;
    
    alpha=Gdelta/olddelta;
    for (int i=1; i<size; i++)
      for (int j=1;j<N; j++) {
	d[i][j] = -r[i][j]+alpha*d[i][j];
      } 
    Tcalc+=MPI_Wtime()-tmptime;
  }
  Ttot=MPI_Wtime()-starttime;
  
  //Sum up all the calculation and passing times on every core and
  //store them on BOSS core
  MPI_Reduce(&Ttot, &TtotG, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
  MPI_Reduce(&Tcalc, &TcalcG, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
  
  //Performance output
  if(mympi.nID == 0){
	printf("N = %d\n", N);
  	printf("itr = %i delta = %g\n",itr,Gdelta);
	printf("Complete in %d iterations \n",itr);
	printf("Ttot = %f\n",TtotG/mympi.nprocs);
	printf("Tcomm = %f\n", (TtotG-TcalcG)/mympi.nprocs);
	printf("Tcalc = %f\n",TcalcG/mympi.nprocs);
  }
  
  /* Output result */
  //output(x,size,mympi);
  output_onefile(x, size, mympi);
  
  //Releasing memory
  free(r[0]);
  free(r);
  free(d[0]);
  free(d);
  free(u[0]);
  free(u);
  free(x[0]);
  free(x);
  MPI_Finalize();
  return 0;
}
