//=====================================================
// Solving the Poisson partial differential equation
// using Conjugate Gradient method.
// This code was originally provided as the basis by Dr
// Denniston.
// 
// Compile:
// gcc -O3 serialCG.c -lm -o serialCG.exe
// 
// Run:
// ./serialCG.exe
//=====================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N  2001
#define Tol  0.0001
#define maxitr 2*N
#define h 0.1

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


void initialize(double ** r, double **x, int size)
{
  /* Subtract off b from a particularly boring set of boundary conditions where x=1 on the boundary */
  for(int i = 1; i < size; i++){
    r[i][1] -= 1;
    r[i][size-1] -= 1;
    x[i][0]=1;    // we won't use the x on the boundary but this will be useful for the output at the end
    x[i][size]=1;
  }
  for(int j = 1; j < size; j++){
    r[1][j] -= 1;
    r[size-1][j] -= 1;
    x[0][j]=1;
    x[size][j]=1;
  }
  x[0][0]=x[0][size]=x[size][0]=x[size][size]=1;    /* fill in the corners */
}


void AD(double **Ad, double **d, int start, int finish)
{ // Compute Ad= A*d using the auxiliary layer around the outside that is all zero 
  // for d to account for boundary */
  int i, j;
  for(i = start; i < finish; i++)
    for(j = 1; j < N; j++){
      Ad[i][j] = -(d[i+1][j] + d[i-1][j] + d[i][j+1] + d[i][j-1]-4.0*d[i][j]);
    }
  return;
}


void output(double **x, int size)
{
  char str[20];
  FILE *fp;

  sprintf(str,"Solution.Txt");
  fp = fopen(str,"wt");

  for(int i = 0; i < size+1; i++){
    fprintf(fp,"\n");
    for(int j = 0; j < size + 1; j++)
      fprintf(fp,"%6.4f ",x[size-i][j]);
  }
  fclose(fp);
}


int main(int argc, char** argv) 
{   
  double delta, lambda, olddelta, alpha;
  clock_t starttime, endtime;
  /* create arrays and initialize to all zeros */
  /* Note: Although these are column vectors in the algorithm, they represent locations */
  /* in 2D space so it is easiest to store them as 2D arrays */
  double **x = matrix(N+1,N+1);
  double **r = matrix(N+1,N+1);
  double **d = matrix(N+1,N+1);
  double **u = matrix(N+1,N+1);
  /* Work out A x(0), note that we make use of the x being zero on the boundary so that */
  /* we don't need the boundaries to be special cases in the iteration equations */ 
  AD(r,x,1,N);
  
  initialize(r,x,N); // setups up boundary conditions
  
  /* Fill in d */
  for (int i=1; i<N; i++)
    for (int j=1; j<N; j++)
      d[i][j]=-r[i][j];
  
  /* Compute first delta */
  delta=0.0;
  for (int i=1; i<N; i++)
    for (int j=1; j<N; j++)
      delta+= r[i][j]*r[i][j];
  printf("itr = %i delta = %g\r",0,delta);
  starttime = clock();
  /* Main loop */
  int itr;
  for (itr=0; itr< maxitr; itr++){
    AD(u,d,1,N);
    
    lambda=0.0;
    for (int i=1; i<N; i++)
      for (int j=1;j<N; j++)
	lambda+= d[i][j]*u[i][j];
    lambda=delta/lambda;
    //for debugging
    //printf("iter %d lambda= %f\n",itr, lambda);   
    for (int i=1; i<N; i++)
      for (int j=1;j<N; j++) {
	x[i][j] += lambda*d[i][j];
	r[i][j] += lambda*u[i][j];
      }
    
    olddelta=delta;
    delta=0.0;
    for (int i=1; i<N; i++)
      for (int j=1; j<N; j++)
	delta+= r[i][j]*r[i][j];
    printf("itr = %i delta = %g\r",itr,delta);
    
    if (sqrt(delta) < Tol)
      break;
    
    alpha=delta/olddelta;
    for (int i=1; i<N; i++)
      for (int j=1;j<N; j++) {
	d[i][j] = -r[i][j]+alpha*d[i][j];
      } 
  }
  endtime = clock();
  double Ttot = ((double)(endtime-starttime))/CLOCKS_PER_SEC;
  printf("\nComplete in %d iterations \n",itr);
  printf("Ttot = %f\n", Ttot);
  /* Output result */
  output(x,N);
  
  free(r[0]);
  free(r);
  free(d[0]);
  free(d);
  free(u[0]);
  free(u);
  free(x[0]);
  free(x);
  return 0;
}
