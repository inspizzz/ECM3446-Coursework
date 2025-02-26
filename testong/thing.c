/******************************************************************
Description: The second loop example from ECM3446/ECMM461 unit 2.5

Purpose: demonstrate the use of an OpenMP for parallelising a for
         loop. The loop calculates the function y=log(x) and is 
         parallelised using a combined pragma for starting the 
         parallel region and parallelising the for loop.

Notes: compile with the appropriate OpenMP flag for
       your compiler e.g. gcc -fopenmp -o hello hello.c

       This is the example used in the first demo in unit 2.5.


******************************************************************/


#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(){

  const int N=12;
  float x[N],y[N];

#pragma omp parallel for
  for (int i=0; i<N; i++){
    printf("Thread %d calculating log(%d)\n", omp_get_thread_num(), i+1);
    x[i]= (float) (i+1);
    y[i]= log(x[i]);
  }
  
  return 0;
}