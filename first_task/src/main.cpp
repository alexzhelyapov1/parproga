#include <stdio.h>
#include "mpi.h"
#include <unistd.h>
#include <iostream>

#define STEPS 10000


double CountPi(unsigned long left, unsigned long right, double step)
{
    double x = 0;
    double sum = 0;

    for (unsigned long i = left; i < right; ++i)
    {
        x = ( i + .5 ) * step;
        sum += 4.0 / (1. + x * x);
    }

    double pi = sum * step;
    std::cout << "I count it: " << pi << std::endl;

    return sum;
}





int main(int argc, char **argv)
{	
    int rank, size;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        for (int to_thread = 1; to_thread < thread_size; to_thread++) {
            MPI_Send(&interval, 2, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
        }
    }
    else {
      
    }
    
    
    MPI_Finalize();
    
    printf("Process: %d, size: %d\n", rank, size);
    
    return 0;
}