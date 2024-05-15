#include <stdio.h>
#include "mpi.h"
#include <unistd.h>
#include <iostream>

#define STEPS 10000000


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

    return sum;
}


int main(int argc, char **argv)
{	
    int rank, threads_number;
    MPI_Status status;
    double sum = 0.0;
    double sum_buf = 0.0;
    double step = 1.0 / static_cast<double>(STEPS);
    unsigned long *borders = new unsigned long[2];

    double cpu_time_start, cpu_time_fini;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &threads_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0) {
        cpu_time_start = MPI_Wtime();

        for (int dest_thread = 1; dest_thread < threads_number; dest_thread++) {
            
            borders[0] = STEPS / (threads_number - 1) * (dest_thread - 1);
            borders[1] = STEPS / (threads_number - 1) * dest_thread - 1;

            if (dest_thread == threads_number - 1) {
                borders[1] = STEPS;
            }

            MPI_Send(borders, 2, MPI_UNSIGNED_LONG, dest_thread, 0, MPI_COMM_WORLD);
        }

        for (int dest_thread = 1; dest_thread < threads_number; dest_thread++) {
            MPI_Recv(&sum_buf, 1, MPI_DOUBLE, dest_thread, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sum += sum_buf;
        }

        double pi = sum * step;

        cpu_time_fini = MPI_Wtime();

        std::cout << "PI: " << pi << std::endl;
        std::cout << "TIME: " << (cpu_time_fini - cpu_time_start) * 1000 << " ms." << std::endl;
    }
    else {
        MPI_Recv(borders, 2, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        sum_buf = CountPi(borders[0], borders[1], step);
        MPI_Send(&sum_buf, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    
    delete[] borders;
    
    return 0;
}

// Best time at 16 nodes. (number of processors in my comp)