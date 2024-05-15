#include <stdio.h>
#include "mpi.h"
#include <unistd.h>
#include <iostream>


int main(int argc, char **argv)
{	
    int rank, threads_number;
    MPI_Status status;

    double cpu_time_start, cpu_time_fini;
    int transactions_number = 10;

    int some_data = 5;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &threads_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (threads_number < 2) {
        throw std::runtime_error("Not enough number of threads.");
    }

    if (rank == 0) {

        cpu_time_start = MPI_Wtime();

        for (int i = 0; i < transactions_number; ++i) {
            MPI_Send(&some_data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&some_data, 1, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }

        cpu_time_fini = MPI_Wtime();
        
        std::cout << "One transaction time: " << (cpu_time_fini - cpu_time_start) * 1000 / transactions_number << " ms." << std::endl;

    }
    else if (rank == 1) {

        for (int i = 0; i < transactions_number; ++i) {

            MPI_Recv(&some_data, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Send(&some_data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        }

    }
    
    MPI_Finalize();
    
    return 0;
}
