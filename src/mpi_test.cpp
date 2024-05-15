#include <iostream>
#include "mpi.h"


void CheckUAddress(const int rank, const int threads_number, const double* u) {
    MPI_Status status;
    int result = 0;

    if (rank == 0) {
        unsigned long long test = 0x55bd26bc34f0;

        for (int dest_thread = 1; dest_thread < threads_number; dest_thread++) {
            MPI_Send(&u, 1, MPI_UNSIGNED_LONG_LONG, dest_thread, 0, MPI_COMM_WORLD);
        }

        for (int dest_thread = 1; dest_thread < threads_number; dest_thread++) {
            MPI_Recv(&result, 1, MPI_INT, dest_thread, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (result != 0) {
                std::cout << "Error! I am 0 thread, thread number " << dest_thread << " failed with wrong U address!\n";

            }
            else {
                std::cout << "Thread number " << dest_thread << ": OK.\n";
            }
        }
    }
    else {
        unsigned long long address_buffer = 0;


        MPI_Recv(&address_buffer, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (address_buffer != reinterpret_cast<unsigned long long>(u)) {
            result = -1;
        }

        std::cout << rank << ". u = " << u << ", u = " << reinterpret_cast<unsigned long long>(u) << "\n";
        std::cout << rank << ". a = " << reinterpret_cast<double *>(address_buffer) << ", a = " << address_buffer << "\n";

        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}



int main(int argc, char **argv)
{
    // int rank, threads_number;

    // MPI_Init(&argc, &argv);
    
    // MPI_Comm_size(MPI_COMM_WORLD, &threads_number);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // double *a = new double[5];
    // for (auto i = 0; i < 5; i++) {
    //     a[i] = i;
    // }

    // a[4] = rank;

    // std::cout << a << "\n";
    // CheckUAddress(rank, threads_number, a);

    // // std::cout << "Threads number: " << threads_number << std::endl;
    // // std::cout << "Rank: " << rank << std::endl;
    // // for (auto i = 0; i < 5; i++) {
    // //     std::cout << a[i] << " ";
    // // }
    // // std::cout << "\n";
    // // std::cout << a << "\n---\n";

    // MPI_Finalize();

}