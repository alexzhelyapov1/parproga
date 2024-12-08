#include <iostream>
#include <omp.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Usage: " + std::string(argv[0]) + " N");
    }

    int N = std::atoi(argv[1]);

    // Due to the peculiarities of adding double values, it is better to add large numbers with large ones, 
    // and small ones with small ones

    double sum = 0.0;
    #pragma omp parallel reduction(+:sum)
    {
        int thread_id = omp_get_thread_num();
        int total_threads = omp_get_num_threads();

        int start = thread_id * N / total_threads + 1;
        int end = (thread_id + 1) * N / total_threads + 1;

        if (start <= 0 || end <= 0) {
            throw std::runtime_error("Bad start or end!");
        }

        for (int i = start; i < end; ++i) {
            sum += 1.0 / i; 
        }
    }

    std::cout << "Sum: " << sum << std::endl;

    double control_sum = 0.0;
    for (int i = 1; i < N + 1; ++i) {
        control_sum += 1.0 / i;
    }

    std::cout << "Control sum: " << control_sum << std::endl;

    return 0;
}