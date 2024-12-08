#include <iostream>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int total_threads = omp_get_num_threads();

        #pragma omp critical
        {
            std::cout << "Thread number: " << thread_id << ". Total threads: " << 
            total_threads << ". Message: Hello World!" << std::endl;
        }
    }

    return 0;
}