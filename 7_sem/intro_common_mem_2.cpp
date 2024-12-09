#include <iostream>
#include <omp.h>

int main() {
    int shared_variable = 0;
    int total_threads = omp_get_max_threads();

    #pragma omp parallel for ordered num_threads(total_threads) shared(shared_variable) schedule(static, 1)
    for (int i = 0; i < total_threads; ++i) {
        int thread_id = omp_get_thread_num();

        #pragma omp ordered
        {
            std::cout << "Thread: " << thread_id << ", shared_variable before action: "
                << shared_variable++ << ", shared_variable after action: " << shared_variable << std::endl;
        }
    }

    return 0;
}
