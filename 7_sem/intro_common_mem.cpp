#include <iostream>
#include <omp.h>
#include <mutex>
#include <condition_variable>

int main() {
    int shared_variable = 0;
    std::mutex mtx;
    std::condition_variable cv;
    int current_thread = 0;


    #pragma omp parallel shared(shared_variable, mtx, cv, current_thread)
    {
        int thread_id = omp_get_thread_num();

        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&](){ return thread_id == current_thread; });

        std::cout << "Thread: " << thread_id << ", shared_variable before action: "
            << shared_variable++ << ", shared_variable after action: " << shared_variable << std::endl;

        current_thread++;
        cv.notify_all();
    }

    return 0;
}