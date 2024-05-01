#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>

std::mutex cout_mutex;

inline double f(double x)
{
    return sin(1 / x);
}

inline double f_deriv(double x)
{
    return - cos(1 / x) / (x * x);
}


double Integral(double from, double to)
{
    assert(from <= to);

    cout_mutex.lock();
    std::cout << "Computing integral from, to: " << from << " " << to << std::endl;
    cout_mutex.unlock();

    unsigned long number_of_points = 1000000;
    double step = (to - from) / number_of_points;
    double sum = 0.0;

    for (int i = 0; i < number_of_points; ++i) {
        sum += f(from + i * step);
    }

    cout_mutex.lock();
    std::cout << "Computed: " << sum * step << std::endl;
    cout_mutex.unlock();

    return sum * step;
}


double LaunchParallelIntegral(int threads_number)
{
    double from = 1.1, to = 5.5;
    double from_local, to_local;

    if (from > to) {
        throw std::runtime_error("from > to");
    }


    std::vector<std::thread> threads;

    for (int i = 0; i < threads_number; i++) {

        from_local = from + (to - from) / threads_number * i;
        to_local = from + (to - from) / threads_number * (i + 1);

        if (i == threads_number - 1) {
            to_local = to;
        }

        threads.push_back(std::thread(Integral, from_local, to_local));
    }

    for (auto& thread : threads) {
        thread.join();
    }


    return 0;
}


int main(int argc, char **argv)
{	
    int threads_number = 4;
    LaunchParallelIntegral(threads_number);

    return 0;
}
