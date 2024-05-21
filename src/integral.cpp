#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>


#define UNREACHABLE throw std::runtime_error("Something went wrong in function '" + std::string(__FUNCTION__) + "', line " + std::to_string(__LINE__) + ".");


std::mutex cout_mutex, result_mutext;
double result = 0.0;


// I made this function only for training with lock_guard instead mutex.
void add_to_result(double sum)
{
    std::lock_guard<std::mutex> guard(result_mutext);
    result += sum;
}


inline double f(double x)
{
    return sin(1 / x);
    // return sin(1 / x);
}


inline double f_deriv(double x)
{
    return - cos(1 / x) / (x * x);
}


void Integral(double from, double to)
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

    add_to_result(sum * step);
}


double LaunchParallelIntegral(int threads_number)
{
    double from = 0.1, to = 10;
    double from_local, to_local;

    if (from > to || from <= 0) {
        UNREACHABLE
    }

    // int first_pic = (1 / from + 3 * M_PI / 2) / (2 * M_PI);
    // int last_pic = (1 / to + 3 * M_PI / 2) / (2 * M_PI);
    // // int number_of_pics = last_pic - first_pic;

    // std::cout << "Pics number: " << first_pic << " " << last_pic << std::endl;

    

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

    std::cout << "Result: " << result << std::endl;

    return 0;
}


int main(int argc, char **argv)
{	
    int threads_number = 1;
    LaunchParallelIntegral(threads_number);

    return 0;
}
