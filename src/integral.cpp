#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>


#define UNREACHABLE throw std::runtime_error("Something went wrong in function '" + std::string(__FUNCTION__) + "', line " + std::to_string(__LINE__) + ".");

using namespace std;

std::mutex cout_mutex, result_mutext;
double result = 0.0;
// or
// std::atomic<double> result;


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


int countPeaks(double a, double b) {
    int peaks = 0;
    for (double x = a + 0.01; x < b; x += 0.01) {
        if (f(x) * f(x - 0.01) < 0) {
            peaks++;
        }
    }
    return peaks;
}


// Функция для разбивки интервала на N интервалов с равномерным распределением пиков
vector<double> createIntervals(double a, double b, int N) {
    vector<double> intervals;
    intervals.push_back(a);
    std::cout << "Helio\n";

    // Определяем количество пиков на всем интервале
    int totalPeaks = countPeaks(a, b);


    if (totalPeaks < N) {
        double step = (b - a) / N;
        for (int i = 1; i < N; i++) {
            intervals.push_back(a + i * step);
        }
    }
    else {
        double currentPosition = a;
        double lastPosition = a;
        int currentPeaks = 0;
        int peaksPerInterval = totalPeaks / N;

        for (int i = 1; i < N; i++) {
        // Находим следующую точку разбиения, пока не достигнем желаемого количества пиков
            while (currentPeaks < i * peaksPerInterval && currentPosition < b) {
                currentPosition += 0.01;
                currentPeaks += countPeaks(currentPosition - 0.01, currentPosition);
                std::cout << "-\n";
            }

            if (currentPosition > b) {
                std::cout << "Ended at " << lastPosition << std::endl;

                int restIntervals = N - i;
                for (int k = i; k < N; k++) {
                    intervals.push_back(lastPosition + (b - a) / restIntervals * k);
                }
                break;
            }

            intervals.push_back(currentPosition);
            lastPosition = currentPosition;
        }

    }

    // // Если пиков меньше, чем интервалов, распределяем интервалы равномерно
    // if (totalPeaks < N) {
    //     double step = (b - a) / N;
    //     for (int i = 1; i < N; i++) {
    //         intervals.push_back(a + i * step);
    //     }
    // } else { // Если пиков больше или равно интервалам, действуем как раньше
    //     // Вычисляем желаемое количество пиков на каждом интервале
    //     int peaksPerInterval = totalPeaks / N;

    //     // Делим интервал на N частей, учитывая количество пиков
    //     double currentPosition = a;
    //     int currentPeaks = 0;
    //     for (int i = 1; i < N; i++) {
    //     // Находим следующую точку разбиения, пока не достигнем желаемого количества пиков
    //         while (currentPeaks < i * peaksPerInterval && currentPosition < b) {
    //             currentPosition += 0.01;
    //             currentPeaks += countPeaks(currentPosition - 0.01, currentPosition);
    //             std::cout << "-\n";
    //         }
    //         intervals.push_back(currentPosition);
    //     }
    // }

  intervals.push_back(b);

  return intervals;
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
    double from = 0.01, to = 10;
    double from_local, to_local;

    if (from > to || from <= 0) {
        UNREACHABLE
    }

    // int first_pic = (1 / from + 3 * M_PI / 2) / (2 * M_PI);
    // int last_pic = (1 / to + 3 * M_PI / 2) / (2 * M_PI);
    // // int number_of_pics = last_pic - first_pic;

    // std::cout << "Pics number: " << first_pic << " " << last_pic << std::endl;

    // vector<double> intervals = createIntervals(from, to, threads_number);

    // for (auto i : intervals) {
    //     std::cout << i << " ";
    // }
    // std::cout << "\n";

    

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
    int threads_number = 8;
    LaunchParallelIntegral(threads_number);

    return 0;
}
