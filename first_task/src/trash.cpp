#include <iostream>

#define STEPS 10000


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
    std::cout << "I count it: " << pi << std::endl;

    return sum;
}


int main() 
{
    double step = 1./static_cast<double>(STEPS);
    double sum = 0;

    for (int i = 0; i < 4; ++i) {
        unsigned long left = STEPS / 4 * i;
        unsigned long right = STEPS / 4 * (i + 1);
        std::cout << left << " " << right << std::endl;
        sum += CountPi(left, right, step);
    }

    double pi = sum * step;
    std::cout << pi << std::endl;
}