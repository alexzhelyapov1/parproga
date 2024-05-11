#include <iostream>
#include <fstream>
#include <cmath>


// Функция f(t, x) = x + t
double f(double t, double x) {
    return x + t;
}

// Функция φ(x) = cos(pi * x)
double phi(double x) {
    return cos(M_PI*x);
}

// Функция ψ(t) = exp(-t)
double psi(double t) {
    return exp(-t); 
}

int main() {
    // Параметры задачи
    double T = 1.0;    // Время
    double X = 1.0;    // Пространство
    int K = 100;        // Количество шагов по времени
    int M = 100;        // Количество шагов по пространству
    double tau = T / K;    // Шаг по времени
    double h = X / M;        // Шаг по пространству
    double a = 1.0;            // Коэффициент переноса

    // Выделение памяти для решения
    double** u = new double*[K + 1];
    for (int k = 0; k <= K; ++k) {
        u[k] = new double[M + 1];
    }

    // Начальное условие
    for (int m = 0; m <= M; ++m) {
        u[0][m] = phi(m * h);
    }

    // Граничное условие
    for (int k = 0; k <= K; ++k) {
        u[k][0] = psi(k * tau);
    }


    // Первый временной слой (k = 0) - схема левый уголок
    for (int m = 1; m <= M; ++m) {
        u[1][m] = u[0][m] - a * tau / h * (u[0][m] - u[0][m - 1]) + tau * f(0, m * h);
    }

    // Остальные временные слои (k = 1, 2, ..., K-1) - схема крест
    for (int k = 1; k < K; ++k) {

        for (int m = 1; m <= M; ++m) {
            u[k + 1][m] = u[k - 1][m] - a * tau / h * (u[k][m + 1] - u[k][m - 1]) + 2 * tau * f(k * tau, m * h);
        }

        // Последний узел по пространству рассчитывается с помощью схемы "левый уголок"
        u[k + 1][M] = u[k][M] - a * tau / h * (u[k][M] - u[k][M - 1]) + tau * f(k * tau, M * h);
    }


    // Вывод результата в файл .csv чтобы потом визуалилировать
    std::ofstream outfile("solution.csv");
    outfile << "t,x,u" << std::endl;
    for (int k = 0; k <= K; ++k) {
        for (int m = 0; m <= M; ++m) {
            outfile << k * tau << "," << m * h << "," << u[k][m] << std::endl;
        }
    }
    outfile.close();


    // Освобождение памяти
    for (int k = 0; k <= K; ++k) {
        delete[] u[k];
    }
    delete[] u;

    return 0;
}