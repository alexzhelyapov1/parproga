#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <mutex>
#include "mpi.h"


#define UNREACHABLE throw std::runtime_error("Something went wrong in function '" + std::string(__FUNCTION__) + "', line " + std::to_string(__LINE__) + ".");
#define K 100   // Количество шагов по времени
#define M 100   // Количество шагов по пространству
#define LEFT_VALUE 1
#define RIGHT_VALUE 2


double **ideal;
double **Ideal();

int rank, threads_number;
MPI_Status status;


// Функция f(t, x) = x + t
double f(double t, double x)
{
    return x + t;
}

// Функция φ(x) = cos(pi * x)
double phi(double x)
{
    return cos(M_PI*x);
}

// Функция ψ(t) = exp(-t)
double psi(double t)
{
    return exp(-t); 
}


double u_teor(double t, double x) {
    if (2 * t <= x) {
        return x * t - std::pow(t, 2) / 2 + cos(M_PI * (2 * t - x));
    }
    return x * t - std::pow(t, 2) / 2 + std::pow((2 * t - x), 2) / 8 + exp(-(t - x / 2));
}


void ThreadCompute(const int m_from, const int m_to,
                            const double a, const double h, const double tau, int rank)
{
    if (m_to < m_from) {
        UNREACHABLE
    }
    int i = 0;

    int m_range = m_to - m_from;
    double u_recvd_left = 0;
    double u_recvd_right = 0;

    // std::cout << "m_from = " << m_from << ", m_to = " << m_to << ", m_range = " << m_range << std::endl;

    // Выделение памяти для решения
    double** u = new double*[K + 1];
    for (int k = 0; k < K + 1; ++k) {
        u[k] = new double[m_range];
    }


    // Начальное условие
    for (int m = 0; m < m_range; ++m) {
        u[0][m] = phi((m + m_from) * h);
    }


    if (m_to != M + 1) {
        // UNREACHABLE
        // double *tmp = u + m_range - 1
        MPI_Send(u[0] + m_range - 1, 1, MPI_DOUBLE, rank + 1, LEFT_VALUE, MPI_COMM_WORLD);
    }


    // Если левая часть сетки
    if (m_from == 0) {
        // Граничное условие
        for (int k = 0; k <= K; ++k) {
            u[k][0] = psi(k * tau);
        }
    }
    else {
        MPI_Recv(&u_recvd_left, 1, MPI_DOUBLE, rank - 1, LEFT_VALUE, MPI_COMM_WORLD, &status);
        u[1][0] = u[0][0] - a * tau / h * (u[0][0] - u_recvd_left) + tau * f(0, m_from * h);
    }



    // Первый временной слой (k = 0) - схема левый уголок
    for (int m = 1; m < m_range; ++m) {
        u[1][m] = u[0][m] - a * tau / h * (u[0][m] - u[0][m - 1]) + tau * f(0, (m + m_from) * h);
    }


    // Остальные временные слои (k = 1, 2, ..., K-1) - схема крест
    for (int k = 1; k < K; ++k) {

        if (m_from != 0) {
            MPI_Send(u[k], 1, MPI_DOUBLE, rank - 1, RIGHT_VALUE, MPI_COMM_WORLD);
        }
        else {
            if (rank != 1) {
                UNREACHABLE // need
            }
            // Граничное условие
            // u[k][0] = psi((k) * tau);
            // u[k + 1][0] = psi((k + 1) * tau);
        }


        if (m_to != M + 1) {
            MPI_Send(u[k] + m_range - 1, 1, MPI_DOUBLE, rank + 1, LEFT_VALUE, MPI_COMM_WORLD);
        }



        for (int m = 1; m < m_range - 1; ++m) {
            u[k + 1][m] = u[k - 1][m] - a * tau / h * (u[k][m + 1] - u[k][m - 1]) + 2 * tau * f(k * tau, (m + m_from) * h);

            double u_t = u_teor((k + 1) * tau, (m + m_from) * h);
            double rel_err = (u[k + 1][m] - u_t) / u_t * 100;

            // if (rel_err > 10 && i < 5) {
            //     std::cout << k * tau << ", " << (m + m_from) * h << ", " << u[k + 1][m] << " " << u_t << std::endl;
            //     std::cout << "Rank = " << rank << ", " << k + 1 << " " << m << "\n---\n";
            //     i++;
            // }
        }



        if (m_to == M + 1) {
            if (rank != threads_number - 1) {
                std::cout << "Rank = " << rank << std::endl;
                UNREACHABLE
            }

            u[k + 1][m_range - 1] = u[k][m_range - 1] - a * tau / h * (u[k][m_range - 1] - u[k][m_range - 2]) + tau * f(k * tau, M * h);
        }
        else {
            MPI_Recv(&u_recvd_right, 1, MPI_DOUBLE, rank + 1, RIGHT_VALUE, MPI_COMM_WORLD, &status);
            u[k + 1][m_range - 1] = u[k - 1][m_range - 1] - a * tau / h * (u_recvd_right - u[k][m_range - 2]) + 2 * tau * f(k * tau, (m_to - 1) * h);
        }

        if (m_from != 0) {
            MPI_Recv(&u_recvd_left, 1, MPI_DOUBLE, rank - 1, LEFT_VALUE, MPI_COMM_WORLD, &status);

            u[k + 1][0] = u[k - 1][0] - a * tau / h * (u[k][1] - u_recvd_left) + 2 * tau * f(k * tau, m_from * h);
        }
    }

    for (int k = 0; k < K + 1; ++k) {
        for (int m = 0; m < m_range; ++m) {
            if (u[k][m] != ideal[k][m + m_from]) {
                UNREACHABLE
            }
        }
    }

    // std::cout << "All is correct!\n";

    double *buffer = new double[(K + 1) * m_range];

    for (int k = 0; k < K + 1; ++k) {
        memcpy(buffer + m_range * k, u[k], m_range * sizeof(double));
    }


    MPI_Send(buffer, (K + 1) * m_range, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);


    // Освобождение памяти
    for (int k = 0; k <= K; ++k) {
        delete[] u[k];
    }

    delete[] u;
    delete[] buffer;
}



int main(int argc, char **argv)
{
    // Параметры задачи
    double T = 1.0;        // Время
    double X = 1.0;        // Пространство
    double tau = T / K;    // Шаг по времени
    double h = X / M;      // Шаг по пространству
    double a = 1.0;        // Коэффициент переноса
    double** u;

    ideal = Ideal();


    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &threads_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0) {
         // Выделение памяти для решения
        u = new double*[K + 1];
        for (int k = 0; k <= K; ++k) {
            u[k] = new double[M + 1];
        }


        double *buffer = new double[(K + 1) * (M + 1)];
        int from_to[2];

        for (int dest_thread = 1; dest_thread < threads_number; ++dest_thread) {
            from_to[0] = (M + 1) / (threads_number - 1) * (dest_thread - 1);
            from_to[1] = (M + 1) / (threads_number - 1) * dest_thread;

            if (dest_thread == threads_number - 1) {
                from_to[1] = M + 1;
            }

            MPI_Send(from_to, 2, MPI_INT, dest_thread, 0, MPI_COMM_WORLD);
        }

        for (int dest_thread = 1; dest_thread < threads_number; ++dest_thread) {
            from_to[0] = (M + 1) / (threads_number - 1) * (dest_thread - 1);
            from_to[1] = (M + 1) / (threads_number - 1) * dest_thread;

            if (dest_thread == threads_number - 1) {
                from_to[1] = M + 1;
            }

            std::cout << "Waiting " << dest_thread << " ...";
            MPI_Recv(buffer, (K + 1) * (from_to[1] - from_to[0]), MPI_DOUBLE, dest_thread, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            std::cout << " got.\n";

            for (int k = 0; k < K + 1; ++k) {
                memcpy(u[k] + from_to[0], buffer + k * (from_to[1] - from_to[0]), (from_to[1] - from_to[0]) * sizeof(double));
            }
        }

        


        // Вывод результата в файл .csv чтобы потом визуалилировать
        std::ofstream outfile("/home/alex/parproga/gen/solution_new_new.csv");
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
        delete[] buffer;

    }
    else {
        int from_to[2];
        MPI_Recv(from_to, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        ThreadCompute(from_to[0], from_to[1], a, h, tau, rank);
    }

    
    MPI_Finalize();

    return 0;
}



double **Ideal() {
    // Параметры задачи
    double T = 1.0;    // Время
    double X = 1.0;    // Пространство
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

    return u;
}