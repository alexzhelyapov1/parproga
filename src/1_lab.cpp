#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include "mpi.h"
#include <mutex>


#define UNREACHABLE throw std::runtime_error("Something went wrong in function '" + std::string(__FUNCTION__) + "'");
#define K 100   // Количество шагов по времени
#define M 100   // Количество шагов по пространству


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


struct TransferData {
    double u[2][M + 1];
    int m_from;
    int m_to;
    int m_max;
    double a;
    double h;
    double tau;
    int k;
    bool end;
};


TransferData* PackData(double** u, const int m_from, const int m_to, const int m_max,
                    const double a, const double h, const double tau, const int k)
{
    TransferData* buffer = new TransferData;
    
    memcpy(buffer->u[0], u[k - 1], (M + 1) * sizeof(double));
    memcpy(buffer->u[1], u[k], (M + 1) * sizeof(double));

    buffer->m_from = m_from;
    buffer->m_to = m_to;
    buffer->m_max = m_max;
    buffer->a = a;
    buffer->h = h;
    buffer->tau = tau;
    buffer->k = k;
    buffer->end = false;

    return buffer;
}


double* ComputeTimeLayer(double u[2][M + 1], const int m_from, const int m_to, const int m_max,
                            const double a, const double h, const double tau, const int k)
{
    if (m_to < m_from || m_max < m_to || u == nullptr || m_from < 1) {
        UNREACHABLE
    }

    double* res = new double[m_to + 1 - m_from];

    // Временные слои (k = 1, 2, ..., K-1) - схема крест
    for (int m = m_from; m <= m_to; ++m) {

        if (m != m_max) {
            res[m - m_from] = u[0][m] - a * tau / h * (u[1][m + 1] - u[1][m - 1]) + 2 * tau * f(k * tau, m * h);
            
        }

        // Последний узел по пространству рассчитывается с помощью схемы "левый уголок"
        else {
            res[m - m_from] = u[1][m] - a * tau / h * (u[1][m] - u[1][m - 1]) + tau * f(k * tau, m * h);
        }

    }

    return res;
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

    
    int rank, threads_number;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &threads_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0) {
        // Выделение памяти для решения
        u = new double*[K + 1];
        for (int k = 0; k <= K; ++k) {
            u[k] = new double[M + 1];
        }

        double res[M + 1];

        // Начальное условие
        for (int m = 0; m <= M; ++m) {
            u[0][m] = phi(m * h);
        }

        // Граничное условие
        for (int k = 0; k <= K; ++k) {
            u[k][0] = psi(k * tau);
        }


        // Первый временной слой (k = 1) - схема левый уголок
        for (int m = 1; m <= M; ++m) {
            u[1][m] = u[0][m] - a * tau / h * (u[0][m] - u[0][m - 1]) + tau * f(0, m * h);
        }



        for (int k = 1; k < K; ++k) {
        // for (int k = 1; k < 2; ++k) {
            // TransferData* buffer = PackData(u, 1, M, M, a, h, tau, k);

            // Send data for calculation
            TransferData* buffer = PackData(u, 0, 0, M, a, h, tau, k);

            for (int dest_thread = 1; dest_thread < threads_number; ++dest_thread) {

                int m_from = 1 + M / (threads_number - 1) * (dest_thread - 1);
                int m_to = M / (threads_number - 1) * dest_thread;

                if (dest_thread == threads_number - 1) {
                    m_to = M;
                }

                buffer->m_from = m_from;
                buffer->m_to = m_to;
                // std::cout << m_from << " " << m_to << "\n";

                MPI_Send(buffer, sizeof(TransferData), MPI_CHAR, dest_thread, 0, MPI_COMM_WORLD);

            }

            delete buffer;

            // Receive results
            for (int dest_thread = 1; dest_thread < threads_number; ++dest_thread) {

                int m_from = 1 + M / (threads_number - 1) * (dest_thread - 1);
                int m_to = M / (threads_number - 1) * dest_thread;

                if (dest_thread == threads_number - 1) {
                    m_to = M;
                }

                MPI_Recv(res, m_to + 1 - m_from, MPI_DOUBLE, dest_thread, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                memcpy(u[k + 1] + m_from, res, sizeof(double) * (m_to + 1 - m_from));
            }
        }

        for (int dest_thread = 1; dest_thread < threads_number; ++dest_thread) {
            TransferData buffer;
            buffer.end = true;
            MPI_Send(&buffer, sizeof(TransferData), MPI_CHAR, dest_thread, 0, MPI_COMM_WORLD);
        }

        // Вывод результата в файл .csv чтобы потом визуалилировать
        std::ofstream outfile("/home/alex/parproga/first_task/gen/solution_new.csv");
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

    }

    else {
        TransferData buffer;
        buffer.end = false;

        while (buffer.end == false) {

            MPI_Recv(&buffer, sizeof(TransferData), MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (buffer.end == true) {
                break;
            }

            double* res = ComputeTimeLayer(buffer.u, buffer.m_from, buffer.m_to, buffer.m_max, buffer.a,
                                                                                    buffer.h, buffer.tau, buffer.k);

            MPI_Send(res, buffer.m_to + 1 - buffer.m_from, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();

    return 0;
}
