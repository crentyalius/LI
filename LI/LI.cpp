// LI.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

int main()
{
    std::cout << "Hello World!\n";
}

double rss(int N, double* Y, double* X, double A, double B) {
    double sum_rss = 0.0;

    for (int i = 0; i < N; ++i) {
        double Y_predicted = A + B * X[i];
        double deviation = std::abs(Y[i] - Y_predicted);
        sum_rss += deviation;
    }

    return sum_rss;
}
