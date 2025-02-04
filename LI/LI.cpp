#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <array>
#include <chrono>
#include <numeric>
using namespace std;


double mean(const std::vector<double>& X) {
    double sum = 0.0;
    for (double x : X) {
        sum += x;
    }
    return sum / X.size();
}


template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<T>& weights)

{
    if (values.empty() || weights.empty()) {
        throw std::invalid_argument("Values and weights must be non-empty");
    }

    if (values.size() != weights.size()) {
        throw std::invalid_argument("Sizes of values and weights must match");
    }

    // Преобразуем веса в кумулятивные суммы
    std::vector<size_t> cumulativeWeights(weights.size());
    std::partial_sum(weights.begin(), weights.end(), cumulativeWeights.begin());

    size_t totalWeight = cumulativeWeights.back();

    // Найти медианный индекс
    size_t medianIndex;
    if (totalWeight % 2 == 0) {
        // Четное число элементов
        medianIndex = totalWeight / 2 - 1;
    }
    else {
        // Нечетное число элементов
        medianIndex = totalWeight / 2;
    }

    // Бинарный поиск медианного значения
    auto it = std::upper_bound(cumulativeWeights.begin(), cumulativeWeights.end(), medianIndex);
    size_t rank = std::distance(cumulativeWeights.begin(), it);

    // Определить медианное значение
    T medianValue = values[rank];

    return medianValue;
}

template<typename T> T median(const std::vector<T>& values)
{
    return weightedMedian(values, values);
}


void generateMatrix(std::vector<std::vector<double>>& matrix, int rows, int cols) {
    static std::default_random_engine generator;
    static std::uniform_real_distribution<double> distribution(-10.0, 10.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = distribution(generator);
        }
    }
}

std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> inverse(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[i][j] = 0.0;
        }
        inverse[i][i] = 1.0;
    }
    for (int k = 0; k < n; ++k) {
        for (int i = k; i < n; ++i) {
            double pivot = matrix[i][k] / matrix[k][k];
            for (int j = k; j < n; ++j) {
                inverse[i][j] -= pivot * matrix[k][j];
            }
            inverse[i][k] = pivot;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                inverse[i][j] /= matrix[i][i];
            }
            else if (i < j) {
                inverse[i][j] = -(inverse[j][i]);
            }
            else {
                inverse[i][j] = inverse[j][i];
            }
        }
    }
    return inverse;
}

void leastSquares(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b, std::vector<double>& a_l, std::vector<double>& b_l) {
    int n = a.size();
    std::vector<std::vector<double>> X(n, std::vector<double>(n)), Y(n, std::vector<double>(n));
    generateMatrix(X, n, n);
    generateMatrix(Y, n, n);
    std::vector<double> b_l_temp;
    bool converged = false;
    do {
        std::vector<double> a_l_temp;
        for (int i = 0; i < n; ++i) {
            if (b_l_temp.size() != Y[i].size()) {
                std::cerr << "Size mismatch between b_l_temp and Y[" << i << "]" << std::endl;
                continue;
            }
            std::vector<double> diff(b_l_temp.size());
            std::transform(b_l_temp.begin(), b_l_temp.end(), Y[i].begin(), diff.begin(), std::minus<>());
            b_l_temp.push_back(median(b[i]) / std::accumulate(diff.begin(), diff.end(), 0.0));
        }
        for (int i = 0; i < n; ++i) {
            if (a_l_temp.size() != X[i].size()) {
                std::cerr << "Size mismatch between a_l_temp and X[" << i << "]" << std::endl;
                continue;
            }
            

            std::vector<double> diff(a_l_temp.size());
            std::transform(a_l_temp.begin(), a_l_temp.end(), Y[i].begin(), diff.begin(), std::minus<>());
            a_l_temp.push_back(median(b[i]) / std::accumulate(diff.begin(), diff.end(), 0.0));
        }
        converged = true;
        for (int i = 0; i < n; ++i) {
            if (!std::equal(a_l_temp.begin(), a_l_temp.begin() + i + 1, a_l.begin())) {
                converged = false;
                break;
            }
            if (!std::equal(b_l_temp.begin(), b_l_temp.begin() + i + 1, b_l.begin())) {
                converged = false;
                break;
            }
        }
        std::swap(a_l, a_l_temp);
        std::swap(b_l, b_l_temp);
    } while (!converged);
}

int main() {
    int n = std::rand() % 10 + 1;
    std::vector<std::vector<double>> a(n, std::vector<double>(n)), b(n, std::vector<double>(n));
    generateMatrix(a, n, n);
    generateMatrix(b, n, n);
    std::vector<double> a_l(n), b_l(n);
    leastSquares(a, b, a_l, b_l);
    std::cout << "a_l: ";
    for (auto el : a_l) {
        std::cout << el << ' ';
    }
    std::cout << '\n';
    std::cout << "b_l: ";
    for (auto el : b_l) {
        std::cout << el << ' ';
    }
    std::cout << '\n';
    return 0;
}



























































/*
* 
* 
* // LI.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//123

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "LI.h"

using namespace std;
template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<size_t>& weights);


#pragma region Utility



// Структура для точки данных
struct Point {
    vector<double> x; // Независимые переменные
    double y;         // Зависимая переменная
};


//демонстрация точек
void printPoints(const vector<Point>& points) {
    cout << "Сгенерированные точки:" << endl;
    for (size_t i = 0; i < points.size(); ++i) {
        cout << "Точка " << i + 1 << ": ";
        cout << "x = [";
        for (size_t j = 0; j < points[i].x.size(); ++j) {
            cout << points[i].x[j];
            if (j < points[i].x.size() - 1) cout << ", ";
        }
        cout << "] ";
        cout << "y = " << points[i].y << endl;
    }
}



//генерация точек
vector<Point> generatePoints(size_t numPoints = 10, size_t numVariables = 10, double min = 0.0, double max = 10.0) {
    vector<Point> points;

    // Инициализация генератора случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min, max);

    for (size_t i = 0; i < numPoints; ++i) {
        Point point;
        point.x.resize(numVariables);

        // Генерация независимых переменных
        for (size_t j = 0; j < numVariables; ++j) {
            point.x[j] = dist(gen);
        }

        // Генерация зависимой переменной (например, линейная зависимость + шум)
        double noise = dist(gen) * 0.1; // Малый шум
        point.y = 0.0;
        for (size_t j = 0; j < numVariables; ++j) {
            point.y += point.x[j] * (j + 1); // Пример линейной зависимости
        }
        point.y += noise; // Добавление шума

        points.push_back(point);
    }

    return points;
}



#pragma endregion

int main()
{



    //vector<Point> points;
    //int numVariables  = points[1].x.size(), numPoints = points.size();


    //std::cout<<"количество измерений" << numVariables << "  \n"<<"количество точек" << numPoints << endl;


std::vector<double> values = { 1.0, 2.0, 3.0, 4.0, 5.0 };
std::vector<size_t> weights = { 1, 2, 3, 4, 5 };

try {
    double result = weightedMedian(values, weights);
    std::cout << "Weighted Median: " << result << std::endl;
}
catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
}

return 0;


}




double gen_B0(const vector<double>& X, const vector<double>& Y, int N)
{
    double Yy = accumulate(Y.begin(), Y.end(), 0.0) / Y.size(), Xx = accumulate(X.begin(), X.end(), 0.0) / X.size();
    double verh = 0, niz = 0;

    for (int i = 0; i < N; i++)
    {
        verh += (X[i] - Xx) * (Yy * Y[i] - Xx * Y[i]);
        niz += (X[i] - Xx) * (X[i] - Xx);
    }


    return verh / niz;

}


double gen_A0(const vector<double>& X, const vector<double>& Y, double b, int N)
{
    vector <double> medMass;
    for (int i = 0; i < N; i++)
    {
        medMass.push_back(abs(X[i]) * ((Y[i] - b) / X[i]));
    }



    double result = weightedMedian(medMass, medMass);
    return result;

}

//сумма модулей. (2)
double rss(int N, double* Y, double* X, double A, double B) {
    double sum_rss = 0.0;

    for (int i = 0; i < N; ++i) {
        sum_rss += Y[i] - A * X[i] - B;
    }

    return sum_rss;
}


//взвешенная медиана
template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<T>& weights)

{
    if (values.empty() || weights.empty()) {
        throw std::invalid_argument("Values and weights must be non-empty");
    }

    if (values.size() != weights.size()) {
        throw std::invalid_argument("Sizes of values and weights must match");
    }

    // Преобразуем веса в кумулятивные суммы
    std::vector<size_t> cumulativeWeights(weights.size());
    std::partial_sum(weights.begin(), weights.end(), cumulativeWeights.begin());

    size_t totalWeight = cumulativeWeights.back();

    // Найти медианный индекс
    size_t medianIndex;
    if (totalWeight % 2 == 0) {
        // Четное число элементов
        medianIndex = totalWeight / 2 - 1;
    }
    else {
        // Нечетное число элементов
        medianIndex = totalWeight / 2;
    }

    // Бинарный поиск медианного значения
    auto it = std::upper_bound(cumulativeWeights.begin(), cumulativeWeights.end(), medianIndex);
    size_t rank = std::distance(cumulativeWeights.begin(), it);

    // Определить медианное значение
    T medianValue = values[rank];

    return medianValue;
}




//
//
//
//vector<Point> points = {
//   {{1.0, 2.0}, 3.0},
//   {{2.0, 1.0}, 4.0},
//   {{3.0, 3.0}, 6.0},
//   {{4.0, 5.0}, 8.0},
//   {{5.0, 4.0}, 9.0}
//};*/