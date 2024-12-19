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


using namespace std;
template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<T>& weights);
double gen_B0(const vector<double>& X, const vector<double>& Y, int N);

double gen_A0(const vector<double>& X, const vector<double>& Y, double b, int N);
vector<int> Indexate(const vector<double>& X, const vector<double>& Y, double a0, double b0, int N);

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


    /*
    std::vector<double> values = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::vector<size_t> weights = { 1, 2, 3, 4, 5 };

    try {
        double result = weightedMedian(values, weights);
        std::cout << "Weighted Median: " << result << std::endl;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }*/

    //инициализация массивов
    vector <double> X;
    vector <double> Y;
    vector <int> Index;
    int N = 100;

    for (int i = 0; i < N; i++)
    {
        X[i] = i;
        Y[i] = i;


    }
    //инициализация первичных a0 и b0
   double b0= gen_B0(X, Y, N);
   double a0 = gen_A0(X, Y, b0, N);
   
   Index = Indexate(X, Y, a0, b0, N);




    return 0;


}

vector<int> Indexate(const vector<double>& X, const vector<double>& Y, double a0, double b0, int N)
{
    vector<int> index;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (a0 == (Y[i] - b0) / X[j])
            {
                index[i] = j;
                break;
            }

        }


    }
    return index;

}

double gen_B0(const vector<double>& X, const vector<double>& Y,  int N)
{
    double Yy, Xx;

    Yy= accumulate(Y.begin(), Y.end(), 0)/Y.size();
    Xx = accumulate(X.begin(), X.end(), 0) / X.size();

    double up =0, down =0;

    for (int i = 0; i < N; i++)
    {
        up += (X[i] - Xx) * (Yy * Y[i] - Xx * Y[i]);
        down += (X[i] - Xx) * (X[i] - Xx);
    }


    return up / down;

}


double gen_A0(const vector<double>& X, const vector<double>& Y,double b, int N)
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

