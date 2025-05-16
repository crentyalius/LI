#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <clocale>
#include <limits>

struct DataPoint {
    double x; // Р·РЅР°С‡РµРЅРёРµ
    double y; // РІРµСЃ

    DataPoint(double xVal, double yVal) : x(xVal), y(yVal) {}
    DataPoint() : x(0.0), y(0.0) {} // РџРѕ СѓРјРѕР»С‡Р°РЅРёСЋ
};

DataPoint weightedMedian(const std::vector<DataPoint>& dataPoints) {
    std::vector<DataPoint> sortedData = dataPoints;
    std::sort(sortedData.begin(), sortedData.end(),
        [](const DataPoint& a, const DataPoint& b) {
            return a.x < b.x;
        });

    double totalWeight = 0.0;
    for (const auto& p : sortedData) {
        totalWeight += p.y;
    }

    double cumulativeWeight = 0.0;
    const double halfWeight = totalWeight / 2.0;

    for (const auto& point : sortedData) {
        cumulativeWeight += point.y;
        if (cumulativeWeight >= halfWeight) {
            return point;
        }
    }

    return sortedData.back();
}

double calculateInitialB(const std::vector<DataPoint>& points) {
    if (points.empty()) {
        throw std::invalid_argument("Data points cannot be empty");
    }

    size_t N = points.size();
    double sum_x = 0.0, sum_y = 0.0;

    for (const auto& point : points) {
        sum_x += point.x;
        sum_y += point.y;
    }

    double x_mean = sum_x / N;
    double y_mean = sum_y / N;

    double numerator = 0.0;
    for (const auto& point : points) {
        numerator += (point.x - x_mean) * (y_mean  - x_mean * point.y);
    }

    double denominator = 0.0;
    for (const auto& point : points) {
        denominator += std::pow(point.x - x_mean, 2);
    }

    if (denominator == 0) {
        throw std::runtime_error("Denominator is zero (all x values are equal)");
    }
    std::cout << "b0 " << numerator / denominator << std::endl;// return 3.7;
    return abs(numerator / denominator);
}

double calculateA0(const std::vector<DataPoint>& points, double b0) {
    std::vector<DataPoint> weightedValues;
    for (const auto& point : points) {
        if (std::abs(point.x) < 1e-10) continue;
        double value = (point.y - b0) / point.x;
        double weight = std::abs(point.x);
        weightedValues.push_back({ value, weight });
    }
    std::cout << "a0 " << weightedMedian(weightedValues).x << std::endl;
    return weightedMedian(weightedValues).x;
}

struct EdgeParameters {
    size_t index;
    double slope;
    double intercept;
};

EdgeParameters findEdgePoint(double a0, double b0, const std::vector<DataPoint>& points) {
    EdgeParameters result;
    double min_diff = std::numeric_limits<double>::max();

    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        if (std::abs(p.x) < 1e-10) continue;

        double current_a = (p.y - b0) / p.x;
        double diff = std::abs(current_a - a0);

        if (diff < min_diff) {
            min_diff = diff;
            result.index = i;
            result.slope = -p.y / p.x;
            result.intercept = p.y;
        }
    }

    return result;
}

std::vector<DataPoint> transformDataSpace(const std::vector<DataPoint>& points, size_t j) {
    double x_j = points[j].x;
    std::vector<DataPoint> transformed;

    for (const auto& p : points) {
        transformed.push_back({ p.x - x_j, p.y });
    }

    return transformed;
}

std::vector<DataPoint> getMedianData(const std::vector<DataPoint>& points, double bk) {
    std::vector<DataPoint> data;

    for (const auto& p : points) {
        if (std::abs(p.x) < 1e-10) continue;
        double value = (p.y - bk) / p.x;
        double weight = std::abs(p.x);
        data.emplace_back(value, weight);
    }

    return data;
}

void ladRegression(const std::vector<DataPoint>& points, double& a, double& b, double tolerance = 1e-6, int maxIterations = 100)
{
    if (points.empty()) {
        throw std::invalid_argument("Data points cannot be empty");
    }

    // РЁР°Рі 0: РЅР°С‡Р°Р»СЊРЅРѕРµ РїСЂРёР±Р»РёР¶РµРЅРёРµ b0 Рё a0
    double b_k = calculateInitialB(points);       // РЅР°С‡Р°Р»СЊРЅС‹Р№ b0
    double a_k = calculateA0(points, b_k);        // РЅР°С‡Р°Р»СЊРЅС‹Р№ a0

    for (int iter = 0; iter < maxIterations; ++iter) {
        // РЁР°Рі 1: РЅР°Р№С‚Рё Р±Р»РёР¶Р°Р№С€СѓСЋ С‚РѕС‡РєСѓ j
        EdgeParameters edge = findEdgePoint(a_k, b_k, points);
        size_t j = edge.index;
        double x_j = points[j].x;
        double y_j = points[j].y;

        // РЁР°Рі 2: СЃРґРІРёРіР°РµРј РґР°РЅРЅС‹Рµ РїРѕ j
        std::vector<DataPoint> shiftedPoints = transformDataSpace(points, j);

        // РЁР°Рі 3: РІС‹С‡РёСЃР»СЏРµРј РЅРѕРІРѕРµ Р·РЅР°С‡РµРЅРёРµ b'k = b_k + a_k * x_j
        double b_k_shifted = b_k + a_k * x_j;

        // РЁР°Рі 4: РјРµРґРёР°РЅРЅР°СЏ РѕС†РµРЅРєР° РґР»СЏ РЅРѕРІРѕРіРѕ РЅР°РєР»РѕРЅР° a
        std::vector<DataPoint> medianData = getMedianData(shiftedPoints, b_k_shifted);
        double a_k1 = weightedMedian(medianData).x;

        // РЁР°Рі 5: РїРµСЂРµСЃС‡С‘С‚ b (РІРѕР·РІСЂР°С‚ Рє РёСЃС…РѕРґРЅРѕР№ СЃРёСЃС‚РµРјРµ РєРѕРѕСЂРґРёРЅР°С‚)
        double b_k1 = y_j - a_k1 * x_j;

        // РџСЂРѕРІРµСЂРєР° СЃС…РѕРґРёРјРѕСЃС‚Рё
        if (std::abs(a_k1 - a_k) < tolerance && std::abs(b_k1 - b_k) < tolerance) {
            a = a_k1;
            b = b_k1;
            return;
        }

        // РћР±РЅРѕРІР»СЏРµРј РїР°СЂР°РјРµС‚СЂС‹
        a_k = a_k1;
        b_k = b_k1;
    }

    // РџРѕСЃР»Рµ РёС‚РµСЂР°С†РёР№, РІРѕР·РІСЂР°С‰Р°РµРј РїРѕСЃР»РµРґРЅРёРµ Р·РЅР°С‡РµРЅРёСЏ
    a = a_k;
    b = b_k;
}


int main() {
    setlocale(LC_ALL, "Russian");

    std::vector<DataPoint> points = {
        {0.45, 4.0},
        {0.5, 5.0},
        {0.6, 7.0},
        {2.0, 10.0},
        {1.2, 10.0}
    };

    //points = {
    //   {4.0, 0.45},
    //   {5.0,0.5 },
    //   { 7.0,0.6},
    //   { 10.0,2.0},
    //   { 10.0,1.2}
    //};


    points = {
       {1.0, 6.66},
       {2.0, 10.0},
       {3.0, 13.33}
    };


    double a = 0.0, b = 0.0;
    ladRegression(points, a, b);

    std::cout << "Slope (a): " << a << std::endl;
    std::cout << "Intercept (b): " << b << std::endl;

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

