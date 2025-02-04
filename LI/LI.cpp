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
    int N = 5;
    vector <double> X= { 0.450, 0.50, 0.60, 2.0, 1.2 };
    vector <double> Y= { 4.0, 5.0, 7.0, 10.0, 10.0 };
    vector <int> Index;
    

   /* for (int i = 0; i < N; i++)
    {
        X.push_back(i);
        Y.push_back(i);


    }*/

    //инициализация первичных a0 и b0
   double b0= gen_B0(X, Y, N);
   double a0 = gen_A0(X, Y, b0, N);// не до конца понятно, корректно ли работает.
   
   Index = Indexate(X, Y, a0, b0, N);

   int z = Index.size();

   for (int i = 0; i < Index.size(); i++)
   {

       printf("%d ", Index[i]);

   }
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
                index.push_back(i);
                //index[i] = j;
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

void sortValuesAndWeights(std::vector<double>& values, std::vector<double>& weights) {
    if (values.empty() || weights.empty() || values.size() != weights.size()) {
        std::cerr << "Ошибка: неверные входные данные!" << std::endl;
        return;
    }

    // Создаем пары (значение, вес)
    std::vector<std::pair<double, double>> data;
    for (size_t i = 0; i < values.size(); ++i) {
        data.push_back({ values[i], weights[i] });
    }

    // Сортируем пары по значениям
    std::sort(data.begin(), data.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
        });

    // Разделяем отсортированные значения и веса
    for (size_t i = 0; i < data.size(); ++i) {
        values[i] = data[i].first;
        weights[i] = data[i].second;
    }
}