// Li&ARCH.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

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
#include <windows.h>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<T>& weights);//взвешивание медианы
void sortValuesAndWeights(std::vector<double>& values, std::vector<double>& weights);//сортировка по значениям

std::string openFileDialog();
bool fillArraysFromFile(const std::string& filename, std::vector<double>& array1, std::vector<double>& array2, int n );

double gen_B0(const vector<double>& X, const vector<double>& Y, int N);//генерация b0
double gen_A0(const vector<double>& X, const vector<double>& Y, double b, int N);//генерация a0


vector<int> Indexate(const vector<double>& X, const vector<double>& Y, double a0, double b0, int N);//вывод индексов




int main()
{


    setlocale(LC_ALL, "ru_RU.UTF-8");
    const double epsilon = 0.001;
    vector <double> bk;
    vector <double> ak;
    //инициализация массивов
    int N = 5, k = 0;
    vector <double> X = { 0.450, 0.50, 0.60, 2.0, 1.2 };
    vector <double> Y = { 4.0, 5.0, 7.0, 10.0, 10.0 };
    vector <int> Index;
    vector <int> IndexStory;






    std::string filename = "arrays.txt";
    
    fillArraysFromFile(filename, X, Y,N);

    N = X.size();


    //инициализация первичных a0 и b0
    bk.push_back(gen_B0(X, Y, N));//работает 

    ak.push_back(gen_A0(X, Y, bk.back(), N));// работает.

    Index = Indexate(X, Y, ak.back(), bk.back(), N);







    // printf("b%d=%f \n", k, bk.back());
   //  printf("a%d=%f \n\n\n", k, ak.back());

    k++;

    double Xj;
    Xj = X[Index.back()];
    // Xj = X[Index[0]];


    for (int i = 0; i < N; i++)
        X[i] -= Xj;


    bk.push_back(bk.back() + ak.back() * Xj);
    ak.push_back(gen_A0(X, Y, bk.back(), N));

    // printf("b%d=%f \n",k, bk.back());
    // printf("a%d=%f \n\n\n",k, ak.back());


    Index.clear();
    Index = Indexate(X, Y, ak.back(), bk.back(), N);
    IndexStory.push_back(Index.back());
    k++;


    for (; abs(ak.back() - ak[ak.size() - 2]) > epsilon; k++)
    {
        Xj = X[Index.back()];
        //Xj = X[Index[0]];


        for (int i = 0; i < N; i++)
            X[i] -= Xj;


        bk.push_back(bk.back() + ak.back() * Xj);
        ak.push_back(gen_A0(X, Y, bk.back(), N));



        Index.clear();
        Index = Indexate(X, Y, ak.back(), bk.back(), N);


        bk.push_back(bk.back() - ak.back() * Xj);
        // IndexStory.push_back(Index.back());


         //std::printf("bLS=%f \n", (gen_B0(X, Y, N)));
        // std::printf("b%d=%f \n", k, bk.back());
        // std::printf("a%d=%f \n\n\n", k, ak.back());
    }
    std::printf("b%d=%f \n", k, bk.back());
    std::printf("a%d=%f \n\n\n", k, ak.back());

    return 0;


}





bool fillArraysFromFile(const std::string& filename, std::vector<double>& array1, std::vector<double>& array2,int n=0) {
    
        ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Failed to open file for reading.\n";
        return false;
    }

    // Reading the values from the file and initializing the
    // array
    for (int i = 0; i < n; ++i) {
        infile >> array1[i]>> array2[i];
    }

    // Closing the file
    infile.close();
    return true;
}


vector<int> Indexate(const vector<double>& X, const vector<double>& Y, double a0, double b0, int N)
{
    double min = abs(a0 - (Y[1] - b0) / X[1]);

    vector<int> index;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {

            if (a0 == (Y[i] - b0) / X[j])
            {
                index.push_back(i);
                //index[i] = j;
                //printf("X[%d]  подходит Y[%d]\n", j,i);
                break;
            }
            //else printf("NO %f \n", a0 -((Y[i] - b0) / X[j]));



        }


    }

    return index;

}

double gen_B0(const vector<double>& X, const vector<double>& Y, int N = 0)
{
    if (N == 0)
        N = X.size();


    double Yy = 0, Xx = 0;
    for (int i = 0; i < N; i++)
    {
        Yy += Y[i];
        Xx += X[i];
    }
    Yy /= (double)Y.size();
    Xx /= (double)X.size();




    double up = 0.0, down = 0.0;
    for (int i = 0; i < N; i++)
    {
        up += (X[i] - Xx) * (Yy * Y[i] - Xx * Y[i]);
        down += (X[i] - Xx) * (X[i] - Xx);
    }


    return up / down;

}

double gen_B0_Alt(const vector<double>& X, const vector<double>& Y, int N = 0)
{
    if (N == 0)
        N = X.size();

    double Xx = 0, Yy = 0, XXX = 0, XY = 0;

    for (int i = 0; i < N; i++)
    {
        Xx += X[i];
        Yy += Y[i];
        XXX += X[i] * X[i];
        XY += X[i] * Y[i];
    }

    return (N * XY - Xx * Yy) / (N * XXX - Xx * Xx);


}


double gen_A0(const vector<double>& X, const vector<double>& Y, double b, int N)
{
    vector <double> medMass, x;
    for (int i = 0; i < N; i++)
    {
        x.push_back(abs(X[i]));
        medMass.push_back(((Y[i] - b) / X[i]));
        //printf("%f ", ((Y[i] - b) / X[i]));
    }

    //sortValuesAndWeights(medMass, x);

    double result = weightedMedian(medMass, x);


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
template<typename T> T weightedMedian(const std::vector<T>& values1, const std::vector<T>& weights1)

{
    if (values1.empty() || weights1.empty()) {
        throw std::invalid_argument("Values and weights must be non-empty");
    }

    if (values1.size() != weights1.size()) {
        throw std::invalid_argument("Sizes of values and weights must match");
    }

    vector<T>values = sortValuesAndWeights(values1, weights1);
    vector<T>weights = sortValuesAndWeights(values1, weights1, 1);
   
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

template<typename T> T weightedMedianRecursed( std::vector<T>& values,  std::vector<T>& weights, T p, T r)
{
    


    if (p == r)
        return p;
    if (r - p = 1)
    {
        if (weights[p] == weights[r])
            return (values[p] + values[r]) / 2;


        if (weights[p] > weights[r])
            return values[p];
        else
            return values[r];
           

    }
    T wl= reduce(weights.begin(), weights[p])/ reduce(weights.begin(), weights.end()), wg= reduce(weights[r], weights.end()/ reduce(weights.begin(), weights.end()));
    if (wl < 0.5 || wg < 0.5)
    {


    }

   

   
}


template<typename T> void exchange(T& A, T& B)
{
    A = A + B;
    B = A - B;
    A = A - B;
    return A;


}


template<typename T> int Partition_Around(std::vector<T>& A, std::vector<T>& W, T p, T r, T x)
{
    int  i = p - 1;
    for (int j = p; j < r; j++)
    {
        if (A[j] <= x)
        {
            i++;
            exchange(A[i], A[j]);
            exchange(W[i], W[j]);
        }
        exchange(A[i + 1], A[r]);
        exchange(W[i + 1], W[r]);
    }
    return i + 1;
}
 
template<typename T> T Select( std::vector<T>& A , std::vector<T>& W, T p, T r, T i)
{
    while ((r - p + 1) % 5 != 0)
    {
        for (int j = p + 1; j <= r; j++)// put the minimum into AŒp�
        {
            if (A[p] > A[j])
            {
                exchange(A[p], A[j]);
                exchange(W[p], W[j]);
            }
        }
        if (i == 1)// If we want the minimum of AŒp W r�, we’re done.
            return A[p];
        p = p + 1;// Otherwise, we want the .i  1/st element of AŒp C 1 W r�
        i = i - 1;
    }
    T g = (r - p + 1) / 5;
    for (int j = p; j < p + g; j++)
    {
        //создание временного массива для сортировки
        vector <T> tempvalues = { A[j],A[j + g] ,A[j + 2 * g],A[j + 3 * g],A[j + 4 * g] };  
        vector <T> tempweights = { W[j],W[j + g] ,W[j + 2 * g],W[j + 3 * g],W[j + 4 * g] };
        // Создаем пары (значение, вес)
        vector <std::pair<T, T>> data;
        for (size_t i = 0; i < tempvalues.size(); ++i) {
            data.push_back({ tempvalues[i], tempweights[i] });
        }
        //окончание  генерации массива для  сортировки
        sort(data.begin(), data.end(), [](const std::pair<T, T>& a, const std::pair<T, T>& b) {
            return a.first < b.first;
            });

        // Разделяем отсортированные значения и веса
        for (size_t i = 0; i < data.size(); ++i) 
        {
            tempvalues[i] = data[i].first;
            tempweights[i] = data[i].second;
        }
        //возврат инициализированных значений
        A[j] = tempvalues[1];
        A[j+g] = tempvalues[2];
        A[j + 2*g] = tempvalues[3];
        A[j +3* g] = tempvalues[4];
        A[j + 4* g] = tempvalues[5];

        W[j] = tempweights[1];
        W[j + g] = tempweights[2];
        W[j + 2 * g] = tempweights[3];
        W[j + 3 * g] = tempweights[4];
        W[j + 4 * g] = tempweights[5];
    }
        T X = Select(A, W, p + 2 * g, p + 3 * g - 1, g / 2);
        T q = Partition_Around(A, W, p, r, X);
        T k = q - p + 1;
        if (i == k)
        {
            return A[q];
        }
        else if (i < k)
        {
            return Select(A, W, p, q - 1, i);
        }
        else
            return Select(A, W, q+1,r,i - k);

}


template<typename T> vector<T>  sortValuesAndWeights(const std::vector<T>& values1, const std::vector<T>& weights1,int rg=0) {
    if (values1.empty() || weights1.empty() || values1.size() != weights1.size()) {
        std::cerr << "Ошибка: неверные входные данные!" << std::endl;
        return false;
    }
    vector<T>values = values1;
    vector<T>weights = weights1;
    // Создаем пары (значение, вес)
    std::vector<std::pair<T, T>> data;
    for (size_t i = 0; i < values.size(); ++i) {
        data.push_back({ values[i], weights[i] });
    }

    // Сортируем пары по значениям
    std::sort(data.begin(), data.end(), [](const std::pair<T, T>& a, const std::pair<T, T>& b) {
        return a.first < b.first;
        });

    // Разделяем отсортированные значения и веса
    for (size_t i = 0; i < data.size(); ++i) {
        values[i] = data[i].first;
        weights[i] = data[i].second;
    }
    if( rg==0)
        return values;
    return weights;
}
