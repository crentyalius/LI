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


using namespace std;
template<typename T> T weightedMedian(const std::vector<T>& values, const std::vector<T>& weights);//взвешивание медианы
void sortValuesAndWeights(std::vector<double>& values, std::vector<double>& weights);//сортировка по значениям


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
    int N = 5,k=0;
    vector <double> X = { 0.450, 0.50, 0.60, 2.0, 1.2 };
    vector <double> Y = { 4.0, 5.0, 7.0, 10.0, 10.0 };
    vector <int> Index;
    vector <int> IndexStory;
     std::string filename = openFileDialog();
     if (filename.empty()) {
        std::cerr << "Файл не выбран!" << std::endl;
        return 1;
    }
    fillArraysFromFile(filename, array1, array2);
    
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

std::string openFileDialog() {
    OPENFILENAME ofn;
    char fileName[MAX_PATH] = "";

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFilter = "Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    ofn.lpstrDefExt = "txt";

    if (GetOpenFileName(&ofn)) {
        return fileName; // Возвращаем путь к выбранному файлу
    }
    return ""; // Если файл не выбран, возвращаем пустую строку
}


bool fillArraysFromFile(const std::string& filename, std::vector<double>& array1, std::vector<double>& array2) {
    std::ifstream file(filename); // Открываем файл для чтения
    if (!file.is_open()) {        // Проверяем, удалось ли открыть файл
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return false;
    }

    double a, b;
    while (file >> a >> b) { // Читаем пары значений
        array1.push_back(a); // Добавляем первое значение в первый массив
        array2.push_back(b); // Добавляем второе значение во второй массив
    }

    file.close(); // Закрываем файл
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


    double Yy=0, Xx=0;
    for (int i = 0; i < N; i++)
    {
        Yy += Y[i];
        Xx += X[i];
    }
    Yy /=  (double)Y.size();
    Xx/= (double)X.size();
    



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

    double Xx=0, Yy=0, XXX=0,XY=0;

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

    double result = weightedMedian( medMass, x );

     
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
