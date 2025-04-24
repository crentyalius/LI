#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

// Структура для хранения пары (X, Y)
struct Point {
    double X;
    double Y;
};

// Функция для вычисления взвешенной медианы
double weightedMedian(const std::vector<double>& samples, const std::vector<double>& weights) {
    // Проверка на одинаковый размер
    if (samples.size() != weights.size()) {
        throw std::invalid_argument("Samples and weights must have the same size.");
    }

    // Создаем вектор индексов для сортировки
    std::vector<size_t> indices(samples.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Сортируем индексы по значениям samples
    std::sort(indices.begin(), indices.end(), [&samples](size_t i, size_t j) {
        return samples[i] < samples[j];
        });

    // Вычисляем половину суммы весов
    double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
    double halfWeight = totalWeight / 2.0;

    // Находим взвешенную медиану
    double cumulativeWeight = 0.0;
    for (size_t i = 0; i < indices.size(); ++i) {
        cumulativeWeight += weights[indices[i]];
        if (cumulativeWeight >= halfWeight) {
            return samples[indices[i]];
        }
    }

    // Если все веса нулевые (маловероятно)
    return samples.back();
}

// Функция для вычисления LAD-регрессии
void ladRegression(const std::vector<Point>& points, double& a, double& b, double tolerance = 1e-6, int maxIterations = 100) {
    if (points.empty()) {
        throw std::invalid_argument("Points vector is empty.");
    }

    // Шаг 1: Инициализация через метод наименьших квадратов (LS)
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
    for (const auto& p : points) {
        sumX += p.X;
        sumY += p.Y;
        sumXY += p.X * p.Y;
        sumX2 += p.X * p.X;
    }

    double N = points.size();
    double meanX = sumX / N;
    double meanY = sumY / N;

    // Вычисляем начальное b (сдвиг) через LS
    b = (sumXY - sumX * meanY) / (sumX2 - sumX * meanX);
    // Вычисляем начальное a (наклон) через LS (альтернативная формула)
    a = meanY - b * meanX;

    // Итеративный процесс
    int iteration = 0;
    int j = 0; // Индекс точки, через которую проходит текущая реберная линия

    while (iteration < maxIterations) {
        // Шаг 2: Обновление параметра a через взвешенную медиану
        std::vector<double> samplesA;
        std::vector<double> weightsA;
        for (const auto& p : points) {
            if (p.X != 0.0) { // Избегаем деления на ноль
                samplesA.push_back((p.Y - b) / p.X);
                weightsA.push_back(std::abs(p.X));
            }
        }

        double aNew = weightedMedian(samplesA, weightsA);

        // Находим индекс j точки, через которую проходит реберная линия
        for (j = 0; j < points.size(); ++j) {
            if (std::abs((points[j].Y - b) / points[j].X - aNew) < tolerance) {
                break;
            }
        }

        // Шаг 3: Преобразование координат и обновление параметра b
        double Xj = points[j].X;
        double Yj = points[j].Y;

        // Сдвигаем координаты
        std::vector<double> shiftedY;
        std::vector<double> weightsB;
        for (const auto& p : points) {
            shiftedY.push_back(p.Y - aNew * p.X);
            weightsB.push_back(1.0); // Веса для медианы (можно адаптировать)
        }

        // Обновляем b через медиану (не взвешенную, так как веса равны)
        double bNew = weightedMedian(shiftedY, weightsB);

        // Проверка сходимости
        if (std::abs(aNew - a) < tolerance && std::abs(bNew - b) < tolerance) {
            a = aNew;
            b = bNew;
            break;
        }

        // Обновляем параметры
        a = aNew;
        b = bNew;
        iteration++;
    }

    if (iteration == maxIterations) {
        std::cerr << "Достигнуто максимальное число итераций. Результат может быть неточным." << std::endl;
    }
}

int main() {
    // Пример данных
    std::vector<Point> points = {
        {0.45, 4},
        {0.5, 5},
        {0.6, 7},
        {2.0, 10},
        {1.2, 10}
    };
    setlocale(LC_ALL, "Russian");
    double a, b;
    ladRegression(points, a, b);

    std::cout << "Наклон (a): " << a << std::endl;
    std::cout << "Сдвиг (b): " << b << std::endl;

    return 0;
}
