///////////////////////////////////////////////////////////
//Разработать программу решения транспортной задачи.
//Составить оптимальный план перевозок между N складами и
//K магазинами, при котором стоимость перевозок будет
//минимальна.Известна потребность в товаре каждым
//магазином, наличие товара на складах и стоимость
//перевозки единицы продукции с каждого склада до каждого
//магазина.
///////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <iomanip>
#include <string>
#include <algorithm>
using namespace std;

int randomInt(int min, int max);
vector<vector<int>> costMatrixGenerate(int n, int m, int cost_min, int cost_max);
void showMatrix(vector<vector<int>> matrix, int value_max);
double randomDouble(double min, double max); // случайное число от min до max 
vector<int> generate_n_nums_s_sum(int N, int Sum); // генерирует N чисел, сумма которых равна S

class AntAlgorithmTransportProblem {
public:
    AntAlgorithmTransportProblem(const vector<vector<int>> costMatrix, 
                                 const vector<int> supplies, 
                                 const vector<int> demands,
                                 const int iterations,
                                 const double evaporation,
                                 const double Q,
                                 const double odorInit)
                        : costMatrix(costMatrix),
                          supplies(supplies), 
                          demands(demands), 
                          n(supplies.size()),
                          m(demands.size()),
		                  iterations(iterations),
		                  evaporation(evaporation),
                          Q(Q),
                          odorInit(odorInit)


    {
        transportationPlanMatrix = vector<vector<int>>(n, vector<int>(m, 0));
        odorMatrix = vector<vector<double>>(n, vector<double>(m, odorInit));
        srand(time(0));
    }

    void Solve() {

        initialAntDistributionAmongCells();
        costOfPlan = calculateCostOfPlan(transportationPlanMatrix);

        cout << "Изначальная матрица плана перевозок:" << endl;
        showMatrix(transportationPlanMatrix, 100);

        for (int i = 0; i < iterations; i++) {
            distributeOdorOfAntsToOdorMatrix(); // распределение ферромона
            antsDistribution();                 // случайное распределение муравьев (составление плана перевозок, учитывая карту запахов)
            evaporationOdorMatrix();            // испарение ферромона
            cout << "Стоимость перевозок: " << calculateCostOfPlan(transportationPlanMatrix) << endl;
        }

        cout << "Получившаяся матрица плана перевозок:" << endl;
        showMatrix(transportationPlanMatrix, 100);
    }

private:
    const vector<vector<int>> costMatrix;
    const vector<int> supplies;
    const vector<int> demands;
    const int iterations;
    const double evaporation;
    const double Q;
    const double odorInit;

    int costOfPlan;

    vector<vector<int>> transportationPlanMatrix;
    //vector<vector<int>> transportationPlanCostMatrix;
    vector<vector<double>> odorMatrix;
    int n, m; // n - строки, m - столбцы
    double bestCost = INT_MAX;

    void initialAntDistributionAmongCells() {
        for (int i = 0; i < supplies.size(); i++) {
            int supplies_current = supplies[i];
            for (int j = 0; j < demands.size(); j++) {
				int sum_m_j = 0;
                for (int k = 0; k < supplies.size(); k++) {
                    sum_m_j += transportationPlanMatrix[k][j];
                }
                while ((sum_m_j < demands[j]) && (supplies_current > 0)) {
                    transportationPlanMatrix[i][j]++;
                    sum_m_j++;
                    supplies_current--;
                }
            }
        }
    }

    // Распределение муравьев по матрице плана перевозок (transportationPlanMatrix) согласно карте запахов (odorMatrix)
    void antsDistribution() {
        vector<vector<int>> transportationPlanMatrixNew = vector<vector<int>>(n, vector<int>(m, 0));

        for (int i = 0; i < supplies.size(); i++) {
            // вычисление вероятности для каждой из клеток demands в строке supplies[i]
            vector<double> probabilities(demands.size());

            // Подсчет общего количество ферромона на строке i
            int sum_odor = 0;
            for (int j = 0; j < demands.size(); j++) {
                sum_odor += odorMatrix[i][j];
            }
            // Вычисление вероятности появления муравья в каждой клетке строки i
            for (int j = 0; j < demands.size(); j++) {
                probabilities[j] = odorMatrix[i][j] / sum_odor;
            }

            

            // распределение муравьев (supplies_current) по строке с учетом demands (ограничения, табу)
            int supplies_current = supplies[i];
            while (supplies_current > 0) {
                // Вычисляем рандомное число от 0.0 до 1.0
				double currentAnt = randomDouble(0.0, 1.0);

                // |--------|-----|----o---|-----|------|
                // 0                                    1


                // Для каждого диапазона рассматриваем, попадает ли сгенерированное случайное число в диапазон
                // Если попадает, значит муравей отправляется на ячейку [i][j]
                for (int j = 0; j < demands.size(); j++) {
                    // Подсчет суммы первых j+1 членов вектора probabilities
                    double sum_probabilities = 0;
                    for (int k = 0; k < j+1; k++) {
                        sum_probabilities += probabilities[k];
                    }

                    // Подсчет суммы ячеек столбца j
                    int sum_j = 0;
                    for (int k = 0; k < supplies.size(); k++) {
                        sum_j += transportationPlanMatrixNew[k][j];
                    }

                    // Попадает ли наше число currentAnt
                    if ((currentAnt < sum_probabilities) && (sum_j < demands[j])) {
                        transportationPlanMatrixNew[i][j]++;
                        supplies_current--;
                        break;
                    }
                }
            }
        }

        if (calculateCostOfPlan(transportationPlanMatrixNew) < calculateCostOfPlan(transportationPlanMatrix)) {
            transportationPlanMatrix = transportationPlanMatrixNew;
        }
    }

    double calculateCostOfPlan(vector<vector<int>> transportationPlanMatrix) {
        double sum = 0.0;
        for (int i = 0; i < supplies.size(); i++) {
            for (int j = 0; j < demands.size(); j++) {
                sum += costMatrix[i][j] * transportationPlanMatrix[i][j];
            }
        }
		return sum;
    }

    // испарение фермента
    void evaporationOdorMatrix() {
        for (int i = 0; i < supplies.size(); i++) {
            for (int j = 0; j < demands.size(); j++) {
                odorMatrix[i][j] *= (1 - evaporation);
                if (odorMatrix[i][j] < odorInit) {
                    odorMatrix[i][j] = odorInit;
                }
            }
        }
    }

    // нанесение фермента
    void distributeOdorOfAntsToOdorMatrix() {
        double L = calculateCostOfPlan(transportationPlanMatrix);
        for (int i = 0; i < supplies.size(); i++) {
            for (int j = 0; j < demands.size(); j++) {
                odorMatrix[i][j] += (transportationPlanMatrix[i][j] * (Q / L));
            }
        }
    }
};

int main() {
    //                    /*demands, m*/
    //    /*supplies, n*/{2, 3, 4, 2, 5},
    //                   {5, 6, 7, 3, 10},
    //                   {6, 3, 3, 1, 6},

    // costMatrix parameters
    const int n = 4; // количество складов
    const int m = 5; // количество магазинов
    const int cost_min = 1; // минимальная стоимость перевозки
    const int cost_max = 10; // максимальная стоимость перевозки
    // supplies and demands parameter
    const int sum_supplies_or_demands = 80; // cумма единиц товара в складах или сумма единиц потребности в магазинах
    // AntAlgorithmTransportProblem parameters
    const int iterations = 100; // количество итераций
    const double evaporation = 0.5; // 0.0 < evaporation < 1.0
    const double Q = 600.0; // "бочка" с ферромонами
    const double odorInit = 1.0; // значение, которым инициализируют карту запахов


    vector<vector<int>> costMatrix = costMatrixGenerate(n, m, cost_min, cost_max);
    cout << "Матрица стоимости перевозок:" << endl;
	showMatrix(costMatrix, cost_max);

    
    vector<int> supplies = vector<int>(n);
    vector<int> demands = vector<int>(m);
    supplies = /* {20, 30, 30}; */ generate_n_nums_s_sum(n, sum_supplies_or_demands);
    demands = /* {10, 20, 30, 15, 5}; */ generate_n_nums_s_sum(m, sum_supplies_or_demands);


    AntAlgorithmTransportProblem AntAlgorithmTransportProblem(costMatrix, supplies, demands, iterations, evaporation, Q, odorInit);
    AntAlgorithmTransportProblem.Solve();


    system("pause");
    return 0;
}



// Случайное целочисленное число от min до max 
int randomInt(int min, int max) {
    return (rand() % (max - min + 1)) + min;
}

double randomDouble(double min, double max) { // случайное число от min до max 
    return (double)(rand()) * (max - min) / RAND_MAX + min;
}


//                    /*demands, m*/
//    /*supplies, n*/{2, 3, 4, 2, 5},
//                   {5, 6, 7, 3, 10},
//                   {6, 3, 3, 1, 6},
vector<vector<int>> costMatrixGenerate(int n, int m, int cost_min, int cost_max) {
    srand(time(0));
    vector<vector<int>> costMatrix;
    for (int i = 0; i < n; i++) {
        vector<int> v;
        for (int j = 0; j < m; j++) {
            v.push_back(randomInt(cost_min, cost_max));
        }
        costMatrix.push_back(v);
    }
    return costMatrix;
}

void showMatrix(vector<vector<int>> matrix, int value_max) {
    for (auto row : matrix) {
        for (auto val : row) {
            std::cout << setw(std::to_string(value_max).length()) << val << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << endl;
}

vector<int> generate_n_nums_s_sum(int N, int Sum) { // генерирует N чисел, сумма которых равна S
    std::vector<int> numbers;

    // Проверка на корректность входных данных
    if (N <= 0 || Sum <= 0 || Sum < N) {
        std::cerr << "Некорректные входные данные." << std::endl;
        return numbers;
    }

    // Инициализация генератора случайных чисел
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    int remainingSum = Sum;

    // Генерируем N-1 случайное уникальное число
    for (int i = 0; i < N - 1; i++) {
        int maxNumber = remainingSum - (N - 1 - i);
        int randomNum = std::rand() % maxNumber + 1;
        numbers.push_back(randomNum);
        remainingSum -= randomNum;
    }

    // Последний элемент массива
    numbers.push_back(remainingSum);

    // Перемешиваем числа, чтобы сделать их случайными
    std::random_shuffle(numbers.begin(), numbers.end());

    return numbers;
}