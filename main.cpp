
#include "triangle.h"
#include "gauss.h"
#include <iomanip> // setprecicion for cout Debug

INITIALIZE_EASYLOGGINGPP

double_t func_calculate_rho(double_t h);

double_t func_calculate_q(double_t t);

double_t* func_multiply_matrix_and_vector(double_t** matrix, double_t* vect, uint_fast32_t dim);

double_t ** func_substract_two_matrices(double_t** matrix1, double_t** matrix2, uint_fast32_t dim, double_t devide_element);

using namespace std;

int main() {
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, "%level: %msg");
    LOG(INFO) << "Starting programm...";


    /*****************************************************************************************************/
    /* Тест вычисления матрицы K и вектора F, а также определителя матрицы A и R с чертой. (стр. 95 учебник)*/
    /*const points test_Point_A = (points){ 8.0, 3 };
    const points test_Point_B = (points){ 8.5, 4 };
    const points test_Point_C = (points){ 7.5, 4 };
    triangles test_triangle = {test_Point_A, test_Point_B, test_Point_C};
    const double_t test_A = 40.0; // heat coeff
    const double_t test_INIT_TEMPERATURE = 300.0;
    const double_t test_H = 0.5;
    double_t test_Q = test_INIT_TEMPERATURE * test_H;

    double_t test_res_A = test_triangle.GetMatrixADeterminant(); // determinant f matrix A
    cout << "2A = " << 2 * test_res_A << endl;

    double_t test_R = test_triangle.GetR();
    cout << "R = " << test_R << endl;

    double_t test_something = (2.0 * PI * test_R * test_A) / (4.0 * test_res_A);
    cout << "2 * pi * R * a / (4 * A) = " << test_something << endl;

    double_t** test_K;
    test_K = new double_t*[3];
    for (int i = 0; i < 3; i++) {
        test_K[i] = new double_t[3];
    }
    test_triangle.Matrix_K(test_K);


    cout << "Matrix K: " << endl;
    for (auto i = 0; i < 3; i ++) {
        for (auto j = 0; j < 3; j ++)
            cout << setprecision(10) << test_K[i][j] << "\t" ;
        cout << endl;
    }

    double_t * test_F = new double[3];
    test_triangle.Column_F(test_F ,150);

    cout << "F = :";
    for (auto i = 0; i < 3; i++)
        cout << test_F[i] << "\t";*/
    /*****************************************************************************************************/


    /*****************************************************************************************************/
    /* Объявление переменных */

    /**
     * Вывод итоговой матрицы жесткости
     */
    ofstream file1("Main_matrix.dat");
    /**
     * Вывод распределения температур на треугольнике во времени
     */
    ofstream file2("Triangle.dat");
    /**
     * Координаты треугольника
     */
    const points Point_A = (points) {0, 0};
    const points Point_B = (points) {1, 0};
    const points Point_C = (points) {1, 2};
    /**
     * Инстанс класса, отвечающий за триангулацию треугольника
     */
    IsoscelesTriangleGrid triangle(Point_A, Point_B, Point_C, STEP_X);
    /**
     * Инстанс класса, отвечающий за решение СЛАУ на каждом элементе методом Гаусса
     */
    GaussMethod equation;
    /**
     * Шаг по временному слою
     */
    double_t tau = TAU;
    /**
     * Число точек
     */
    uint_fast32_t n = triangle.GetGreed(false) + 1;
    /**
     * Число узлов
     */
    uint_fast32_t k = triangle.triangles_array.size();
    /**
     * Тепловой поток, заданный на границе
     */
    double_t q = func_calculate_q(tau);
    /**
     * Столбец правых частей для элемента
     */
    double_t *F = new double_t[3];
    /**
     * Значения температуры в узлах
     */
    double_t *Temp = new double_t[n];
    /**
     * Cюда записываем те температуры с предыдыщуего слоя которые нас интересуют для данного конечного элемента
     */
    double_t Fi[3];
    /**
     * Промежуточная переменная для записи произведения матрицы и столбца
     */
    double_t *Resultic = new double_t[3];
    /**
     * Вектор правых частей итоговой системы
     */
    double_t *Result = new double_t[n];
    /**
     * Итоговая матрица жесткости системы
     */
    double_t **Main_Matrix= new double_t *[n];
    /**
     * Матрица коэффициентов элемента K
     */
    double_t **K= new double_t *[3];
    /**
     * Матрица коэффициентов элемента C
     */
    double_t **C = new double_t *[3];
    /**
     * Массив номеров строк в итоговом векторе правых частей Result,
     * к которому нужно будет прибавить полученные правые части F на элементе
     */
    uint_fast32_t ind[3];
    /**
     * Результат вычитания двух матриц. Вспомогательная переменная
     */
    double_t ** substracted_matrix = new double_t* [3];
    /*****************************************************************************************************/



    /* Debug */
    q = 100000;

    LOG(INFO) << "Number of dots = " << n;
    LOG(INFO) << "Number of triangles = " << k;
    LOG(INFO) << "Heat flux q = " << q;
    /* End debug */

    /* Задаем начальное условие распределения температуры */
    for (uint_fast32_t i = 0; i < n; i++) {
        Temp[i] = INITIAL_TEMPERATURE;
    }

    for (uint_fast32_t i = 0; i < 3; i++) {
        substracted_matrix[i] = new double_t[3];
        for (uint_fast32_t j = 0; j < 3; j++)
            substracted_matrix[i][j] = 0;
    }
    /* Ищем решения по временным слоям */
    for (uint_fast32_t global_tau =0; global_tau < 10; global_tau++) {

        /* Инициализация матрицы жесткости и вектора правой части итоговой системы */
        for (uint_fast32_t i = 0; i < n; i++) {
            Main_Matrix[i] = new double[n];
            Result[i] = 0;
            for (int j = 0; j < n; j++) {
                Main_Matrix[i][j] = 0;
            }
        }

        /* Инициализация матрц K, C и вектора правых частей F для элемента */
        for (uint_fast32_t i = 0; i < 3; i++) {
            K[i] = new double_t[3];
            C[i] = new double_t[3];
            F[i] = 0;
        }

        /* Идем по всем элементам */
        for (uint_fast32_t i = 0; i < k; i++) {
            /* Вычисляем матрицы коэффициентов K, C и вектор правых частей F для элемента k */
            triangle.triangles_array[i].Matrix_K(K);
            triangle.triangles_array[i].Matrix_C(C);
            triangle.triangles_array[i].Column_F(F, q);

            ind[0] = triangle.triangles_array[i].first_point.point_num;
            ind[1] = triangle.triangles_array[i].second_point.point_num;
            ind[2] = triangle.triangles_array[i].third_point.point_num;

            /* В Fi записываем температуры Temp на текущем временном слое тех узлов, которые входят в элемент k */
            for (uint_fast32_t j = 0; j < 3; j++) {
                Fi[j] = Temp[ind[j]];
            }

            Resultic = func_multiply_matrix_and_vector(C, Fi, 3);
            substracted_matrix = func_substract_two_matrices(C, K, 3, 2.0/tau);
            /* Заполняем итоговую матрицу жесткости и вектор правой части полученными на элементе значениями */
            for (uint_fast32_t j = 0; j < 3; j++) {
                Result[ind[j]] += (Resultic[j] / tau) + F[j];
                for (uint_fast32_t l = 0; l < 3; l++) {
                    Main_Matrix[ind[j]][ind[l]] += (C[j][l] / tau) + K[j][l];
                }
            }

            /* Debug */
/*            file1 << k << endl;
            for (int iii = 0; iii < n; iii++) {
                for (int j = 0; j < n; j++) {
                    file1 << Main_Matrix[iii][j] << " ";
                }
                file1 << endl;
            }
            file1 << " ------------------------------------------" << endl;*/
            /* End debug */
        }

        /* Debug */
/*        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                file1 << Main_Matrix[i][j] << " ";
            }
            file1 << endl;
        }
        file1.close();*/
        /* End debug */

        Temp = equation.GaussSolve(Main_Matrix, Result, n);

        /* Debug */
        file2 << "----------------------\n";
        file2 << "Triangle temperature: \n";
        int tmp = 0;
        char linebr[10];
        for (auto i = 0; i < 11; i++)
            file2 << setw(15) << setprecision(9) << Temp[i] << "\t";
        file2 << "\n";
        int nn = 11;
        for (auto i = 1; i < 11; i++) {
            nn -= 1;
            linebr[tmp] = ' ';
            tmp++;
            for (auto ii = 0; ii < tmp; ii++) {
                file2 << setw(15) << linebr[ii] << "\t";
            }
            for (auto j = 0; j < nn; j++) {
                file2 << setw(15) << setprecision(9) << Temp[j + nn] << "\t";
            }
            file2 << "\n";
        }
        /* End debug */
    }
    return 0;
}

double_t func_calculate_rho(double_t h)
{
    return (1.2 * exp(-0.00013 * h));
}

double_t func_calculate_q(double_t t)
{
    // Тепловой поток при скорости 3000 м/c на расстоянии 1 м от Земли
    double_t S = 1;
    double_t rho = func_calculate_rho(1);
    double_t V = 3000;
    double_t pi = 3.1415926;
    double_t r = 1;
    double_t H = 2;

    return ((S * rho * V * V * V) / (Thermal_Conductivity * pi * r * sqrt(r * r + H * H)));
}

double_t* func_multiply_matrix_and_vector(double_t** matrix, double_t* vect, uint_fast32_t dim)
{
    double_t* result_vector;
    result_vector = new double_t[3];
    for (auto i = 0; i < 3; i++)
        result_vector[i] = 0;

    for (auto i = 0; i < dim; i++) {
        for (auto j = 0; j < dim; j++) {
            result_vector[i] += vect[j] * matrix[j][i];
        }
    }
    return result_vector;
}

double_t ** func_substract_two_matrices(double_t** matrix1, double_t** matrix2, uint_fast32_t dim, double_t devide_element) {
    double_t ** result_matrix = new double_t* [dim];
    for (auto i = 0; i < dim; i++)
        result_matrix[i] = new double_t[dim];

    for (auto i = 0; i < dim; i++)
        for (auto j = 0; j < dim; j++)
            result_matrix[i][j] = (devide_element * matrix1[i][j]) - matrix2[i][j];

    return result_matrix;
}