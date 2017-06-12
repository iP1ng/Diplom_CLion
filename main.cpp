
#include "triangle.h"
#include "gauss.h"
#include <iomanip> // setprecicion for cout Debug
#include "gauss_12062017.h" // another gauss Debug

INITIALIZE_EASYLOGGINGPP

double func_calculate_rho(double h);

double func_calculate_q(double t);

void func_multiply_matrix_and_vector(double* result_vector, double** matrix, double* vect, uint_fast32_t dim);

void func_substract_two_matrices(double** result_matrix, double** matrix1, double** matrix2, uint_fast32_t dim, double devide_element);

void print(vector< vector<double> > A) {
    ofstream file_slau("SLAU.dat", std::ofstream::out | std::ofstream::app);
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            file_slau << A[i][j] << "\t";
            if (j == n-1) {
                file_slau << "| ";
            }
        }
        file_slau << "\n";
    }
    cout << endl;
    file_slau.close();
}

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
    const double test_A = 40.0; // heat coeff
    const double test_INIT_TEMPERATURE = 300.0;
    const double test_H = 0.5;
    double test_Q = test_INIT_TEMPERATURE * test_H;

    double test_res_A = test_triangle.GetMatrixADeterminant(); // determinant f matrix A
    cout << "2A = " << 2 * test_res_A << endl;

    double test_R = test_triangle.GetR();
    cout << "R = " << test_R << endl;

    double test_something = (2.0 * PI * test_R * test_A) / (4.0 * test_res_A);
    cout << "2 * pi * R * a / (4 * A) = " << test_something << endl;

    double** test_K;
    test_K = new double*[3];
    for (int i = 0; i < 3; i++) {
        test_K[i] = new double[3];
    }
    test_triangle.Matrix_K(test_K);


    cout << "Matrix K: " << endl;
    for (auto i = 0; i < 3; i ++) {
        for (auto j = 0; j < 3; j ++)
            cout << setprecision(10) << test_K[i][j] << "\t" ;
        cout << endl;
    }

    double * test_F = new double[3];
    test_triangle.Column_F(test_F ,150);

    cout << "F = :";
    for (auto i = 0; i < 3; i++)
        cout << test_F[i] << "\t";*/
    /*****************************************************************************************************/


    /*****************************************************************************************************/
    /* Объявление переменных */

    /**
     * Debug сообщения
     */
    boolean debug = false;
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
     * Число точек
     */
    uint_fast32_t n = triangle.GetGreed(debug) +1;
    /**
     * Число узлов
     */
    uint_fast32_t k = triangle.triangles_array.size();
    /**
     * Тепловой поток, заданный на границе
     */
    double q = -func_calculate_q(TAU);
    //q = 0;
    /**
     * Столбец правых частей для элемента
     */
    double *F = new double[DIMENSION];
    /**
     * Значения температуры в узлах
     */
    double *Temp = new double[n];
    /**
     * Cюда записываем те температуры с предыдыщуего слоя которые нас интересуют для данного конечного элемента
     */
    double Fi[3];
    /**
     * Промежуточная переменная для записи произведения матрицы и столбца
     */
    double *Resultic = new double[DIMENSION];
    /**
     * Вектор правых частей итоговой системы
     */
    double *Result = new double[n];
    /**
     * Итоговая матрица жесткости системы
     */
    double **Main_Matrix= new double *[n];
    /**
     * Матрица коэффициентов элемента K
     */
    double **K= new double *[DIMENSION];
    /**
     * Матрица коэффициентов элемента C
     */
    double **C = new double *[DIMENSION];
    /**
     * Массив номеров строк в итоговом векторе правых частей Result,
     * к которому нужно будет прибавить полученные правые части F на элементе
     */
    uint_fast32_t ind[DIMENSION];
    /**
     * Результат вычитания двух матриц. Вспомогательная переменная
     */
    double ** substracted_matrix = new double* [DIMENSION];
    /**
     * Еще один Гаусс
     */
    vector< vector<double> > equation1_A(n, std::vector<double>(n+1));
    gauss_12062017 equation1;
    /*****************************************************************************************************/



    /* Debug */
    LOG(INFO) << "Number of dots = " << n;
    LOG(INFO) << "Number of triangles = " << k;
    LOG(INFO) << "Heat flux q = " << q;
    LOG(INFO) << "Space step h = " << STEP_X;
    LOG(INFO) << "Time step tau = " << TAU;
    cout << endl;
    /* End debug */

    /* Задаем начальное условие распределения температуры */
    for (uint_fast32_t i = 0; i < n; i++) {
        Temp[i] = INITIAL_TEMPERATURE;
    }

    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        substracted_matrix[i] = new double[DIMENSION];
    }

    /* Инициализация матрц K, C и вектора правых частей F для элемента */
    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        K[i] = new double[DIMENSION];
        C[i] = new double[DIMENSION];
    }

    /* Ищем решения по временным слоям */
    for (uint_fast32_t global_tau =0; global_tau < 5; global_tau++) {

        for (uint_fast32_t i = 0; i < DIMENSION; i++) {
            F[i] = 0;
            Resultic[i] = 0;
            Fi[i] = 0;
            for (uint_fast32_t j = 0; j < DIMENSION; j++) {
                K[i][j] = 0;
                C[i][j] = 0;
                substracted_matrix[i][j] = 0;
            }
        }
        /* Debug */
        cout << "Temperature on step " << global_tau << endl;
        for (uint_fast32_t i = 0; i < n; i++) {
            cout << Temp[i] << "\t";
        }
        cout << endl;
        /* End debug */

        /* Инициализация матрицы жесткости и вектора правой части итоговой системы */
        for (uint_fast32_t i = 0; i < n; i++) {
            Main_Matrix[i] = new double[n];
            Result[i] = 0;
            for (int j = 0; j < n; j++) {
                Main_Matrix[i][j] = 0;
            }
        }


        double **R1; R1 = new double *[n];
        double *FF; FF  = new double [n];
        for (int kk = 0; kk<n; kk++){
            R1[kk] = new double[n];
            FF[kk] = 0;
            for (int j = 0; j< n; j++)
                R1[kk][j] = 0;
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
            for (uint_fast32_t j = 0; j < DIMENSION; j++) {
                Fi[j] = Temp[ind[j]];
                cout<<Fi[j]<<" ";
            }
            cout<<endl;

            func_substract_two_matrices(substracted_matrix, C, K, 3, 2.0/TAU);


           /* cout<<"Fi"<<endl;
            for (int ii = 0; ii<3; ii++){
                cout<<Fi[ii]<<" ";
            }
            cout<<endl;

            cout<<"substracted_matrix"<<endl;
            for (int ii = 0; ii<3; ii++){
                for (int j = 0; j<3; j++)
                cout<<substracted_matrix[ii][j]<<" ";
            }
            cout<<endl;*/

            for (uint_fast32_t ii = 0; ii < DIMENSION; ii++)
                Resultic[ii] = 0;
            for (int ix = 0; ix < DIMENSION; ix++) {
                Resultic[ix] = 0;
                for (int jx = 0; jx < DIMENSION; jx++)
                    Resultic[ix] += substracted_matrix[ix][jx] * Fi[jx];
            }

            //func_multiply_matrix_and_vector(Resultic, substracted_matrix, Fi, 3);
            /*cout<<"Resultic"<<endl;
            for (int ii = 0; ii<3; ii++){
                cout<<Resultic[ii]<<" ";
            }
            cout<<endl;*/



            /* Заполняем итоговую матрицу жесткости и вектор правой части  полученными на элементе значениями */
            cout<<"F[j]"<<endl;
            for (uint_fast32_t j = 0; j < 3; j++) {
                //Result[ind[j]] += (Resultic[j]) - 2 * F[j];
                FF[ind[j]] += - 2 * F[j];
                //cout<<"FF[ind[j]] = "<<FF[ind[j]]<<" ";
               //cout<<"F[j] = "<<F[j]<<" ";
                //Result[ind[j]] += F[j];
                for (uint_fast32_t l = 0; l < 3; l++) {
                    Main_Matrix[ind[j]][ind[l]] += (C[j][l] * 2 / TAU) + K[j][l];
                    R1[ind[j]][ind[l]] += substracted_matrix[j][l];
                    //Main_Matrix[ind[j]][ind[l]] += K[j][l];
                }
            }
            cout<<endl;


            /* Debug*/
            /*if ( i == (k - 2 ) || i == (k - 1)) {
                cout << "C" << endl;
                for (auto m = 0; m < 3; m++)
                    for (auto p = 0; p < 3; p++)
                        cout << C[m][p] << "\t";
                cout << endl;

                cout << "K" << endl;
                for (auto m = 0; m < 3; m++)
                    for (auto p = 0; p < 3; p++)
                        cout << K[m][p] << "\t";
                cout << endl;

                cout << "Resultic" << endl;
                for (auto m = 0; m < 3; m++)
                    cout << Resultic[m] << "\t";
                cout << endl;

                cout << "F" << endl;
                for (auto m = 0; m < 3; m++)
                    cout << F[m] << "\t";
                cout << endl;
            }*/
            /* End debug*/
        }
        for (int ix = 0; ix < n; ix++) {
            Result[ix] = 0;
            for (int jx = 0; jx <n; jx++)
                Result[ix] += R1[ix][jx] * Temp[jx];
        }
        cout<<"Result 1"<<endl;
        for (int ii = 0; ii<n; ii++){
            cout<<Result[ii]<<" ";
            Result[ii]+= FF[ii];
        }
        cout<<endl;
        cout<<"Result 2"<<endl;
        for (int ii = 0; ii<n; ii++){
            cout<<Result[ii]<<" ";
        }
        cout<<endl;
        /* Debug */
        for (int p = 0; p < n; p++) {
            for (int j = 0; j < n; j++) {
                file1 << Main_Matrix[p][j] << " ";
            }
            file1 << endl;
        }
        file1 << " ------------------------------------------" << endl;
        file1 << "Result" << endl;
        for (int ii = 0; ii < n; ii++)
            file1 << Result[ii] << " " << endl;
        file1 << " --------********---------------------------" << endl;
        /* End debug */


        for (int p=0; p<n; p++) {
            for (int j=0; j<n; j++) {
                equation1_A[p][j] = Main_Matrix[p][j];
            }
        }
        for (int p=0; p<n; p++) {
            equation1_A[p][n] = Result[p];
        }

        //print(equation1_A);


        equation1.gauss(Temp, equation1_A);
        //equation.GaussSolve(Temp, Main_Matrix, Result, n);

        /* Debug */
        cout << "After GaussSolve: " << endl;
        cout << "Temperature " << "\t" << "Result" << endl;
        for (auto p = 0; p < n; p++) {
            cout << FF[p] <<"\t"<<Result[p]<<endl;
            //cout << Result[i] << endl;
        }
        cout << "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << endl;
        /* End debug */

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
    file1.close();
    file2.close();
    return 0;
}

double func_calculate_rho(double h)
{
    return (1.2 * exp(-0.00013 * h));
}

double func_calculate_q(double t)
{
    // Тепловой поток при скорости 3000 м/c на расстоянии 1 м от Земли
    double S = 1;
    double rho = func_calculate_rho(1);
    double V = 3000;
    double pi = 3.1415926;
    double r = 1;
    double H = 2;

    return ((S * rho * V * V * V) / (Thermal_Conductivity * pi * r * sqrt(r * r + H * H)));
}

void func_multiply_matrix_and_vector(double* result_vector, double** matrix, double* vect, uint_fast32_t dim) {
    for (uint_fast32_t i = 0; i < dim; i++)
        result_vector[i] = 0;

    for (uint_fast32_t i = 0; i < dim; i++) {
        for (uint_fast32_t j = 0; j < dim; j++) {
            result_vector[i] += vect[j] * matrix[j][i];
        }
    }
}

void func_substract_two_matrices(double ** result_matrix, double** matrix1, double** matrix2, uint_fast32_t dim, double devide_element) {

    for (auto i = 0; i < dim; i++)
        for (auto j = 0; j < dim; j++)
            result_matrix[i][j] = (devide_element * matrix1[i][j]) - matrix2[i][j];

}