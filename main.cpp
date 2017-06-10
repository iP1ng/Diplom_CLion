
#include "triangle.h"
#include "gauss.h"

INITIALIZE_EASYLOGGINGPP

double_t func_calculate_rho(double_t h);
double_t func_calculate_q(double_t t);

using namespace std;
int main()
{
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, "%level: %msg");
    LOG(INFO) << "Starting programm...";

    GaussMethod equation;

    ofstream file1("Main_matrix.dat");

    const points Point_A = (points){ 0, 0 };
    const points Point_B = (points){ 1, 0 };
    const points Point_C = (points){ 1, 2 };

    IsoscelesTriangleGrid triangle(Point_A, Point_B, Point_C, STEP_X);
    double tau = TAU;
    // chislo tochek
    int n = triangle.GetGreed() + 1;
    LOG(INFO) << "Number of dots = " << n;
    // chislo treugolnilov
    int k = triangle.triangles_array.size();
    LOG(INFO) << "Number of triangles = " << k;

    double q = func_calculate_q(tau);
    double* Result;
    Result = new double[n];
    double* F;
    F = new double[3];
    double* Temp;
    Temp = new double[n];

    double** Main_Matrix;
    Main_Matrix = new double*[n];

    for (int i = 0; i < n; i++) {
        Main_Matrix[i] = new double[n];

        for (int j = 0; j < n; j++) {
            Main_Matrix[i][j] = 0;
        }
    }
    for (int i =0; i<n; i++){
        Result[i] = 0;
        Temp[i] = INITIAL_TEMPERATURE;

    }
    double_t** K;
    K = new double_t*[3];
    for (int i = 0; i < 3; i++) {
        K[i] = new double_t[3];
        F[i] = 0;
    }
    double_t** C;
    C = new double_t*[3];
    for (int i = 0; i < 3; i++) {
        C[i] = new double_t[3];
    }

    for (int i = 0; i < k; i++) {
        triangle.triangles_array[i].Matrix_K(K);
        triangle.triangles_array[i].Matrix_C(C);
        triangle.triangles_array[i].Column_F(F, q);

        cout << "Matrix K: " << endl;

        for (int iii = 0; iii < 3; iii++) {
            for (int j = 0; j < 3; j++) {
                cout << K[iii][j] << " ";
            }
            cout << endl;
        }

        int ind[3];
        ind[0] = triangle.triangles_array[i].first_point.point_num;
        ind[1] = triangle.triangles_array[i].second_point.point_num;
        ind[2] = triangle.triangles_array[i].third_point.point_num;

        cout<<"-------------"<<endl;
        for (int ii = 0; ii < 3; ii++) {
            Result[ind[ii]] += -F[ii];
            for (int j = 0; j < 3; j++){
                Main_Matrix[ind[ii]][ind[j]] += C[ii][j] / tau + 0.5 * K[ii][j];
                cout<<ind[ii]<<" "<<ind[j]<<endl;
                //Main_Matrix[ind[ii]][ind[j]] += 0.5 * K[ii][j];
            }
        }

        file1<< k << endl;
        for (int iii = 0; iii < n; iii++) {
            for (int j = 0; j < n; j++) {
                file1 << Main_Matrix[iii][j] << " ";
            }
            file1 << endl;
        }

        file1<<" ------------------------------------------"<< endl;
    }


    for (int jj = 0; jj < n; jj++){
        Result[jj] += Temp[jj] / tau;
        cout<<Result[jj]<<endl;

    }



    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file1 << Main_Matrix[i][j] << " ";
        }
        file1 << endl;
    }
    file1.close();
    Temp = equation.GaussSolve(Main_Matrix, Result, n);
    for (int i = 0; i < n; i++)
        cout << Temp[i] << endl;
    return 0;
}

double_t func_calculate_rho(double_t h){
    return  (1.2 * exp(- 0.00013 * h));
}

double_t func_calculate_q(double_t t){
    // Тепловой поток при скорости 3000 м/c на расстоянии 1 м от Земли
    double_t S = 1;
    double_t rho = func_calculate_rho(1);
    double_t V = 3000;
    double_t pi = 3.1415926;
    double_t r = 1;
    double_t H = 2;

    return ((S * rho * V * V * V) / (Thermal_Conductivity * pi * r * sqrt(r * r + H * H)));
}
