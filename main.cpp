


#include "triangle.h"
#include "gauss.h"

// Test
INITIALIZE_EASYLOGGINGPP

using namespace std;
int main()
{
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, "%level  (%func): %msg");
    LOG(INFO) << "Starting programm...";

    GaussMethod equation;
    ofstream file("triangles.dat");
    ofstream file1("Main_matrix.dat");

    IsoscelesTriangleGrid triangle(A, B, C, STEP_X);

    double tau = TAU;

    // chislo tochek
    int n = triangle.GetGreed();
    LOG(INFO) << "Number of dots = " << n;
    // chislo treugolnilov
    int k = triangle.triangles_array.size();
    LOG(INFO) << "Number of triangles = " << k;

    file << k << endl;
    double q = 200;
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
        Temp[i] = 300;

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
        file << triangle.triangles_array[i].first_point.x << " " << triangle.triangles_array[i].first_point.y << " ";
        file << triangle.triangles_array[i].second_point.x << " " << triangle.triangles_array[i].second_point.y << " ";
        file << triangle.triangles_array[i].third_point.x << " " << triangle.triangles_array[i].third_point.y;
        file << endl;
    }
    file.close();

    for (int i = 0; i < k; i++) {
        triangle.triangles_array[i].Matrix_K(K);
        triangle.triangles_array[i].Matrix_C(C);
        triangle.triangles_array[i].Column_F(F, q);

        int ind[3];
        ind[0] = triangle.triangles_array[i].first_point.point_num;
        ind[1] = triangle.triangles_array[i].second_point.point_num;
        ind[2] = triangle.triangles_array[i].third_point.point_num;
        cout<<"-------------"<<endl;
        for (int ii = 0; ii < 3; ii++) {
            Result[ind[ii]] += -F[ii];
            for (int j = 0; j < 3; j++){
                //Main_Matrix[ind[ii]][ind[j]] += C[ii][j] / tau + 0.5 * K[ii][j];
                Main_Matrix[ind[ii]][ind[j]] += 0.5 * K[ii][j];
                cout<<ind[ii]<<" "<<ind[j]<<endl;
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
        //Result[jj] += Temp[jj] / tau;
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
    //Temp = gauss(Main_Matrix, Result, n);
    for (int i = 0; i < n; i++)
        cout << Temp[i] << endl;
    return 0;
}

