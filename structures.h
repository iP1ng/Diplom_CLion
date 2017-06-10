
#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <iostream>

const double step_x = 0.5;
const double Pi = 3.1415926;

/*
 * Структура, описывающая координаты и номер узла, а также радиус-вектор узла.
 */
struct points {
    double_t x;
    double_t y;
    double_t rad_vector() { return x;  }
    uint_fast32_t point_num;
};

/*
 * Структура, описывающая координаты точек и номер треугольника,
 * а также вычисляющая коэффициенты функции формы узлов, входящих в треугольник,
 * а также вычисляюая площадь треугольника и матрицы элементов [K], [C], [F].
 */
struct triangles {
    points first_point;
    points second_point;
    points third_point;

    void coef_a(double_t* a)
    {
        a[0] = second_point.x * third_point.y - third_point.x * second_point.y;
        a[1] = third_point.x * first_point.y - first_point.x * third_point.y;
        a[2] = first_point.x * second_point.y - second_point.x * first_point.y;
    }
    void coef_b(double_t* b)
    {
        b[0] = second_point.y - third_point.y;
        b[1] = third_point.y - first_point.y;
        b[2] = first_point.y - second_point.y;
    }
    void coef_c(double_t* c)
    {
        c[0] = third_point.x - second_point.x;
        c[1] = first_point.x - third_point.x;
        c[2] = second_point.x - first_point.x;
    }
    double_t GetSquareTriangleArea()
    {
        return step_x * 2 * step_x * 0.5;
    }
    double_t GetMatrixADeterminant()
    {
        return 0.5 * (second_point.x * third_point.y
                - third_point.x * second_point.y
                - first_point.x * third_point.y
                + first_point.x * second_point.y
                + third_point.x * first_point.y
                - second_point.x * first_point.y);
    }
    void Matrix_K(double** K)
    {
        points p[3];
        double_t a[3];
        double_t b[3];
        double_t c[3];
        coef_a(a);
        coef_b(b);
        coef_c(c);

        p[0].x = first_point.x;
        p[0].y = first_point.y;
        p[1].x = second_point.x;
        p[1].y = second_point.y;
        p[2].x = third_point.x;
        p[2].y = third_point.y;
        double_t R = (1/12) * ((2 * p[0].rad_vector()
                     + p[1].rad_vector()
                     + p[2].rad_vector()) * p[0].rad_vector()
                     + (p[0].rad_vector()
                     + 2 * p[1].rad_vector()
                     + p[2].rad_vector()) * p[1].rad_vector()
                     + (p[0].rad_vector()
                     + p[1].rad_vector()
                     + 2 * p[2].rad_vector()) * p[2].rad_vector());
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                K[i][j] = (2 * Pi * R / (4 * this->GetMatrixADeterminant())) * (b[i] * b[j] + c[i] * c[j]);
            }
        }
    }
    void Matrix_C(double** C)
    {
        points p[3];
        p[0].x = first_point.x;
        p[0].y = first_point.y;
        p[1].x = second_point.x;
        p[1].y = second_point.y;
        p[2].x = third_point.x;
        p[2].y = third_point.y;

        double_t D = 2 * Pi * this->GetMatrixADeterminant() / 180;
        double_t R[3];
        for (int i = 0; i < 3; i++) {
            R[i] = p[i].rad_vector();
        }

        C[0][0] = D * (12 * R[0] * R[0]
                       + 2 * R[1] * R[1]
                       + 2 * R[2] * R[2]
                       + 6 * R[0] * R[1]
                       + 6 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[0][1] = D * (3 * R[0] * R[0]
                       + 3 * R[1] * R[1]
                       + R[2] * R[2]
                       + 4 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[0][2] = D * (3 * R[0] * R[0]
                       + 1 * R[1] * R[1]
                       + 3 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 4 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[1][0] = C[0][1];
        C[2][0] = C[0][2];
        C[1][1] = D * (2 * R[0] * R[0]
                       + 12 * R[1] * R[1]
                       + 2 * R[2] * R[2]
                       + 6 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 6 * R[1] * R[2]);
        C[1][2] = D * (1 * R[0] * R[0]
                       + 3 * R[1] * R[1]
                       + 3 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 4 * R[1] * R[2]);
        C[2][2] = D * (2 * R[0] * R[0]
                       + 2 * R[1] * R[1]
                       + 12 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 6 * R[0] * R[2]
                       + 6 * R[1] * R[2]);
        C[2][1] = C[1][2];
    }
    void Column_F(double* F, double q)
    {
        double L = sqrt(step_x * step_x + 4 * step_x * step_x);
        double k = L * q * 2 * Pi / 6;
        F[0] = 0;
        F[1] = 0;
        F[2] = 0;
        if (fabs(2 * first_point.x - first_point.y) < 0.01 && fabs(2 * third_point.x - third_point.y) < 0.01) {
            F[0] = k * (2 * first_point.rad_vector() + third_point.rad_vector());
            F[2] = k * (first_point.rad_vector() + 2 * third_point.rad_vector());
            //std::cout << "Hi!" << std::endl;
        }
    }
};

#endif //STRUCTURES_H
