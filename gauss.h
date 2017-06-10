
#ifndef DIPLOM_GAUSS_H
#define DIPLOM_GAUSS_H
//#include <cstdlib>
//Если использовать cstdlib std:abs, то в ответе нет чисел с отрицательным знаком и в целом ответ полностью меняется.
// статья на тему https://stackoverflow.com/questions/21392627/abs-vs-stdabs-what-does-the-reference-say
// https://stackoverflow.com/questions/3118165/when-do-i-use-fabs-and-when-is-it-sufficient-to-use-stdabs
// <cmath> и std::abs дает такие же результаты
#include <math.h>

class GaussMethod {
public:
    double* GaussSolve(double** a, double* y, int n);
};


#endif //DIPLOM_GAUSS_H
