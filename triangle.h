
#ifndef TRIANGLE_H
#define TRIANGLE_H

/* Подключаем быстрые int фиксированного размера */
#include <cstdint>
#include <array>

/* Подключаем double_t тип для оптимальной работы на любой архитектуре компьютера. */
#include <math.h>

/* Мои хедеры */
#include "constants.h"
#include "structures.h"


/* Debug */
#include <iostream>
#include <fstream>
#include "EasyLogging/easylogging++.h"

/*
 * Построение треугольной сетки на равнобедренном треугольнике.
 * Левая координата треугольника должна совпадать с началом координат.
 */
class IsoscelesTriangleGrid {

private:
    /*
     * m_a - левая точка треугольника
     * m_b - правая точка треугольника
     * m_c - центральная точка треугольника
     * m_hx - шаг по оси x
     */
    points m_a;
    points m_b;
    points m_c;
    double_t m_hx;
    triangles triangle;

public:
    std::vector<triangles> triangles_array;

    /*
     * Конструктор класса IsoscelesTriangleGrid
     */
    IsoscelesTriangleGrid(points a, points b, points c, double_t hx)
    {
        m_a = a;
        m_b = b;
        m_c = c;
        m_hx = hx;
    }

    double_t GetStep() { return m_hx; }

    /*
     * Метод построения сетки
     * TODO массив задавать не хардкодом, а вычислять число треугольников заранее
     * На выход число узлов
     */
    uint_fast32_t GetGreed();

    /*
     * Функция, вычисляющая значение координаты y по заданной координате x на левой грани треугольника
     * или на прямой, сдвинутой относительно левой грани на step
     */
    double_t LineFunction_ab(double_t x, double_t step);

    /*
     * Функция, вычисляющая значение координаты y по заданной координате x на правой грани треугольника
     * (граница фиксированная)
     */
    double_t LineFunction_bc(double_t x);
};

#endif //TRIANGLE_H
