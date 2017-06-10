
#ifndef DIPLOM_CONSTANTS_H
#define DIPLOM_CONSTANTS_H

const int TRIANGLE_ARRAY_DIMENSION = 300;
const int GRID_ARRAY_DIMENSION = 300;
const points A = (points){ 0, 0 };
const points B = (points){ 1, 0 }; // осталось с прошлой версии проги. Переделать, сейчас по факту тут в два раза больше
// указывается в x, т.к. раньше треугольник зеркально отображался
const points C = (points){ 1, 2 };
const double_t STEP_X = 0.5;
const double_t TRIANGLE_HEIGHT = 2;
const double_t EPS_T = 0.01;
const double_t TAU = 0.01;

#endif //DIPLOM_CONSTANTS_H
