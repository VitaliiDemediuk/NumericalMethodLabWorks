#ifndef LAB1_CODE_EQUATION_H
#define LAB1_CODE_EQUATION_H

#include <vector>

using Matrix = std::vector<std::vector<long double>>;
using Vector = std::vector<long double>;

//Ax = b
Vector GaussianElimination(Matrix A, Vector b);
Vector TridiagonalMatrixAlgorithm(Matrix A, Vector b);

#endif //LAB1_CODE_EQUATION_H