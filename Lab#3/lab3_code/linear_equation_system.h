#ifndef LAB1_CODE_EQUATION_H
#define LAB1_CODE_EQUATION_H

#include <vector>

using Matrix = std::vector<std::vector<long double>>;
using Vector = std::vector<long double>;

Vector operator- (const Vector& v1, const Vector& v2);

//Ax = b
Vector JacobiMethod(const Matrix& A, const Vector& b);
Vector SeidelMethod(const Matrix& A, const Vector& b);

#endif //LAB1_CODE_EQUATION_H