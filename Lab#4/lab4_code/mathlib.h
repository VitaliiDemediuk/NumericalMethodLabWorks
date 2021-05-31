#ifndef LAB4_CODE_MATHLIB_H
#define LAB4_CODE_MATHLIB_H

#include <functional>
#include <ostream>

namespace mathlib{

    using Function = std::function<long double (long double)>;
    using Matrix = std::vector<std::vector<long double>>;
    using Vector = std::vector<long double>;
    using VectorFunction = std::function<long double (const Vector&)>;
    using EquationSystem = std::vector<VectorFunction>;
    using JacobianMatrix = std::vector<std::vector<VectorFunction>>;

    // integral of f(x)dx from a to b, n - segments
    Vector GaussianElimination(Matrix A, Vector b);
    long double TrapezoidalRule(const Function& f, long double a, long double b, int n);
    long double MaximalEigenvalues(const Matrix& A, Vector x_0, long double eps);
    Vector NewtonMethod(const EquationSystem& F, const JacobianMatrix& A, const Vector& x, long double eps);
}

std::ostream& operator<<(std::ostream& out, const mathlib::Matrix& A);

#endif //LAB4_CODE_MATHLIB_H
