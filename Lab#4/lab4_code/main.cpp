#include <iostream>
#include <iomanip>
#include <cmath>
#include <exception>
#include "mathlib.h"

using std::cout;

void Task1 () {
    std::cout << "------------------------------\n";
    std::cout << "Trapezoidal rule:\n";
    mathlib::Function f = [] (long double x) { return x*x*x*x + 2*x*x + x; };
    long double result = mathlib::TrapezoidalRule(f, 0, 2.5, 321);
    cout << "Integral of x^4 + 2*x^2 + x from 0 to 2.5 equal " << result << "\n";
}

void Task2 () {
    cout << "------------------------------\n";
    cout << "Solve equation system {sin(x) + 2*y - 1.6 = 0, cos(y - 1) - 1 = 0}" << "\n";
    mathlib::EquationSystem F{
        [] (const mathlib::Vector& x) {
            if (x.size() != 2) { throw std::invalid_argument("x.size() != 2"); }
            return std::sin(x[0]) + 2*x[1] - 1.6l;
        },
        [] (const mathlib::Vector& x) {
            if (x.size() != 2) { throw std::invalid_argument("x.size() != 2"); }
            return std::cos(x[1] - 1.l) - 1.l;
        }};
    mathlib::JacobianMatrix A {
            {[] (const mathlib::Vector& x) { return std::cos(x[0]); }, [] (const mathlib::Vector& x) { return 2; }},
            {[] (const mathlib::Vector& x) { return 0; }, [] (const mathlib::Vector& x) { return -std::sin(x[1]-1); }}
    };
    const mathlib::Vector res = mathlib::NewtonMethod(F, A, {1, 0}, 1e-6);
    cout << "x = " << res[0] << "\n";
    cout << "y = " << res[1] << "\n";
}

void Task3 () {
    std::cout << "------------------------------\n";
    std::cout << "Maximal eigenvalues:\n";
    mathlib::Matrix A = {{1, 2, 3},
                         {2, 3, 4},
                         {3, 4, 5}};
    long double result = mathlib::MaximalEigenvalues(A, {1, 1, 1}, 0.0001);
    cout << A;
    cout << "Maximal eigenvalues equal " << std::fixed << std::setprecision(5) << result << "\n";
}

int main() {
    Task1();
    Task2();
    Task3();
    return 0;
}