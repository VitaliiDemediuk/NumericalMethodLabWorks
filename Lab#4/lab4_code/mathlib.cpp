#include "mathlib.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <exception>

namespace {
    using namespace mathlib;

    long double Norm(const mathlib::Vector& x) {
        long double sum = std::accumulate(x.begin(), x.end(), 0.l,
                                          [] (long double part_sum, long double a) { return part_sum + a*a; });
        return std::sqrt(sum);
    }

    Matrix JacobianMatrixOfVector (const JacobianMatrix& A, const Vector& v) {
        const size_t n = A.size();
        if (n != v.size()) {
            throw std::invalid_argument("A.size() != v.size()");
        }

        Matrix res(n, Vector(n));
        for (size_t i = 0; i < n; ++i) {
            if (n != A[i].size()) {
                throw std::invalid_argument("A[i].size() != v.size()");
            }
            for(size_t j = 0; j < n; ++j){
                res[i][j] = A[i][j](v);
            }
        }

        return res;
    }

    Vector EquationSystemOfVector (const EquationSystem& F, const Vector& v) {
        const size_t n = F.size();
        if (n != v.size()) {
            throw std::invalid_argument("F.size() != v.size()");
        }

        Vector res(n);
        for (size_t i = 0; i < n; ++i) {
            res[i] = F[i](v);
        }

        return res;
    }

    Vector &operator*=(Vector &x, long double a) {
        for_each(begin(x), end(x), [a] (long double &el) { el *= a; });
        return x;
    }

    Vector operator*(const Vector &x, long double a) {
        Vector res = x;
        res *= a;
        return res;
    }

    long double operator*(const Vector &x1, const Vector &x2) {
        if (x1.size() != x2.size()) {
            throw std::invalid_argument("x1.size() != x2.size()");
        }

        long double res = 0;
        for (size_t i = 0; i < x1.size(); ++i) {
            res += x1[i] * x2[i];
        }
        return res;
    }

    Vector& operator-=(Vector& x1, const Vector& x2) {
        if (x1.size() != x2.size()) {
            throw std::invalid_argument("x1.size() != x2.size()");
        }

        for (size_t i = 0; i < x1.size(); ++i) {
            x1[i] -= x2[i];
        }
    }

    Vector operator*(const Matrix &A, const Vector &x) {
        if (A.size() != x.size()) {
            throw std::invalid_argument("A.size() != x.size()");
        }

        Vector res(A.size(), 0);

        for (int i = 0; i < A.size(); ++i) {
            res[i] = A[i] * x;
        }

        return res;
    }

}

namespace mathlib {
    using namespace std;

    Vector GaussianElimination (Matrix A, Vector b) {
        const size_t n = A.size();
        if (n != b.size()){
            throw std::invalid_argument("A.size() != b.size()");
        }

        for (size_t i = 0; i < n; ++i) {
            // choose the largest element
            int max_idx = i;
            for (size_t j = i; j < n; ++j) {
                if (abs(A[j][i]) > abs(A[max_idx][i])) {
                    max_idx = j;
                }
            }
            // check for degeneracy matrix
            if (A[max_idx][i] == 0) {
                throw std::invalid_argument("Matrix A is degenerate");
            }
            // swap i-th and j-th rows
            if (max_idx != i) {
                std::swap(A[i], A[max_idx]);
                std::swap(b[i], b[max_idx]);
            }
            // divide i-th row by A[i][i]
            for (size_t j = i+1; j < n; ++j) {
                A[i][j] /= A[i][i];
            }
            b[i] /= A[i][i];
            A[i][i] = 1;
            // Subtract i-th row from j-th row (j = 1, 2, ..., n and j != i)
            for (size_t j = 0; j < n; ++j) {
                if (j != i) {
                    for (size_t k = i + 1; k < n; ++k) {
                        A[j][k] -= A[i][k] * A[j][i];
                    }
                    b[j] -= b[i] * A[j][i];
                    A[j][i] = 0;
                }
            }
        }
        return b;
    }

    long double TrapezoidalRule(const Function &f, long double a, long double b, int n) {
        long double d = b - a; //delta
        long double h = d / n;
        long double res = f(a) * h / 2;
        long double x = a;
        for (int i = 1; i < n; ++i) {
            x += h;
            res += f(x) * h;
        }
        res += f(b) * h / 2;
        return res;
    }

    long double MaximalEigenvalues(const Matrix &A, Vector x_k, long double eps) {
        long double u_k, u_k_next;
        u_k_next = MAXFLOAT;
        int i = 0;
        do {
            ++i;
            Vector e_k = x_k * (1 / Norm(x_k));
            x_k = A * e_k;
            u_k = u_k_next;
            u_k_next = x_k * e_k;
        } while (std::abs(u_k_next - u_k) > eps);
        return u_k_next;
    }

    Vector NewtonMethod(const EquationSystem& F, const JacobianMatrix& A, const Vector& x, long double eps){
        Vector z_k;
        Vector x_k = x;

        do {
            auto A_k = JacobianMatrixOfVector(A, x_k);
            auto F_k = EquationSystemOfVector(F, x_k);
            z_k = GaussianElimination(A_k, F_k);
            x_k -= z_k;
        } while (Norm(z_k) > eps);

        return x_k;
    }

}

std::ostream &operator<<(std::ostream &out, const mathlib::Matrix &A) {
    for (const auto &v : A) {
        for (auto el : v) {
            out << el << ' ';
        }
        out << std::endl;
    }
    return out;
}