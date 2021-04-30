#include "linear_equation_system.h"
#include <cmath>

Vector operator- (const Vector& v1, const Vector& v2){
    Vector v3;
    if (v1.size() == v2.size()) {
        const size_t n = v1.size();
        v3.resize(n);
        for(int i = 0; i < n; ++i){
            v3[i] = v1[i] - v2[i];
        }
    }
    return v3;
}

static long double Norma(const Vector& v){
    long double norma = 0;
    for(auto x : v){
        norma = std::max(norma, std::abs(x));
    }
    return norma;
}

Vector JacobiMethod(const Matrix& A, const Vector& b){
    const size_t n = A.size();
    Vector prev_res(n);
    Vector res(n, 0);
    int k = 0;
    do{
        std::swap(prev_res, res);
        for(int i = 0; i < n; ++i){
            res[i] = 0;
            for(int j = 0; j < i; ++j){
                res[i] -= (A[i][j]*prev_res[j])/A[i][i];
            }
            for(int j = i+1; j < n; ++j){
                res[i] -= (A[i][j]*prev_res[j])/A[i][i];
            }
            res[i] += b[i]/A[i][i];
        }
        ++k;
    } while(Norma(res - prev_res) > 1e-9);
    return res;
}

Vector SeidelMethod(const Matrix& A, const Vector& b){
    const size_t n = A.size();
    Vector prev_res(n);
    Vector res(n, 0);
    do {
        std::swap(prev_res, res);
        for(int i = 0; i < n; ++i){
            res[i] = 0;
            for(int j = 0; j < i; ++j){
                res[i] -= (A[i][j]*res[j])/A[i][i];
            }
            for(int j = i+1; j < n; ++j){
                res[i] -= (A[i][j]*prev_res[j])/A[i][i];
            }
            res[i] += b[i]/A[i][i];
        }
    } while(Norma(res - prev_res) > 1e-9);
    return res;
}