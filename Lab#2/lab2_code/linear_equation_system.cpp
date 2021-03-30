#include "linear_equation_system.h"
#include <iostream>

Vector GaussianElimination(Matrix A, Vector b){
    const int n = A.size();
    bool have_res = true;
    for(int i = 0; i < n; ++i){
        // choose the largest element
        int max_idx = i;
        for(int j = i; j < n; ++j){
            if(abs(A[j][i]) > abs(A[max_idx][i])){
                max_idx = j;
            }
        }
        // check for degeneracy matrix
        if(A[max_idx][i] == 0){
            have_res = false;
            break;
        }
        // swap i-th and j-th rows
        if(max_idx != i){
            std::swap(A[i], A[max_idx]);
            std::swap(b[i], b[max_idx]);
        }
        // divide i-th row by A[i][i]
        for(int j = i+1; j < n; ++j){
            A[i][j] /= A[i][i];
        }
        b[i] /= A[i][i];
        A[i][i] = 1;
        // Subtract i-th row from j-th row (j = 1, 2, ..., n and j != i)
        for(int j = 0; j < n; ++j){
            if(j != i) {
                for (int k = i + 1; k < n; ++k) {
                    A[j][k] -= A[i][k] * A[j][i];
                }
                b[j] -= b[i] * A[j][i];
                A[j][i] = 0;
            }
        }
    }
    return have_res ? b : Vector{};
}

Vector TridiagonalMatrixAlgorithm(Matrix A, Vector d){
    const int n = A.size();
    if(n <= 2){
        return n != 0 ? d : Vector{};
    }else{
        A[0][1] /= A[0][0]; //c'_1 = c_1/b_1
        d[0] /= A[0][0]; //d'_1 = d_1/b_1
        for(int i = 1; i < n-1; ++i){
            // c'_i = c_i/(b_i - c'_{i-1}*a_i)
            A[i][i+1] /= A[i][i] - A[i][i-1]*A[i-1][i];
        }
        for(int i = 1; i < n; ++i){
            // c'_i = (d_i - d'_{i-1}*a_i)/(b_i - c'_{i-1}*a_i)
            d[i] = (d[i]-A[i][i-1]*d[i-1])/(A[i][i] - A[i][i-1]*A[i-1][i]);
        }

        Vector x(n);
        x[n-1] = d[n-1]; // x_n = d'_n
        for(int i = n-2; i >= 0; --i){
            // x_i = d'_i - c'_i*x_{i+1}
            x[i] = d[i] - A[i][i+1]*x[i+1];
        }
        return x;
    }
}