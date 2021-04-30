#include <iostream>
#include "linear_equation_system.h"

using std::cout;

void PrintLinearEquationsSystem(const Matrix &A, const Vector &b){
    const int n = A.size();
    for(int i = 0; i < n; ++i){
        bool is_first_val = true;
        for(int j = 0; j < n; ++j){
            if(is_first_val){
                if(A[i][j] != 0) {
                    cout << A[i][j] << "x_" << j+1;
                    is_first_val = false;
                }
            }else {
                if (A[i][j] < 0) {
                    cout << " - ";
                    if(A[i][j] != -1){
                        cout << -1*A[i][j];
                    }
                    cout << "x_" << j+1;
                } else if (A[i][j] > 0) {
                    cout << " + ";
                    if(A[i][j] != 1){
                        cout << A[i][j];
                    }
                    cout << "x_" << j+1;
                }
            }
        }
        if(is_first_val){
            cout << 0;
        }
        cout << " = " << b[i] << std::endl;
    }
}

void PrintResult(const Vector &x){
    if(!x.empty()) {
        for (int i = 0; i < x.size(); ++i) {
            cout << "x_" << i+1 << " = " << x[i];
            cout << (i != x.size() - 1 ? " , " : ".\n");
        }
    }else{
        cout << "Have no result!";
    }
}

void Task1(){
    std::cout << "------------------------------\n";
    std::cout << "Jacobi method\n";
    // Ax = b
    const Matrix A {{4, 0, 1, 1},
                    {0, 3, 0, 1},
                    {1, 0, 2, 0},
                    {1, 1, 0, 5}};
    const Vector b {11, 10, 7, 23};
    PrintLinearEquationsSystem(A, b);
    Vector x = JacobiMethod(A, b);
    cout << "Result:\n";
    PrintResult(x);
}

void Task2(){
    std::cout << "------------------------------\n";
    std::cout << "Seidel method\n";
    // Ax = b
    const Matrix A {{3, 0, 0, 1},
                    {0, 6, 2, 0},
                    {0, 2, 3, 0},
                    {1, 0, 0, 4}};
    const Vector b {7, 18, 13, 17};
    PrintLinearEquationsSystem(A, b);
    Vector x = SeidelMethod(A, b);
    cout << "Result:\n";
    PrintResult(x);
}

int main() {
    Task1();
    Task2();
    return 0;
}