#include <iostream>
#include <vector>
#include "equation.h"

void task1(){
    std::cout << "------------------------------\n";
    std::cout << "Relaxation method" << std::endl;
    std::cout << "x^3 - 6*x^2 + 5*x + 12 = 0" << std::endl;
    std::vector<long double> solutions;
    auto f = [](long double x){ return x*x*x - 6*x*x + 5*x + 12; };

    solutions.push_back(equation::solve_relaxation(f, -0.043, 0, "solve_relaxation.log"));
    std::cout << "x*0 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_relaxation(f, 0.28, 2, "solve_relaxation.log"));
    std::cout << "x*1 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_relaxation(f, -0.097, 5, "solve_relaxation.log"));
    std::cout << "x*2 = " << solutions.back() << std::endl;

    long double min_solution = *std::min_element(solutions.begin(), solutions.end());
    if(min_solution < 0){
        std::cout << "Min negative solution: " << min_solution << std::endl;
    }else{
        std::cout << "There are no negative solutions" << std::endl;
    }
}

void task2(){
    std::cout << "------------------------------\n";
    std::cout << "Newton method" << std::endl;
    std::cout << "x^3 + 3x^2 - x - 3 = 0" << std::endl;
    std::vector<long double> solutions;
    auto f = [](long double x){ return x*x*x + 3*x*x - x - 3; };
    auto df = [](long double x){ return 3*x*x + 6*x - 1; };


    solutions.push_back(equation::solve_newton(f, df, -3.5, "solve_newton.log"));
    std::cout << "x*0 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_newton(f, df, -1.5, "solve_newton.log"));
    std::cout << "x*1 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_newton(f, df, 2.5, "solve_newton.log"));
    std::cout << "x*2 = " << solutions.back() << std::endl;

    long double max_solution = *std::max_element(solutions.begin(), solutions.end());
    if(max_solution > 0){
        std::cout << "Max positive solution: " << max_solution << std::endl;
    }else{
        std::cout << "There are no positive solutions" << std::endl;
    }
}

void task3(){
    std::cout << "------------------------------\n";
    std::cout << "Secant method" << std::endl;
    std::cout << "x^3 + x^2 - 4x - 4 = 0" << std::endl;
    std::vector<long double> solutions;
    auto f = [](long double x){ return x*x*x + x*x - 4*x - 4; };

    solutions.push_back(equation::solve_secant(f, -4,-3.5, "solve_secant.log"));
    std::cout << "x*0 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_secant(f, 0.5, 0, "solve_secant.log"));
    std::cout << "x*1 = " << solutions.back() << std::endl;

    solutions.push_back(equation::solve_secant(f, 1, 1.5, "solve_secant.log"));
    std::cout << "x*2 = " << solutions.back() << std::endl;

    long double max_solution = *std::max_element(solutions.begin(), solutions.end());
    if(max_solution > 0){
        std::cout << "Max positive solution: " << max_solution << std::endl;
    }else{
        std::cout << "There are no positive solutions" << std::endl;
    }

}

int main() {
    //task1();
    //task2();
    task3();
    return 0;
}