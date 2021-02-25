#include "equation.h"
#include <iostream>
#include <fstream>
#include <iomanip>

long double equation::solve_relaxation(const Function& f, long double const_t,
                                       long double x, const std::string& log_file_name, int n){
    auto old_rdbuf = std::clog.rdbuf();
    std::ofstream logfile;
    if(!log_file_name.empty()){
        logfile.open(log_file_name, std::ofstream::app);
        std::clog.rdbuf(logfile.rdbuf());
    }

    std::clog << "------------------------------\n";
    std::clog << "x0 = " << x << std::endl;
    for(int i = 0; i < n and std::abs(f(x)) > 1e-6; ++i){
        x = x+const_t*f(x);
        std::clog << std::fixed << std::setprecision(6) << "x" << i+1 << " = " << x << std::endl;
    }

    if(!log_file_name.empty()){
        std::clog.rdbuf(old_rdbuf);
    }
    return x;
}

long double equation::solve_newton(const Function& f, const Function& df,
                         long double x, const std::string& log_file_name, int n){
    auto old_rdbuf = std::clog.rdbuf();
    std::ofstream logfile;
    if(!log_file_name.empty()){
        logfile.open(log_file_name, std::ofstream::app);
        std::clog.rdbuf(logfile.rdbuf());
    }

    std::clog << "------------------------------\n";
    std::clog << "x0 = " << x << std::endl;
    for(int i = 0; i < n and std::abs(f(x)) > 1e-6; ++i){
        x = x-f(x)/df(x);
        std::clog << std::fixed << std::setprecision(6) << "x" << i+1 << " = " << x << std::endl;
    }

    if(!log_file_name.empty()){
        std::clog.rdbuf(old_rdbuf);
    }
    return x;
}

long double equation::solve_secant(const Function& f, long double x_prev,
                                   long double x_next, const std::string& log_file_name, int n){
    auto old_rdbuf = std::clog.rdbuf();
    std::ofstream logfile;
    if(!log_file_name.empty()){
        logfile.open(log_file_name, std::ofstream::app);
        std::clog.rdbuf(logfile.rdbuf());
    }

    std::clog << "------------------------------\n";
    std::clog << "x0 = " << x_prev << std::endl;
    std::clog << "x1 = " << x_next << std::endl;
    long double f_prev, f_next = f(x_prev);
    for(int i = 0; i < n and std::abs(f(x_next)) > 1e-6; ++i){
        f_prev = f_next;
        f_next = f(x_next);
        long double x_prev_prev = x_prev;
        x_prev = x_next;
        x_next = x_prev - ((x_prev-x_prev_prev)*f_next)/(f_next-f_prev);
        std::clog << std::fixed << std::setprecision(6) << "x" << i+1 << " = " << x_next << std::endl;
    }

    if(!log_file_name.empty()){
        std::clog.rdbuf(old_rdbuf);
    }
    return x_next;
}