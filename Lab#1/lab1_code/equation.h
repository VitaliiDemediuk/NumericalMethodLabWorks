#ifndef LAB1_CODE_EQUATION_H
#define LAB1_CODE_EQUATION_H
#include <functional>

namespace equation{
    using Function = std::function<long double (long double)>;

    long double solve_relaxation(const Function& f, long double const_t,
                                 long double x0, const std::string& log_file_name = "", int n = 100);

    long double solve_newton(const Function& f, const Function& df,
                                 long double x0, const std::string& log_file_name = "", int n = 100);
    long double solve_secant(const Function& f, long double x0, long double x1,
                             const std::string& log_file_name = "", int n = 100);
}

#endif //LAB1_CODE_EQUATION_H
