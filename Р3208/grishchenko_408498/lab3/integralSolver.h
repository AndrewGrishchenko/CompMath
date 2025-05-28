#ifndef _INTEGRAL_SOLVER_H
#define _INTEGRAL_SOLVER_H

#include <iostream>
#include <sstream>
#include <functional>
#include <vector>
#include <numeric>

class Equation {
    public:
        Equation(std::string repr, std::function<double(double)> expr);

        double eval(double x) const;

        friend std::ostream& operator<<(std::ostream& os, const Equation& eq);

        const std::string repr;
        const std::function<double(double)> expr;
};

enum class IntegralMethod {
    leftRectangle,
    rightRectangle,
    midRectangle,
    trapezoid,
    simpson
};

template<typename T>
struct MethodChoice {
    T method;
    std::string name;

    friend std::ostream& operator<<(std::ostream& os, const MethodChoice& mc) {
        return os << mc.name;
    }
};

class IntegralSolver {
    public:
        IntegralSolver();

        double solve(const Equation& eq, IntegralMethod method);

        int needed_n;
    private:
        double solve_left_rectangle(const Equation& eq);
        double solve_mid_rectangle(const Equation& eq);
        double solve_right_rectangle(const Equation& eq);
        double solve_trapezoid(const Equation& eq);
        double solve_simpson(const Equation& eq);

        std::vector<double> make_divide(const std::function<double(double)> f, double a, double b, double h);

        const double n = 4;
        const int max_iter = 1000;

        double eps, a, b, h;
};

#endif