#ifndef _EQUATIONS_SOLVER_H
#define _EQUATIONS_SOLVER_H

#include <iostream>
#include <sstream>
#include <functional>
#include <vector>
#include <variant>
#include "math.h"
#include <matplot/matplot.h>

class Equation {
    public:
        Equation(std::string repr, std::function<double(double)> expr);

        double eval(double x) const;

        friend std::ostream& operator<<(std::ostream& os, const Equation& eq);

        const std::string repr;
        const std::function<double(double)> expr;
};

class EquationSystem {
    public:
        EquationSystem(std::string repr, std::vector<std::function<double(double, double)>> exprs);

        std::vector<double> eval(double x, double y) const;

        friend std::ostream& operator<<(std::ostream& os, const EquationSystem& eq_system);
    
        const std::string repr;
        std::vector<std::function<double(double, double)>> exprs;
};

enum class EquationMethod {
    Bisection,
    Newton,
    SimpleIteration
};

enum class EquationSystemMethod {
    Newton
};

template<typename T>
struct MethodChoice {
    T method;
    std::string name;

    friend std::ostream& operator<<(std::ostream& os, const MethodChoice& mc) {
        return os << mc.name;
    }
};

class EquationSolver {
    public:
        EquationSolver();

        double solve(const Equation& eq, EquationMethod method);
        std::pair<double, double> solve_system(const EquationSystem& eq_system, EquationSystemMethod method);

        std::string get_log() const { return log_iss.str(); }
        void make_plot(size_t samples = 500);
        void make_system_plot(size_t samples = 500);
    private:
        double solve_bisection(const Equation& eq);
        double solve_newton(const Equation& eq);
        double solve_simple_iteration(const Equation& eq);

        std::pair<double, double> solve_system_newton(const EquationSystem& eq_system);

        std::function<double(double)> numerical_derivative(std::function<double(double)> f, double h = 1e-6);
        std::function<double(double)> build_phi(std::function<double(double)> f, double lambda);

        int iter_count;
        std::ostringstream log_iss;

        double graph_a, graph_b, graph_x, graph_y;
        std::function<double(double)> graph_f;
        std::function<double(double, double)> graph_f1, graph_f2;
};

#endif