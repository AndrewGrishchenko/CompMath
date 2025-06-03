#ifndef _INTERPOLATION_SOLVER_H
#define _INTERPOLATION_SOLVER_H

#include <iostream>
#include <sstream>
#include <vector>
#include <functional>

template<typename T>
class Equation {
    public:
        Equation(std::string repr, std::function<T(T)> expr)
            : repr(repr), expr(expr) { }

        T eval(T val) { return expr(val); }

        const std::string repr;
        const std::function<T(T)> expr;

        friend std::ostream& operator<<(std::ostream& os, const Equation<T>& eq) {
            os << eq.repr;
            return os;
        }
};

enum class InterpolationMethod {
    Lagrange,
    NetwonDivided,
    Gauss
};

template<typename T>
struct MethodChoice {
    T method;
    std::string name;

    friend std::ostream& operator<<(std::ostream& os, const MethodChoice& mc) {
        return os << mc.name;
    }
};

struct InterpolationResult {
    double value;
    std::vector<std::string> diffTableHeader;
    std::vector<std::vector<double>> diffTable;
    std::string name;
};

class InterpolationSolver {
    public:
        InterpolationSolver();

        InterpolationResult solve(const std::vector<double>& x, const std::vector<double>& y, double x0, InterpolationMethod method);

    private:
        InterpolationResult solveLagrange(const std::vector<double>& x, const std::vector<double>& y, double x0);
        InterpolationResult solveNewtonDivided(const std::vector<double>& x, const std::vector<double>& y, double x0);
        InterpolationResult solveGauss(const std::vector<double>& x, const std::vector<double>& y, double x0);

        std::vector<std::vector<double>> dividedDifferencesTable(const std::vector<double>& x, const std::vector<double>& y);
        std::vector<std::vector<double>> finiteDifferencesTable(const std::vector<double>& x, const std::vector<double>& y);
};

#endif