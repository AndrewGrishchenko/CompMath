#ifndef _ODE_SOLVER_H
#define _ODE_SOLVER_H

#include <iostream>
#include <sstream>
#include <functional>
#include <vector>

#include "utils.hpp"

template<typename T>
class Equation {
    public:
        Equation(std::function<T(T)> expr)
            : repr(""), expr(expr) { }
        Equation(std::string repr, std::function<T(T)> expr)
            : repr(repr), expr(expr) { }

        T eval(T val) const { return expr(val); }

        const std::string repr;
        const std::function<T(T)> expr;

        friend std::ostream& operator<<(std::ostream& os, const Equation<T>& eq) {
            os << eq.repr;
            return os;
        }
};

template<typename T>
class Equation2D {
    public:
        Equation2D(std::function<T(T,T)> expr)
            : repr(""), expr(expr) { }
        Equation2D(std::string repr, std::function<T(T, T)> expr)
            : repr(repr), expr(expr) { }
        
        T eval(T val1, T val2) const { return expr(val1, val2); }

        const std::string repr;
        const std::function<T(T, T)> expr;

        friend std::ostream& operator<<(std::ostream& os, const Equation2D<T>& eq) {
            os << eq.repr;
            return os;
        }
};

template<typename T>
class ODEEquation {
    public:
        ODEEquation(std::string repr, std::function<T(T, T)> expr, std::function<T(T)> exact)
            : repr(repr), expr(Equation2D<T>(expr)), exact(Equation<T>(exact)) { }
        ODEEquation(std::string repr, Equation2D<T> expr, Equation<T> exact)
            : repr(repr), expr(expr), exact(exact) { }

        const std::string repr;
        const Equation2D<T> expr;
        const Equation<T> exact;

        friend std::ostream& operator<<(std::ostream& os, const ODEEquation<T>& eq) {
            os << eq.repr;
            return os;
        }
};

enum class ODEMethod {
    ImprovedEuler,
    Runge4,
    Millne
};

template<typename T>
struct MethodChoice {
    T method;
    std::string name;

    friend std::ostream& operator<<(std::ostream& os, const MethodChoice& mc) {
        return os << mc.name;
    }
};

struct ODEResult {
    std::vector<std::string> tableHeader;
    std::vector<std::vector<double>> table;
    double result;
    std::vector<std::pair<double, double>> points;
    std::vector<std::pair<double, double>> exactPoints;
};

class ODESolver {
    public:
        ODESolver() { }

        ODEResult solve(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h, double eps, ODEMethod method);
    private:
        ODEResult solveImprovedEuler(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h);
        ODEResult solveRunge4(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h);
        ODEResult solveMillne(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double eps, double h);
};

#endif