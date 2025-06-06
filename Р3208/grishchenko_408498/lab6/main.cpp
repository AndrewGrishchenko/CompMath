#include <iostream>
#include <iomanip>
#include <cmath>

#include "utils.hpp"
#include "ODESolver.h"

int main() {
    const ODEEquation<double> eqs[] = {
        {"y' = y", [](double x, double y) { return y; }, [](double x) { return exp(x); }}, //y(0) = 1
        {"y' = x", [](double x, double y) { return x; }, [](double x) { return pow(x, 2) / 2; }}, //y(0) = 0
        {"y' = x*y", [](double x, double y) { return x * y; }, [](double x) { return exp(pow(x, 2) / 2); }} // y(0) = 1
    };
    const MethodChoice<ODEMethod> methods[] = {
        { ODEMethod::ImprovedEuler, "Усовершенствованный метод Эйлера" },
        { ODEMethod::Runge4, "Метод Рунге-Кутта 4-го порядка" },
        { ODEMethod::Millne, "Метод Милна" }
    };

    ODESolver solver;
    ODEEquation<double> eq = promptOptions("Выберите уравнение", eqs);
    ODEMethod method = promptOptions("Выберите метод", methods).method;

    std::pair<double, double> interval = prompt<std::pair<double, double>>("Введите интервал:");
    double h = prompt<double>("Введите h:");
    double y0 = prompt<double>("Введите y(" + std::to_string(interval.first) + "):");
    double eps = prompt<double>("Введите eps (1e-3):", 1e-3);

    ODEResult result = solver.solve(eq, y0, interval, h, eps, method);

    if (method == ODEMethod::ImprovedEuler) {
        double p = 1;
        double h_new = h;
        while(true) {
            double y_old = result.result;
            h_new /= 2;
            result = solver.solve(eq, y0, interval, h_new, eps, method);
            if (std::abs(y_old - result.result) / (pow(2, p) - 1) < eps) break;
            y_old = result.result;
        }
    } else if (method == ODEMethod::Runge4) {
        double p = 4;
        double h_new = h;
        while(true) {
            double y_old = result.result;
            h_new /= 2;
            result = solver.solve(eq, y0, interval, h_new, eps, method);
            if (std::abs(y_old - result.result) / (pow(2, p) - 1) < eps) break;
            y_old = result.result;
        }
    }

    printTable(result.tableHeader, result.table);

    make_plot({}, {interval}, {result.points, result.exactPoints}, {"Приближенное решение", "Точное решение"});
}