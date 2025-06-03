#include <iostream>
#include <iomanip>
#include <cmath>

#include "utils.hpp"
#include "interpolationSolver.h"

int main() {
    std::vector<double> x, y;

    const std::string modes[] = { "Вручную", "Из файла", "С помощью функции" };
    const Equation<double> eqs[] = {
        {"sin(x)", [](double x) { return sin(x); }},
        {"1/(1+x^2)", [](double x) { return 1 / (1 + pow(x, 2)); }},
        {"x^3 - 2x + 1", [](double x) { return pow(x, 3) - 2 * x + 1; }}
    };
    const MethodChoice<InterpolationMethod> ims[] = {
        { InterpolationMethod::Lagrange, "Многочлен Лагранжа" },
        { InterpolationMethod::NetwonDivided, "Многочлен Ньютона с разделенными разностями" },
        { InterpolationMethod::Gauss, "Многочлен Гаусса" }
    };

    std::vector<std::function<double(double)>> fs;
    std::string f_name;
    std::vector<std::pair<double, double>> intervals;
    std::vector<std::vector<std::pair<double, double>>> points;
    std::vector<std::pair<double, double>> func_points;
    std::pair<double, double> interpolation_point;

    InterpolationMethod method = promptOptions("Выберите метод", ims).method;
    int mode = promptIndexOptions("Выберите режим", modes);

    if (mode == 0) {
        x = prompt<std::vector<double>>("Введите x:");
        y = prompt<std::vector<double>>("Введите y:");

        if (x.size() != y.size()) exit_error("Размеры x и y не совпадают");
        points.push_back(vec2vecpair<double, double>(x, y));
    } else if (mode == 1) {
        file_mode();
        x = prompt<std::vector<double>>("Введите x:");
        y = prompt<std::vector<double>>("Введите y:");
        points.push_back(vec2vecpair<double, double>(x, y));
    } else {
        Equation<double> eq = promptOptions("Выберите уравнение", eqs);
        std::pair<double, double> interval = prompt<std::pair<double, double>>("Введите интервалы:");
        if (interval.first >= interval.second) exit_error("invalid interval");
        int num = prompt<int>("Ввелите количество точек:");
        
        double step = (interval.second - interval.first) / (num - 1);
        double val = interval.first;
        for (int i = 0; i < num; i++) {
            x.push_back(val);
            y.push_back(eq.eval(val));
            val += step;
        }

        fs.push_back(eq.expr);
        f_name = eq.repr;
        intervals.push_back({interval.first, interval.second});
    }

    double x0 = prompt<double>("Введите x0:");
    
    InterpolationSolver solver;
    InterpolationResult result = solver.solve(x, y, x0, method);

    std::cout << "f(x) ≈ " << result.value << std::endl;
    std::cout << result.name << std::endl;
    printTable(result.diffTableHeader, result.diffTable);

    points.push_back({{x0, result.value}});

    make_plot(fs, intervals, points, {fs.empty() ? "points" : f_name, "interpolation point"});
}