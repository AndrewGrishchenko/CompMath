#include <iostream>
#include <iomanip>
#include "math.h"

#include "utils.hpp"
#include "equationsSolver.h"

int main() {
    const Equation eqs[] = {
        Equation("x^3 - x + 4", [](double x) { return pow(x, 3) - x + 4; }),
        Equation("x^3 - x^2 - 25x + 2", [](double x) { return pow(x, 3) - pow(x, 2) - 25 * x + 2; }),
        Equation("arctg(x)", [](double x) { return atan(x); })
    };

    const EquationSystem eq_systems[] = {
        EquationSystem("{x^2 + y^2 = 4; y = 3x^2}", {
            [](double x, double y) { return pow(x, 2) + pow(y, 2) - 4; },
            [](double x, double y) { return y - 3 * pow(x, 2); }
        }),
        EquationSystem("{x^2 + y^2 = 4; y = sin(x)}", {
            [](double x, double y) { return pow(x, 2) + pow(y, 2) - 4; },
            [](double x, double y) { return y - sin(x); }
        })
    };

    const MethodChoice<EquationMethod> eq_methods[] = {
        { EquationMethod::Bisection, "Метод половинного деления" },
        { EquationMethod::Newton, "Метод Ньютона" },
        { EquationMethod::SimpleIteration, "Метод простой итерации" }
    };

    const MethodChoice<EquationSystemMethod> eq_system_methods[] = {
        { EquationSystemMethod::Newton, "Метод Ньютона" }
    };

    EquationSolver solver;

    int mode = prompt_multi("Выберите что нужно решить:", "Нелинейное уравнение", "Система нелинейных уравнений");
    
    if (mode == 1) {
        Equation eq = prompt_multi("Выберите уравнение:", eqs);
        auto method = prompt_multi("Выберите метод:", eq_methods);
        select_input_mode();
        std::cout << std::setprecision(17) << "Найденный корень: " 
                  << solver.solve(eq, method.method) << std::endl;
        std::cout << solver.get_log();
        solver.make_plot();
    } else {
        EquationSystem eq_system = prompt_multi("Выберите систему уравнений:", eq_systems);
        auto method = prompt_multi("Выберите метод:", eq_system_methods);
        select_input_mode();
        auto ans = solver.solve_system(eq_system, method.method);
        std::cout << std::setprecision(17) << "Найденный корень: " 
                  << "{" << ans.first << ", " << ans.second << "}" << std::endl;
        std::cout << solver.get_log();
        solver.make_system_plot();
    }
}