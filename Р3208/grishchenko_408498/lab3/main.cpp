#include <iostream>
#include <iomanip>
#include <cmath>

#include "utils.hpp"
#include "integralSolver.h"
#include <numeric>

int main() {
    const Equation ieqs[] = {
        Equation("sin(x)", [](double x) { return sin(x); }),
        Equation("cos(x)", [](double x) { return cos(x); }),
        Equation("x^2", [](double x) { return pow(x, 2); }),
        Equation("e^x", [](double x) { return exp(x); }),
        Equation("-3x^3 - 5x^2 + 4x - 2", [](double x) { return -3 * pow(x, 3) - 5 * pow(x, 2) + 4 * x - 2; }),
        Equation("5x^2 + 3x", [](double x) { return 5 * pow(x, 2) + 3 * x; })
    };

    const MethodChoice<IntegralMethod> ieq_methods[] = {
        { IntegralMethod::leftRectangle, "Метод прямоугольников (левых)" },
        { IntegralMethod::midRectangle, "Метод прямоугольников (средних)" },
        { IntegralMethod::rightRectangle, "Метод прямоугольников (правых)" },
        { IntegralMethod::trapezoid, "Метод трапеций" },
        { IntegralMethod::simpson, "Метод Симпсона" }
    };

    IntegralSolver solver;


    Equation ieq = prompt_multi("Выберите подынтегральную функцию:", ieqs);
    auto method = prompt_multi("Выберите метод:", ieq_methods);

    double ans = solver.solve(ieq, method.method);
    std::cout << "Значение интеграла: " << std::setprecision(17) << ans << std::endl;
    std::cout << "Необходимый n: " << solver.needed_n << std::endl;
}