#include "integralSolver.h"
#include "utils.hpp"

Equation::Equation(std::string repr, std::function<double(double)> expr)
    :  repr(repr), expr(std::move(expr)) { }

double Equation::eval(double x) const {
    return expr(x);
}

std::ostream& operator<<(std::ostream& os, const Equation& eq) {
    os << eq.repr;
    return os;
}


IntegralSolver::IntegralSolver() { }

double IntegralSolver::solve(const Equation& eq, IntegralMethod method) {
    eps = prompt_double_default("Введите точность (1e-2): ", 1e-2);
    a = prompt_double("Введите нижнюю границу интегрирования: ");
    b = prompt_double("Введите верхнюю границу интегрирования: ");
    h = (b - a) / n;
    
    double current_h = h;

    int p = 1;
    switch (method) {
        case IntegralMethod::leftRectangle:
        case IntegralMethod::rightRectangle:
        case IntegralMethod::midRectangle:
            p = 1; break;
        case IntegralMethod::trapezoid:
            p = 2; break;
        case IntegralMethod::simpson:
            p = 4; break;
        default:
            throw std::invalid_argument("Unknown method");
    }

    double I_h, I_h2, R;

    auto compute_integral = [&](double step) {
        h = step;
        switch(method) {
            case IntegralMethod::leftRectangle:
                return solve_left_rectangle(eq);
            case IntegralMethod::midRectangle:
                return solve_mid_rectangle(eq);
            case IntegralMethod::rightRectangle:
                return solve_right_rectangle(eq);
            case IntegralMethod::trapezoid:
                return solve_trapezoid(eq);
            case IntegralMethod::simpson:
                return solve_simpson(eq);
            default:
                throw std::invalid_argument("Unknown method");
        }
    };

    I_h = compute_integral(current_h);
    I_h2 = compute_integral(current_h / 2);

    R = std::abs(I_h2 - I_h) / (std::pow(2.0, p) - 1);

    int i = 1;
    while (R > eps) {
        current_h /= 2;
        I_h = I_h2;
        I_h2 = compute_integral(current_h / 2);
        R = std::abs(I_h2 - I_h) / (std::pow(2.0, p) - 1);

        i++;
        if (i > max_iter) exit_error("Unable to achieve " + std::to_string(eps) + " accuracy within " + std::to_string(max_iter) + " iterations");
    }

    needed_n = n * pow(2, i);
    return I_h2;
}

double IntegralSolver::solve_left_rectangle(const Equation& eq) {
    std::vector<double> vec = make_divide(eq.expr, a, b, h);
    return h * std::accumulate(vec.begin(), vec.end() - 1, 0.0);
}

double IntegralSolver::solve_mid_rectangle(const Equation& eq) {
    std::vector<double> vec = make_divide(eq.expr, (a + h) / 2, b, h);
    return h * std::accumulate(vec.begin(), vec.end(), 0.0);
}

double IntegralSolver::solve_right_rectangle(const Equation& eq) {
    std::vector<double> vec = make_divide(eq.expr, a, b, h);
    return h * std::accumulate(vec.begin() + 1, vec.end(), 0.0);
}

double IntegralSolver::solve_trapezoid(const Equation& eq) {
    std::vector<double> vec = make_divide(eq.expr, a, b, h);
    return h * ((vec.front() + vec.back()) / 2 + std::accumulate(vec.begin() + 1, vec.end() - 1, 0.0));
}

double IntegralSolver::solve_simpson(const Equation& eq) {
    std::vector<double> vec = make_divide(eq.expr, a, b, h);
    double even_sum = 0, odd_sum = 0;
    for (int i = 1; i < vec.size() - 1; i++) {
        if (i % 2 == 0) even_sum += vec[i];
        else odd_sum += vec[i];
    }
    return (h / 3) * (vec.front() + vec.back() + 4 * odd_sum + 2 * even_sum);
}

std::vector<double> IntegralSolver::make_divide(const std::function<double(double)> f, double a, double b, double h) {
    std::vector<double> vec;
    for (double i = a; i <= b; i += h) {
        vec.push_back(f(i));
    }
    return vec;
}