#include "equationsSolver.h"
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


EquationSystem::EquationSystem(std::string repr, std::vector<std::function<double(double, double)>> exprs)
    : repr(repr), exprs(exprs) { }

std::vector<double> EquationSystem::eval(double x, double y) const {
    std::vector<double> results;
    for (const auto& expr : exprs) {
        results.push_back(expr(x, y));
    }
    return results;
}

std::ostream& operator<<(std::ostream& os, const EquationSystem& eq_system) {
    os << eq_system.repr;
    return os;
}


EquationSolver::EquationSolver() { }

double EquationSolver::solve(const Equation& eq, EquationMethod method) {
    switch(method) {
        case EquationMethod::Bisection:
            return solve_bisection(eq);
        case EquationMethod::Newton:
            return solve_newton(eq);
        case EquationMethod::SimpleIteration:
            return solve_simple_iteration(eq);
        default:
            throw std::invalid_argument("Unknown method");
    }
}

std::pair<double, double> EquationSolver::solve_system(const EquationSystem& eq_system, EquationSystemMethod method) {
    switch(method) {
        case EquationSystemMethod::Newton:
            return solve_system_newton(eq_system);
        default:
            throw std::invalid_argument("Unknown method");
    }
}

double EquationSolver::solve_bisection(const Equation& eq) {
    std::cout << "Bisection sovler" << std::endl;
    std::cout << "Eq: " << eq << std::endl;

    double eps = prompt_double_default("Введите точность (1e-6): ", 1e-6);
    double a = prompt_double("Введите нижнюю границу корня: ");
    double b = prompt_double("Введите верхнюю границу корня: ");
    graph_a = a, graph_b = b;
    graph_f = eq.expr;

    if (eq.eval(a) * eq.eval(b) > 0)
        exit_error("На интервале нет корня или их несколько");

    log_iss.clear();
    log_iss << std::fixed << std::setprecision(6);
    log_iss << std::setw(10) << "Итерация" 
            << std::setw(15) << "a"
            << std::setw(15) << "b"
            << std::setw(15) << "mid"
            << std::setw(15) << "f(mid)"
            << '\n';
    log_iss << std::string(70, '-') << '\n';

    double mid;
    int iter = 0;
    while ((b - a) / 2 > eps) {
        mid = (a + b) / 2;

        double fmid = eq.eval(mid);
        log_iss << std::setw(10) << iter
                << std::setw(15) << a
                << std::setw(15) << b
                << std::setw(15) << mid
                << std::setw(15) << fmid
                << '\n';

        if (eq.eval(mid) == 0.0)
            return mid;
        else if (eq.eval(a) * eq.eval(mid) < 0)
            b = mid;
        else
            a = mid;
        
        iter++;
    }

    graph_x = (a + b) / 2;
    graph_y = eq.eval(graph_x);
    return graph_x;
}

double EquationSolver::solve_newton(const Equation& eq) {
    std::cout << "Newton solver" << std::endl;
    std::cout << "Eq: " << eq << std::endl;

    std::function<double(double)> f_prime = numerical_derivative(eq.expr);
    int max_iter = 1000;
    double eps = prompt_double_default("Введите точность (1e-6): ", 1e-6);
    double x0 = prompt_double("Введите x0: ");
    graph_a = x0 - 5, graph_b = x0 + 5;
    graph_f = eq.expr;

    log_iss.clear();
    log_iss << std::fixed << std::setprecision(6);
    log_iss << std::setw(10) << "Итерация"
            << std::setw(15) << "x"
            << std::setw(15) << "f(x)"
            << std::setw(15) << "f'(x)"
            << std::setw(15) << "x_new"
            << '\n';
    log_iss << std::string(70, '-') << '\n';

    double x = x0;
    int iter = 0;
    for (int i = 0; i < max_iter; ++i) {
        double fx = eq.eval(x);
        double fpx = f_prime(x);
        if (std::abs(fpx) < 1e-12)
            exit_error("Производная близка к нулю");

        double x_next = x - fx / fpx;

        log_iss << std::setw(10) << iter
                << std::setw(15) << x
                << std::setw(15) << fx
                << std::setw(15) << fpx
                << std::setw(15) << x_next
                << '\n';

        if (std::abs(x_next - x) < eps) {
            graph_x = x_next;
            graph_y = eq.eval(graph_x);
            return graph_x;
        }
        x = x_next;
        iter++;
    }
    exit_error("Метод Ньютона не сошелся за указанное число итераций");
    return 0;
}

double EquationSolver::solve_simple_iteration(const Equation& eq) {
    std::cout << "Simple iteration solver" << std::endl;
    std::cout << "Eq: " << eq << std::endl;

    std::function<double(double)> phi = build_phi(eq.expr, 0.1);

    double eps = prompt_double_default("Введите точность (1e-6): ", 1e-6);
    int max_iter = 1000;
    double x0 = prompt_double("Введите x0: ");
    graph_a = x0 - 5, graph_b = x0 + 5;
    graph_f = eq.expr;

    log_iss.clear();
    log_iss << std::fixed << std::setprecision(6);
    log_iss << std::setw(10) << "Итерация"
            << std::setw(15) << "x_old"
            << std::setw(15) << "x_new"
            << std::setw(15) << "|x_new - x_old|"
            << '\n';
    log_iss << std::string(60, '-') << '\n';

    double x_old = x0;
    double x_new = phi(x_old);
    int iter = 0;
    while (iter < max_iter) {
        double diff = fabs(x_new - x_old);

        log_iss << std::setw(10) << iter
                << std::setw(15) << x_old
                << std::setw(15) << x_new
                << std::setw(15) << diff
                << '\n';

        if (diff < eps)
            break;

        x_old = x_new;
        x_new = phi(x_old);
        iter++;
    }

    graph_x = x_new;
    graph_y = eq.eval(graph_x);
    return graph_x;
}

void EquationSolver::make_plot(size_t samples) {
    using namespace matplot;
    
    std::vector<double> x(samples);
    std::vector<double> y(samples);
    double step = (graph_b - graph_a) / (samples - 1);
    for (size_t i = 0; i < samples; ++i) {
        x[i] = graph_a + i * step;
        y[i] = graph_f(x[i]);
    }

    auto f = figure(false);
    f->backend()->run_command("unset warnings");
    f->ioff();
    f->size(1600, 900);

    auto l = plot(x, y)->line_width(2);
    hold(on);
    grid(on);

    auto ax = gca();
    auto x_limits = ax->xlim();
    auto y_limits = ax->ylim();
    xlim(x_limits);
    ylim(y_limits);

    std::vector<double> x_axis_x = {x_limits[0], x_limits[1]};
    std::vector<double> x_axis_y = {0, 0};
    auto x_line = plot(x_axis_x, x_axis_y, "k-");
    x_line->line_width(2.0);
    
    std::vector<double> y_axis_x = {0, 0};
    std::vector<double> y_axis_y = {y_limits[0], y_limits[1]};
    auto y_line = matplot::plot(y_axis_x, y_axis_y, "k-");
    y_line->line_width(2.0);
    
    std::vector<double> point_x = {graph_x};
    std::vector<double> point_y = {graph_y};
    
    auto point = plot(point_x, point_y, "-o");
    point->color("red")
        .marker_size(10)
        .marker_face_color("red")
        .line_width(0);

    show();
}

void EquationSolver::make_system_plot(size_t samples) {
    using namespace matplot;
    
    double step = (graph_b - graph_a) / (samples - 1); 
    std::vector<double> x_vals, y_vals;
    for (size_t i = 0; i < samples; i++) {
        x_vals.push_back(graph_a + i * step);
        y_vals.push_back(graph_a + i * step);
    }

    std::vector<std::vector<double>> Z1(y_vals.size(), std::vector<double>(x_vals.size()));
    std::vector<std::vector<double>> Z2(y_vals.size(), std::vector<double>(x_vals.size()));
    
    for (size_t i = 0; i < y_vals.size(); ++i) {
        for (size_t j = 0; j < x_vals.size(); ++j) {
            double x = x_vals[j];
            double y = y_vals[i];
            Z1[i][j] = graph_f1(x, y);
            Z2[i][j] = graph_f2(x, y);
        }
    }

    int nx = x_vals.size();
    int ny = y_vals.size();

    std::vector<std::vector<double>> X(ny, std::vector<double>(nx));
    std::vector<std::vector<double>> Y(ny, std::vector<double>(nx));

    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            X[i][j] = x_vals[j];
            Y[i][j] = y_vals[i];
        }
    }

    auto f = figure(false);
    f->backend()->run_command("unset warnings");
    f->ioff();
    f->size(1600, 900);

    hold(on);
    grid(on);

    auto c1 = contour(X, Y, Z1);
    c1->levels({0.0});
    c1->line_width(2);

    auto c2 = contour(X, Y, Z2);
    c2->levels({0.0});
    c2->line_width(2);

    gca()->color_box(false); 

    auto ax = gca();
    auto x_limits = ax->xlim();
    auto y_limits = ax->ylim();
    xlim(x_limits);
    ylim(y_limits);

    std::vector<double> x_axis_x = {x_limits[0], x_limits[1]};
    std::vector<double> x_axis_y = {0, 0};
    auto x_line = plot(x_axis_x, x_axis_y, "k-");
    x_line->line_width(2.0);
    
    std::vector<double> y_axis_x = {0, 0};
    std::vector<double> y_axis_y = {y_limits[0], y_limits[1]};
    auto y_line = matplot::plot(y_axis_x, y_axis_y, "k-");
    y_line->line_width(2.0);
    
    std::vector<double> point_x = {graph_x};
    std::vector<double> point_y = {graph_y};
    
    auto point = plot(point_x, point_y, "-o");
    point->color("red")
        .marker_size(10)
        .marker_face_color("red")
        .line_width(0);

    show();
}

std::function<double(double)> EquationSolver::numerical_derivative(std::function<double(double)> f, double h) {
    return [f, h](double x) {
        return (f(x + h) - f(x - h)) / (2 * h);
    };
}

std::function<double(double)> EquationSolver::build_phi(std::function<double(double)> f, double lambda) {
    return [f, lambda](double x) {
        return x - lambda * f(x);
    };
}

std::pair<double, double> EquationSolver::solve_system_newton(const EquationSystem& eq_system) {
    int max_iter = 100;
    double eps = prompt_double_default("Введите точность (1e-6): ", 1e-6);
    double x0 = prompt_double("Введите x0: ");
    double y0 = prompt_double("Введите y0: ");
    graph_a = x0 - 5, graph_b = x0 + 5;
    graph_f1 = eq_system.exprs[0], graph_f2 = eq_system.exprs[1];

    log_iss << std::fixed << std::setprecision(6);
    log_iss << std::setw(10) << "Итерация"
            << std::setw(15) << "x_old"
            << std::setw(15) << "y_old"
            << std::setw(15) << "x_new"
            << std::setw(15) << "y_new"
            << std::setw(20) << "||delta||"
            << '\n';
    log_iss << std::string(90, '-') << '\n';

    double x = x0, y = y0;
    for (int i = 0; i < max_iter; ++i) {
        auto F = eq_system.eval(x, y);

        double h = 1e-8;
        double dF1dx = (eq_system.eval(x + h, y)[0] - F[0]) / h;
        double dF1dy = (eq_system.eval(x, y + h)[0] - F[0]) / h;
        double dF2dx = (eq_system.eval(x + h, y)[1] - F[1]) / h;
        double dF2dy = (eq_system.eval(x, y + h)[1] - F[1]) / h;

        double J[2][2] = {
            {dF1dx, dF1dy},
            {dF2dx, dF2dy}
        };

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (std::abs(det) < 1e-12) {
            exit_error("Якобиан вырожден. Решение невозможно");
        }

        double invJ[2][2] = {
            { J[1][1] / det, -J[0][1] / det },
            { -J[1][0] / det, J[0][0] / det }
        };

        double dx = -(invJ[0][0] * F[0] + invJ[0][1] * F[1]);
        double dy = -(invJ[1][0] * F[0] + invJ[1][1] * F[1]);

        double x_old = x, y_old = y;
        x += dx;
        y += dy;

        double delta = std::sqrt(dx * dx + dy * dy);

        log_iss << std::setw(10) << i + 1
                << std::setw(15) << x_old
                << std::setw(15) << y_old
                << std::setw(15) << x
                << std::setw(15) << y
                << std::setw(20) << delta
                << '\n';

        if (delta < eps)
            break;
    }

    graph_x = x, graph_y = y;
    return {x, y};
}