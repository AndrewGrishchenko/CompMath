#include "ODESolver.h"

ODEResult ODESolver::solve(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h, double eps, ODEMethod method) {
    switch (method) {
        case ODEMethod::ImprovedEuler:
            return solveImprovedEuler(eq, y0, interval, h);
        case ODEMethod::Runge4:
            return solveRunge4(eq, y0, interval, h);
        case ODEMethod::Millne:
            return solveMillne(eq, y0, interval, h, eps);
        default:
            throw std::runtime_error("unknown method");
    }
}

ODEResult ODESolver::solveImprovedEuler(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h) {
    //i xi yi f(xi, yi) exact
    std::vector<std::vector<double>> table;
    std::vector<std::pair<double, double>> points;
    std::vector<std::pair<double, double>> exactPoints;

    table.push_back({0, interval.first, y0, eq.expr.eval(interval.first, y0), eq.exact.eval(interval.first)});
    points.push_back({interval.first, y0});
    exactPoints.push_back({interval.first, eq.exact.eval(interval.first)});

    int i = 1;
    for (double xi = interval.first + h; xi <= interval.second; xi += h) {
        double yi = table[i-1][2] + (h / 2) * (table[i-1][3] + eq.expr.eval(xi, table[i-1][2] + h * table[i-1][3]));
        table.push_back({(double) i, xi, yi, eq.expr.eval(xi, yi), eq.exact.eval(xi)});
        points.push_back({xi, yi});
        exactPoints.push_back({xi, eq.exact.eval(xi)});
        i++;
    }

    std::vector<std::string> header = {"i", "xi", "yi", "f(xi, yi)", "exact value"};
    
    return {header, table, table.front()[2], points, exactPoints};
}

ODEResult ODESolver::solveRunge4(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h) {
    //i xi yi f(xi, yi) exact
    std::vector<std::vector<double>> table;
    std::vector<std::pair<double, double>> points;
    std::vector<std::pair<double, double>> exactPoints;

    table.push_back({0, interval.first, y0, eq.expr.eval(interval.first, y0), eq.exact.eval(interval.first)});
    points.push_back({interval.first, y0});
    exactPoints.push_back({interval.first, eq.exact.eval(interval.first)});

    int i = 1;
    for (double xi = interval.first + h; xi <= interval.second; xi += h) {
        double k1 = h * table[i-1][3];
        double k2 = h * eq.expr.eval(table[i-1][1] + h / 2, table[i-1][2] + k1 / 2);
        double k3 = h * eq.expr.eval(table[i-1][1] + h / 2, table[i-1][2] + k2 / 2);
        double k4 = h * eq.expr.eval(table[i-1][1] + h, table[i-1][2] + k3);
        double yi = table[i-1][2] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        table.push_back({(double) i, xi, yi, eq.expr.eval(xi, yi), eq.exact.eval(xi)});
        points.push_back({xi, yi});
        exactPoints.push_back({xi, eq.exact.eval(xi)});
        i++;
    }
    
    std::vector<std::string> header = {"i", "xi", "yi", "f(xi, yi)", "exact"};

    return {header, table, table.front()[2], points, exactPoints};
}

ODEResult ODESolver::solveMillne(const ODEEquation<double>& eq, double y0, std::pair<double, double> interval, double h, double eps) {
    //i xi yi f(xi, yi) exact y_predict y_corr
    std::vector<std::vector<double>> table;
    std::vector<std::pair<double, double>> points;
    std::vector<std::pair<double, double>> exactPoints;

    points.push_back({interval.first, y0});
    exactPoints.push_back({interval.first, eq.exact.eval(interval.first)});

    ODEResult runge4 = solveRunge4(eq, y0, interval, h);
    for (int i = 0; i < 4; i++) {
        table.push_back({runge4.table[i][0], runge4.table[i][1], runge4.table[i][2], runge4.table[i][3], runge4.table[i][4]});
    }

    int i = 4;
    for (double xi = interval.first + 4 * h; xi <= interval.second; xi += h) {
        double y_predict = table[i-4][2] + (4.0 * h / 3.0) * (2 * table[i-3][3] - table[i-2][3] + 2 * table[i-1][3]);
        double y_corr = table[i-2][2] + (h / 3) * (table[i-2][3] + 4 * table[i-1][3] + eq.expr.eval(xi, y_predict));
        while (std::abs(y_corr - y_predict) > eps) {
            y_predict = y_corr;
            y_corr = table[i-2][2] + (h / 3) * (table[i-2][3] + 4 * table[i-1][3] + eq.expr.eval(xi, y_predict));
        }
        table.push_back({(double) i, xi, y_corr, eq.expr.eval(xi, y_corr), eq.exact.eval(xi), y_predict, y_corr});
        points.push_back({xi, y_corr});
        exactPoints.push_back({xi, eq.exact.eval(xi)});
        i++;
    }

    std::vector<std::string> header = {"i", "xi", "yi", "f(xi, yi)", "exact", "y_predict", "y_corr"};
    
    return {header, table, table.front()[2], points, exactPoints};
}