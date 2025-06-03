#include "interpolationSolver.h"

InterpolationSolver::InterpolationSolver() { }

InterpolationResult InterpolationSolver::solve(const std::vector<double>& x, const std::vector<double>& y, double x0, InterpolationMethod method) {
    switch (method) {
        case InterpolationMethod::Lagrange:
            return solveLagrange(x, y, x0);
        case InterpolationMethod::NetwonDivided:
            return solveNewtonDivided(x, y, x0);
        case InterpolationMethod::Gauss:
            return solveGauss(x, y, x0);
        default:
            throw std::runtime_error("Unknown method");
    }
}

InterpolationResult InterpolationSolver::solveLagrange(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    int n = x.size();
    double result = 0.0;

    for (int i = 0; i < n; ++i) {
        double term = y[i];
        for (int j = 0; j < n; ++j) {
            if (j != i)
                term *= (x0 - x[j]) / (x[i] - x[j]);
        }
        result += term;
    }

    std::vector<std::string> header = {"x", "y"};
    std::vector<std::vector<double>> table(x.size(), std::vector<double>(2));
    for (int i = 0; i < x.size(); i++) {
        table[i][0] = x[i];
        table[i][1] = y[i];
    }

    return {result, header, table, "Таблица значений"};
}

InterpolationResult InterpolationSolver::solveNewtonDivided(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    auto table = dividedDifferencesTable(x, y);
    int n = x.size();
    double result = table[0][0];

    double product = 1.0;
    for (int order = 1; order < n; ++order) {
        product *= (x0 - x[order-1]);
        result += table[0][order] * product;
    }

    std::vector<std::string> header = {"x", "y"};
    for (int i = 1; i < x.size(); i++) {
        header.push_back("f[...,...]");
    }

    return {result, header, table, "Таблица разделенных разностей"};
}

InterpolationResult InterpolationSolver::solveGauss(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    int n = x.size();

    double h = x[1] - x[0];
    for (int i = 2; i < n; ++i) {
        if (abs((x[i] - x[i-1]) - h) > 1e-12) {
            std::cerr << "Error: non-uniform step for Gauss interpolation\n";
            return {0, {}, {}};
        }
    }

    auto table = finiteDifferencesTable(x, y);

    double t = (x0 - x[0]) / h;

    double result = y[0];
    double term = 1.0;

    double factorial = 1.0;

    for (int order = 1; order < n; order++) {
        term *= (t - (order - 1));
        factorial *= order;
        result += (term / factorial) * table[0][order];
    }

    std::vector<std::string> header = {"x", "y"};
    for (int i = 1; i < x.size(); i++) {
        header.push_back("Δ^" + std::to_string(i));
    }

    return {result, header, table, "Таблица конечных разностей"};
}

std::vector<std::vector<double>> InterpolationSolver::dividedDifferencesTable(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    std::vector<std::vector<double>> table(n);



    for (int i = 0; i < n; ++i) {
        table[i].push_back(x[i]);
        table[i].push_back(y[i]);
    }

    for (int order = 1; order < n; ++order) {
        for (int i = 0; i < n - order; ++i) {
            double numerator = table[i+1][order-1] - table[i][order-1];
            double denominator = x[i+order] - x[i];
            table[i].push_back(numerator / denominator);
        }
    }
    return table;
}

std::vector<std::vector<double>> InterpolationSolver::finiteDifferencesTable(const std::vector<double>& x, const std::vector<double>& y) {
    int n = y.size();
    std::vector<std::vector<double>> table(n);

    for (int i = 0; i < n; ++i) {
        table[i].push_back(x[i]);
        table[i].push_back(y[i]);
    }

    for (int order = 1; order < n; ++order) {
        for (int i = 0; i < n - order; ++i) {
            double diff = table[i+1][order-1] - table[i][order-1];
            table[i].push_back(diff);
        }
    }
    return table;
}

