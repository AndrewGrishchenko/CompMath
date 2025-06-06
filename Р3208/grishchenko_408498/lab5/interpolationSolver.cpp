#include "interpolationSolver.h"

InterpolationSolver::InterpolationSolver() { }

InterpolationResult InterpolationSolver::solve(const std::vector<double>& x, const std::vector<double>& y, double x0, InterpolationMethod method) {
    switch (method) {
        case InterpolationMethod::Lagrange:
            return solveLagrange(x, y, x0);
        case InterpolationMethod::NewtonDivided:
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

    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
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
    double result = 0;

    for (int i = 0; i < table.size(); i++) {
        double part = table[i][0];
        for (int j = 0; j < i; j++) {
            part *= (x0 - x[j]);
        }
        result += part;
    }

    std::vector<std::string> header = {"x", "y"};
    for (int i = 1; i < x.size(); i++) {
        header.push_back("f[...,...]");
    }

    table.insert(table.begin(), x);

    return {result, header, transposeMatrix<double>(table), "Таблица разделенных разностей"};
}

InterpolationResult InterpolationSolver::solveGauss(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    size_t n = x.size();
    auto table = finiteDifferencesTable(y);

    int center = 0;
    if (n % 2 == 0) {
        int lmid = n / 2 - 1;
        center = std::abs(x[lmid] - x0) < std::abs(x[lmid+1] - x0) ? lmid + 1: lmid;
    } else {
        center = n / 2;
    }
    
    double h = x[1] - x[0];
    double t = (x0 - x[center]) / h;
    double result = table[0][center];

    for (int i = 1; i < table.size(); i++) {
        int idx;
        double coeff;

        if (x0 >= x[center]) {
            if (i % 2 == 1)
                idx = center - (i / 2);
            else
                idx = center - (i / 2);

            if (idx < 0 || idx >= table[i].size()) break;

            coeff = 1.0;
            for (int j = 0; j < i; j++) {
                if (j % 2 != 0) coeff *= (t - (j / 2 + 1));
                else coeff *= (t + j / 2);
            }
        } else {
            if (i % 2 != 0)
                idx = center - (i / 2) - 1;
            else
                idx = center - (i / 2);

            if (idx < 0 || idx >= table[i].size()) break;

            coeff = 1.0;
            for (int j = 0; j < i; j++) {
                if (j % 2 != 0) coeff *= (t + j / 2 + 1);
                else coeff *= (t - j / 2);
            }
        }
        result += coeff * table[i][idx] / fact(i);
    }


    std::vector<std::string> header = {"x", "y"};
    for (int i = 1; i < x.size(); i++) {
        header.push_back("Δ^" + std::to_string(i));
    }

    table.insert(table.begin(), x);

    return {result, header, transposeMatrix<double>(table), "Таблица конечных разностей"};
}

int InterpolationSolver::fact(int n) {
    return (n <= 1) ? 1 : n * fact(n - 1);
}

std::vector<std::vector<double>> InterpolationSolver::dividedDifferencesTable(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    std::vector<std::vector<double>> table(n);

    table[0] = y;

    for (size_t i = 1; i < n; i++) {
        for (size_t j = 0; j < n - i; j++) {
            table[i].push_back((table[i-1][j+1] - table[i-1][j]) / (x[i+j] - x[j]));
        }
    }

    return table;
}

int InterpolationSolver::findClosestIndex(const std::vector<double>& x, double x0) {
    if (x.size() % 2 == 0) {
        int lmid = x.size() / 2 - 1;
        return std::abs(x[lmid] - x0) < std::abs(x[lmid+1] - x0) ? lmid + 1: lmid;
    } else {
        return x.size() / 2;
    }

}

std::vector<std::vector<double>> InterpolationSolver::finiteDifferencesTable(const std::vector<double>& y) {
    size_t n = y.size();
    std::vector<std::vector<double>> table(n);

    table[0] = y;

    for (size_t i = 1; i < n; i++) {
        for (size_t j = 0; j < n - i; j++) {
            table[i].push_back(table[i-1][j+1] - table[i-1][j]);
        }
    }

    return table;
}

