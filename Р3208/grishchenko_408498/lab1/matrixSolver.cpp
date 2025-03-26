#include "matrixSolver.h"

MatrixSolver::MatrixSolver() { }

MatrixSolver::~MatrixSolver() { }

double** MatrixSolver::copy_augmented(int n, double** matrix) {
    double** copy = new double*[n];
    for (int i = 0; i < n; i++) {
        copy[i] = new double[n+1];
        for (int j = 0; j < n + 1; j++) copy[i][j] = matrix[i][j];
    }
    return copy;
}

std::pair<double, double**> MatrixSolver::calc_triangular(int n, double** augmented_matrix) {
    double determinant_sign = 1.0;

    double** triangular = new double*[n];
    for (int i = 0; i < n; i++) {
        triangular[i] = new double[n+1];
        for (int j = 0; j <= n; j++) triangular[i][j] = augmented_matrix[i][j];
    }

    for (int j = 0; j < n; ++j) {
        int pivot_row = j;
        for (int i = j + 1; i < n; ++i) {
            if (fabs(triangular[i][j]) > fabs(triangular[pivot_row][j])) {
                pivot_row = i;
            }
        }

        if (triangular[pivot_row][j] == 0) {
            std::cout << "Матрица вырожденная, решение не может быть найдено однозначно." << std::endl;
            exit(1);
        }

        if (pivot_row != j) {
            std::swap(triangular[j], triangular[pivot_row]);
            determinant_sign *= -1.0;
        }

        for (int i = j + 1; i < n; ++i) {
            double factor = triangular[i][j] / triangular[j][j];
            for (int k = j; k <= n; ++k) {
                triangular[i][k] -= factor * triangular[j][k];
            }
        }
    }
    return std::pair<double, double**>(determinant_sign, triangular);
}

double MatrixSolver::calc_determinant(int n, double** augmented_matrix) {
    auto [determinant_sign, triangular] = calc_triangular(n, augmented_matrix);
    
    double determinant = determinant_sign;
    for (int i = 0; i < n; ++i) {
        determinant *= triangular[i][i];
    }
    return determinant;
}

double* MatrixSolver::solve(int n, double** augmented_matrix) {
    double** triangular = calc_triangular(n, augmented_matrix).second;

    double* solution = new double[n];
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += triangular[i][j] * solution[j];
        }
        solution[i] = (triangular[i][n] - sum) / triangular[i][i];
    }
    return solution;
}

std::pair<double, double*> MatrixSolver::eigen_solve(int n, double** augmented_matrix) {
    eigen_a.resize(n, n);
    eigen_b.resize(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            eigen_a(i, j) = augmented_matrix[i][j];
        }
        eigen_b(i) = augmented_matrix[i][n];
    }

    eigen_determinant = eigen_a.determinant();
    eigen_solution = eigen_a.colPivHouseholderQr().solve(eigen_b);
    return std::pair<double, double*>(eigen_determinant, eigen_solution.data());
}

double* MatrixSolver::calc_residual(int n, double** augmented_matrix, double* solution) {
    double* residual = new double[n];
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) {
            s += augmented_matrix[i][j] * solution[j];
        }
        residual[i] = s - augmented_matrix[i][n];
    }

    return residual;
}

void MatrixSolver::print_system(int n, double** augmented_matrix) {
    std::cout << std::fixed << std::setprecision(1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(8) << (j == 0 ? augmented_matrix[i][j] : fabs(augmented_matrix[i][j])) << "x" << j + 1;
            if (j != n - 1) {
                std::cout << std::setw(4) << (augmented_matrix[i][j+1] > 0 ? "+" : "-");
            }
        }
        std::cout << (augmented_matrix[i][n] > 0 ? " =  " : " = ") << augmented_matrix[i][n] << std::endl;
    }
}