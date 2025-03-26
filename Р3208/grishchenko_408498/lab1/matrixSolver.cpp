#include "matrixSolver.h"

MatrixSolver::MatrixSolver(int n, double** matrix, double* b_vector)
    : n(n), matrix(matrix), b_vector(b_vector)
{
    augmented_matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        augmented_matrix[i] = new double[n+1];
        for (int j = 0; j < n; j++) augmented_matrix[i][j] = matrix[i][j];
        augmented_matrix[i][n] = b_vector[i];
    }
}

MatrixSolver::MatrixSolver(int n, double** augmented_matrix)
    : n(n), augmented_matrix(augmented_matrix)
{
    matrix = new double*[n];
    b_vector = new double[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++) matrix[i][j] = augmented_matrix[i][j];
        b_vector[i] = augmented_matrix[i][n];
    }
}

MatrixSolver::~MatrixSolver()
{
    if (matrix != nullptr) {
        for (int i = 0; i < n; ++i) {
            delete[] matrix[i];
        }
        delete[] matrix;
        matrix = nullptr;
    }

    if (b_vector != nullptr) {
        delete[] b_vector;
        b_vector = nullptr;
    }

    if (augmented_matrix != nullptr) {
        for (int i = 0; i < n; ++i) {
            delete[] augmented_matrix[i];
        }
        delete[] augmented_matrix;
        augmented_matrix = nullptr;
    }

    if (triangular != nullptr) {
        for (int i = 0; i < n; ++i) {
            delete[] triangular[i];
        }
        delete[] triangular;
        triangular = nullptr;
    }

    if (solution != nullptr) {
        delete[] solution;
        solution = nullptr;
    }
}

void MatrixSolver::parse(int n, double** matrix, double* b_vector) {
    this->n = n;
    this->matrix = matrix;
    this->b_vector = b_vector;

    augmented_matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        augmented_matrix[i] = new double[n+1];
        for (int j = 0; j < n; j++) augmented_matrix[i][j] = matrix[i][j];
        augmented_matrix[i][n] = b_vector[i];
    }

    triangular = nullptr;
    solution = nullptr;
    eigen_solution.resize(0);
}

void MatrixSolver::parse(int n, double** augmented_matrix) {
    this->n = n;
    this->augmented_matrix = augmented_matrix;

    matrix = new double*[n];
    b_vector = new double[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++) matrix[i][j] = augmented_matrix[i][j];
        b_vector[i] = augmented_matrix[i][n];
    }

    triangular = nullptr;
    solution = nullptr;
    eigen_solution.resize(0);
}

double** MatrixSolver::copy_augmented(int n, double** matrix) {
    double** copy = new double*[n];
    for (int i = 0; i < n; i++) {
        copy[i] = new double[n+1];
        for (int j = 0; j < n + 1; j++) copy[i][j] = matrix[i][j];
    }
    return copy;
}

double** MatrixSolver::calc_triangular() {
    determinant_sign = 1.0;

    triangular = copy_augmented(n, augmented_matrix);

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
    return triangular;
}

double MatrixSolver::calc_determinant() {
    if (triangular == nullptr) calc_triangular();
    
    determinant = determinant_sign;
    for (int i = 0; i < n; ++i) {
        determinant *= triangular[i][i];
    }
    return determinant;
}

double* MatrixSolver::solve() {
    if (triangular == nullptr) calc_triangular();

    solution = new double[n];
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += triangular[i][j] * solution[j];
        }
        solution[i] = (triangular[i][n] - sum) / triangular[i][i];
    }
    return solution;
}

std::pair<double, double*> MatrixSolver::eigen_solve() {
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

double* MatrixSolver::calc_error() {
    if (solution == nullptr) solve();
    if (eigen_solution.size() == 0) eigen_solve();

    double* error = new double[n];
    for (int i = 0; i < n; i++) {
        error[i] = solution[i] - eigen_solution(i); 
    }
    return error;
}

void MatrixSolver::print_system() {
    std::cout << std::fixed << std::setprecision(1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(4) << (j == 0 ? matrix[i][j] : fabs(matrix[i][j])) << "x" << j + 1;
            if (j != n - 1) {
                std::cout << (matrix[i][j+1] > 0 ? " +" : " -");
            }
        }
        std::cout << (b_vector[i] > 0 ? " =  " : " = ") << b_vector[i] << std::endl;
    }
    std::cout << std::defaultfloat;
}