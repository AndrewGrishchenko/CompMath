#ifndef _MATRIX_SOLVER_H
#define _MATRIX_SOLVER_H

#include <iostream>
#include <iomanip>
#include "math.h"
#include <vector>

#include <Eigen/Dense>

class MatrixSolver {
    public:
        MatrixSolver(int n, double** matrix, double* b_vector);
        MatrixSolver(int n, double** augmented_matrix);
        ~MatrixSolver();

        void parse(int n, double** matrix, double* b_vector);
        void parse(int n, double** augmented_matrix);
        void print_system();

        double** calc_triangular();
        double calc_determinant();
        double* solve();

        std::pair<double, double*> eigen_solve();

        double* calc_error();

    private:
        int n;
        double** matrix = nullptr;
        double* b_vector = nullptr;

        double** augmented_matrix = nullptr;
        double** triangular = nullptr;

        double determinant_sign;
        double determinant;
        double* solution = nullptr;

        double** copy_augmented(int n, double** source);

        Eigen::MatrixXd eigen_a;
        Eigen::VectorXd eigen_b;
        double eigen_determinant;
        Eigen::VectorXd eigen_solution;
};

#endif