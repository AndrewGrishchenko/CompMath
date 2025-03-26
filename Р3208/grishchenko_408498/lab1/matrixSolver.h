#ifndef _MATRIX_SOLVER_H
#define _MATRIX_SOLVER_H

#include <iostream>
#include <iomanip>
#include "math.h"
#include <vector>

#include <Eigen/Dense>

class MatrixSolver {
    public:
        MatrixSolver();
        ~MatrixSolver();

        void print_system(int n, double** augmented_matrix);

        std::pair<double, double**> calc_triangular(int n, double** augmented_matrix);
        double calc_determinant(int n, double** augmented_matrix);
        double* solve(int n, double** augmented_matrix);

        std::pair<double, double*> eigen_solve(int n, double** augmented_matrix);

        double* calc_residual(int n, double** augmented_matrix, double* solution);

    private:
        double determinant_sign;
        
        double** copy_augmented(int n, double** source);

        Eigen::MatrixXd eigen_a;
        Eigen::VectorXd eigen_b;
        double eigen_determinant;
        Eigen::VectorXd eigen_solution;
};

#endif