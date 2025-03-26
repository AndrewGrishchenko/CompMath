#include <iostream>

#include "utils.h"
#include "matrixSolver.h"

int main() {
    std::cout << "Решение СЛАУ методом Гаусса" << std::endl;

    select_input_mode();
    
    int n = prompt_int("Введите n: ");
    if (n <= 0) exit_error("n должно быть >= 1");

    double** augmented_matrix = prompt_augmented_matrix("Введите построчно расширенную матрицу", n);

    MatrixSolver* solver = new MatrixSolver(n, augmented_matrix);
    std::cout << std::endl << "СЛАУ:" << std::endl;
    solver->print_system();

    double** triangular = solver->calc_triangular();
    std::cout << std::endl << "Треугольная матрица:" << std::endl;
    print_augmented_matrix(n, triangular);

    double determinant = solver->calc_determinant();
    std::cout << std::endl << "Определитель: " << determinant << std::endl;

    double* solution = solver->solve();
    std::cout << "Решение системы методом Гаусса" << std::endl;
    print_vector("x", n, solution);

    std::pair<double, double*> eigen_solution = solver->eigen_solve();
    std::cout << std::endl << "Определитель через Eigen: " << eigen_solution.first << std::endl;
    
    std::cout << std::endl << "Решение системы через Eigen" << std::endl;
    print_vector("x", n, eigen_solution.second);

    double* error = solver->calc_error();
    std::cout << std::endl << "Невязка:" << std::endl;
    print_vector("r", n, error);

    delete solver;
    return 0;
}