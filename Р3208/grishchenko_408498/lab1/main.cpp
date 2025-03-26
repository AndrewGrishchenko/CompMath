#include <iostream>

#include "utils.h"
#include "matrixSolver.h"

int main() {
    std::cout << "Решение СЛАУ методом Гаусса" << std::endl;

    select_input_mode();
    
    int n = prompt_int("Введите n: ");
    if (n <= 0) exit_error("n должно быть >= 1");

    double** augmented_matrix = prompt_augmented_matrix("Введите построчно расширенную матрицу", n);

    MatrixSolver* solver = new MatrixSolver();
    std::cout << std::endl << "СЛАУ:" << std::endl;
    solver->print_system(n, augmented_matrix);

    double** triangular = solver->calc_triangular(n, augmented_matrix).second;
    std::cout << std::endl << "Треугольная матрица:" << std::endl;
    print_augmented_matrix(n, triangular);

    double determinant = solver->calc_determinant(n, augmented_matrix);
    std::cout << std::endl << "Определитель: " << determinant << std::endl;

    double* solution = solver->solve(n, augmented_matrix);
    std::cout << "Решение системы методом Гаусса" << std::endl;
    print_vector("x", n, solution);

    auto [eigen_determinant, eigen_solution] = solver->eigen_solve(n, augmented_matrix);
    std::cout << std::endl << "Определитель через Eigen: " << eigen_determinant << std::endl;
    
    std::cout << std::endl << "Решение системы через Eigen" << std::endl;
    print_vector("x", n, eigen_solution);

    double* error = solver->calc_error(n, solution, eigen_solution);
    std::cout << std::endl << "Невязка:" << std::endl;
    print_vector("r", n, error);

    delete solver;
    return 0;
}