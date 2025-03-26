#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "math.h"

void exit_error(std::string msg);
bool file_exists(std::string filename);
double** generate_augmented_matrix(int n);

void select_input_mode();

std::string read_string();
double read_double();
int read_int();

std::string prompt_string(std::string msg);
int prompt_int(std::string msg, bool newline = false);
double** prompt_augmented_matrix(std::string msg, int n);

void print_augmented_matrix(int n, double** matrix);
void print_vector(std::string variable, int n, double* vector);

#endif