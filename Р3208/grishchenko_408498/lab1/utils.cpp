#include "utils.h"

int READ_MODE;

void exit_error(std::string msg) {
    std::cerr << "Ошибка: " << msg << std::endl;
    exit(1);
}

bool file_exists(std::string filename) {
    if (FILE* file = fopen(filename.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void select_input_mode() {
    std::cout << "Выберите режим ввода матрицы:" << std::endl;
    std::cout << "1 - ввод из консоли" << std::endl;
    std::cout << "2 - ввод из файла" << std::endl;
    READ_MODE = read_int();
    if (READ_MODE != 1 && READ_MODE != 2) {
        exit_error("Неверный режим ввода матрицы");
    } else if (READ_MODE == 1) std::cout << "Выбран ввод из консоли" << std::endl;
    else {
        std::cout << "Выбран ввод из файла" << std::endl;
        
        std::string filename = "input.txt";
        if (!file_exists(filename)) {
            while (true) {
                std::cout << "Введите название файла: ";
                std::cin >> filename;
                if (file_exists(filename)) break;
                std::cout << "Файл не найден" << std::endl;
            }
        }

        std::cout << "Используется файл ввода \"" << filename << "\"" << std::endl;
        freopen(filename.c_str(), "r", stdin);
    }
}

std::string read_string() {
    std::string str;
    std::cin >> str;
    if (READ_MODE == 2) std::cout << str;
    return str;
}

double read_double() {
    std::string str;
    std::getline(std::cin, str);

    std::replace(str.begin(), str.end(), ',', '.');
    std::stringstream ss(str);
    
    double num;
    ss >> num;

    if (!ss.eof()) exit_error("Слишком много аргументов");

    if (READ_MODE == 2) std::cout << num;
    return num;
}

int read_int() {
    double num = read_double();
    if (num != static_cast<int>(num)) exit_error("Expected int");
    return num;
}

std::string prompt_string(std::string msg) {
    std::cout << msg;
    return read_string();
}

int prompt_int(std::string msg) {
    std::cout << msg;
    
    int n = read_int();
    if (READ_MODE == 2) std::cout << std::endl;
    return n;
}

double** prompt_augmented_matrix(std::string msg, int n) {
    std::cout << msg << std::endl;
    
    bool flag = true;
    std::string str;

    double** matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n+1];

        std::getline(std::cin, str);
        if (str.empty() && flag) {
            delete matrix[i], matrix;
            return generate_augmented_matrix(n);
        }

        std::replace(str.begin(), str.end(), ',', '.');
        std::stringstream ss(str);
        for (int j = 0; j < n + 1; j++) {
            if (!(ss >> matrix[i][j])) exit_error("Не хватает элементов строки");
        }
        if (!ss.eof()) exit_error("Слишком много элементов строки");
        if (READ_MODE == 2) std::cout << str << std::endl;
        flag = false;
    }
    return matrix;
}

double** generate_augmented_matrix(int n) {
    std::cout << "Генерация случайной матрицы:" << std::endl;
    std::srand(std::time(0));

    double** matrix = new double*[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n+1];
        for (int j = 0; j <= n; j++) 
            matrix[i][j] = (std::rand() % 101) + (std::rand() % 1000) / 1000.0;
    }

    print_augmented_matrix(n, matrix);

    return matrix;
}

void print_augmented_matrix(int n, double** matrix) {
    std::cout << std::fixed << std::setprecision(2);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(8) << matrix[i][j];
        }
        std::cout << std::setw(4) << " | " << std::setw(4) << matrix[i][n] << std::endl;
    }
}

void print_vector(std::string variable, int n, double* vector) {
    for (int i = 0; i < n; i++) {
        std::cout << variable << std::left << std::setw(3) << std::to_string(i + 1);
        std::cout << (vector[i] > 0 ? "=  " : "= -") << fabs(vector[i]) << std::endl;
    }
    std::cout << std::right;
}