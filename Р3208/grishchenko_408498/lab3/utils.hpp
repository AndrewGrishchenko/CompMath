#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "math.h"

inline int READ_MODE = 1;
inline std::ifstream in;

inline void exit_error(std::string msg) {
    std::cerr << "Ошибка: " << msg << std::endl;
    exit(1);
}

inline bool file_exists(std::string filename) {
    if (FILE* file = fopen(filename.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

inline std::string read_string() {
    std::string str;
    if (READ_MODE == 2) {
        in >> str;
        std::cout << str;
    } else std::cin >> str;
    return str;
}

inline double read_double() {
    std::string str;
    std::getline(READ_MODE == 2 ? in : std::cin, str);

    std::replace(str.begin(), str.end(), ',', '.');
    std::stringstream ss(str);
    
    double num;
    
    if (!(ss >> num)) exit_error("Ошибка чтения");

    if (!ss.eof()) exit_error("Слишком много аргументов");

    if (READ_MODE == 2) std::cout << num << std::endl;
    return num;
}

inline double read_double_default(double value) {
    std::string str;
    std::getline(READ_MODE == 2 ? in : std::cin, str);
    
    std::replace(str.begin(), str.end(), ',', '.');
    std::stringstream ss(str);

    double num;
    
    if (!(ss >> num)) {
        return value;
    }

    if (!ss.eof()) {
        return value;
    }

    if (READ_MODE == 2) std::cout << num << std::endl;
    return num;
}

inline double prompt_double(std::string msg) {
    std::cout << msg;
    return read_double();
}

inline double prompt_double_default(std::string msg, double value) {
    std::cout << msg;
    return read_double_default(value);
}

inline int read_int() {
    double num = read_double();
    if (num != static_cast<int>(num)) exit_error("Expected int");
    return num;
}

inline void select_input_mode() {
    std::cout << "Выберите режим ввода:" << std::endl;
    std::cout << "1 - ввод из консоли" << std::endl;
    std::cout << "2 - ввод из файла" << std::endl;
    READ_MODE = read_int();
    if (READ_MODE != 1 && READ_MODE != 2) {
        exit_error("Неверный режим ввода");
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
        // freopen(filename.c_str(), "r", stdin);
        in.open(filename.c_str());
    }
}

inline std::string prompt_string(std::string msg) {
    std::cout << msg;
    return read_string();
}

template<typename... Args>
inline int prompt_multi(std::string msg, const Args&... opts) {
    std::cout << msg << std::endl;
    
    int index = 1;
    auto print = [&](const auto& opt) {
        std::cout << index++ << ". " << opt << std::endl;
    };

    (print(opts), ...);

    int prompted = read_int();
    if (READ_MODE == 2) std::cout << prompted << std::endl;
    if (prompted >= index) exit_error("Неверный аргумент");
    return prompted;
}

template<typename T, size_t N>
inline T prompt_multi(std::string msg, const T (&arr)[N]) {
    std::cout << msg << std::endl;

    auto i = 0;
    while (i < N) {
        std::cout << (i + 1) << ". " << arr[i] << std::endl;
        i++;
    }

    int prompted = read_int();
    if (READ_MODE == 2) std::cout << prompted << std::endl;
    if (prompted >= i+1) exit_error("Неверный аргумент");
    return arr[prompted - 1];
}

inline int prompt_int(std::string msg) {
    std::cout << msg;
    
    int n = read_int();
    if (READ_MODE == 2) std::cout << std::endl;
    return n;
}