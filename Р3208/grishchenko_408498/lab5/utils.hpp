#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <utility>
#include <optional>
#include "math.h"

#include <termios.h>
#include <unistd.h>

#include <sciplot/sciplot.hpp>

inline struct termios orig_termios;

//==interactive==//

inline void disableRawMode() {
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &orig_termios);
}

inline void enableRawMode() {
    tcgetattr(STDIN_FILENO, &orig_termios);
    atexit(disableRawMode);

    struct termios raw = orig_termios;
    raw.c_lflag &= ~(ECHO | ICANON);
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);
}

inline char readKey() {
    char c;
    read(STDIN_FILENO, &c, 1);
    if (c == '\033') {
        char seq[2];
        read(STDIN_FILENO, &seq[0], 1);
        read(STDIN_FILENO, &seq[1], 1);
        if (seq[0] == '[') {
            if (seq[1] == 'A') return 'u';
            if (seq[1] == 'B') return 'd';
        }
        return 0;
    } else {
        return c;
    }
}

inline void clearOptions(size_t N) {
    for (int i = 0; i < N; i++) {
        std::cout << "\033[F";
        std::cout << "\033[2K";
    }
}

template<typename T, size_t N>
inline void printOptions(const T (&opts)[N], int selected) {
    for (auto i = 0; i < N; i++) {
        if (i == selected) {
            std::cout << "\033[1;32m" << "-> " << opts[i] << "\033[0m\n";
        } else {
            std::cout << "   " << opts[i] << "\n";
        }
    }
}

template<typename T, size_t N>
inline int promptIndexOptions(std::string msg, const T (&opts)[N]) {
    enableRawMode();
    std::cout << msg << std::endl;

    int selected = 0;
    while (true) {
        printOptions(opts, selected);

        char key = readKey();
        if (key == 'u') {
            selected = (selected - 1 + N) % N;
        } else if (key == 'd') {
            selected = (selected + 1) % N;
        } else if (key == '\n') {
            break;
        }

        clearOptions(N);
    }

    disableRawMode();
    return selected;
}

template<typename T, size_t N>
inline T promptOptions(std::string msg, const T (&opts)[N]) {
    return opts[promptIndexOptions(msg, opts)];
}

//==non-interactive==//

inline int READ_MODE = 1;
inline std::ifstream file_in;

inline void exit_error(std::string msg) {
    std::cerr << "Ошибка: " << msg << std::endl;
    exit(1);
}

template<typename T>
constexpr bool is_vector_v = false;

template<typename U, typename Alloc>
constexpr bool is_vector_v<std::vector<U, Alloc>> = true;

template<typename T>
constexpr bool is_pair_v = false;

template<typename A, typename B>
constexpr bool is_pair_v<std::pair<A, B>> = true;

template<typename T>
T read_cli(std::optional<T> default_value = std::nullopt) {
    std::istream* in = READ_MODE == 1 ? &std::cin : &file_in;
    
    std::string line;
    std::getline(*in, line);
    if (line.empty()) {
        if (default_value.has_value()) return default_value.value();
        else throw std::runtime_error("empty input");
    }
    if (READ_MODE == 2) std::cout << line << std::endl;
    std::replace(line.begin(), line.end(), ',', '.');
    std::stringstream ss(line);

    if constexpr (std::is_same_v<T, std::string>) {
        return line;
    } else if constexpr (std::is_same_v<T, int>) {
        T val;
        ss >> val; //TODO: check double, int btw
        return val;
    } else if constexpr (std::is_same_v<T, double>) {
        T val;
        ss >> val;
        return val;
    } else if constexpr (is_vector_v<T>) {
        T vec;
        typename T::value_type num;

        while (ss >> num) {
            vec.push_back(num);
        }
        return vec;
    } else if constexpr (is_pair_v<T>) {
        typename T::first_type t1;
        typename T::second_type t2;

        ss >> t1;
        ss >> t2;

        return std::make_pair(t1, t2);
    } else {
        static_assert(!sizeof(T), "Unsupported type for read()");
    }
}

template<typename T>
T prompt(const std::string& msg, std::optional<T> default_value = std::nullopt) {
    std::cout << msg << " ";
    return read_cli<T>(default_value);
}

void printTable(const std::vector<std::string>& header, const std::vector<std::vector<double>>& data, int colWidth = 15) {
    for (const auto& h : header) {
        std::cout << std::setw(colWidth) << h;
    }
    std::cout << '\n';

    for (const auto& row : data) {
        for (size_t i = 0; i < header.size(); ++i) {
            if (i < row.size()) {
                std::cout << std::setw(colWidth) << std::fixed << std::setprecision(6) << row[i];
            } else {
                std::cout << std::setw(colWidth) << "";
            }
        }
        std::cout << '\n';
    }
}

//==plot==//

void make_plot(std::vector<std::function<double(double)>> f, std::vector<std::pair<double, double>> intervals, std::vector<std::vector<std::pair<double, double>>> points, std::vector<std::string> legend) {
    using namespace sciplot;

    int POINT_COUNT = 1000;

    Plot2D plot;
    plot.size(1920, 1080);
    plot.fontSize(3);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    int legendIndexOffset = 0;
    double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
    if (f.size() != 0) {
        if (f.size() != intervals.size()) throw std::runtime_error("f size doesnt meet intervals size");

        std::vector<std::vector<double>> x(f.size());
        std::vector<std::vector<double>> y(f.size());

        std::vector<double> steps;
        for (int i = 0; i < f.size(); i++) {
            if (intervals[i].first >= intervals[i].second) throw std::runtime_error("invalid interval");
            
            steps.push_back((intervals[i].second - intervals[i].first) / (POINT_COUNT - 1));
            double val = intervals[i].first;
            for (int j = 0; j < POINT_COUNT; j++) {
                x[i].push_back(val);
                y[i].push_back(f[i](val));

                xmin = std::min(xmin, val);
                xmax = std::max(xmax, val);
                ymin = std::min(ymin, y[i].front());
                ymax = std::max(ymax, y[i].front());

                val += steps[i];
            }
        }

        for (int i = 0; i < x.size(); i++) {
            plot.drawCurve(x[i], y[i]).label(legend[i]).lineWidth(2);
        }
        legendIndexOffset = x.size();
    } else {
        for (auto row : points) {
            for (auto& val : row) {
                xmin = std::min(xmin, val.first);
                xmax = std::max(xmax, val.first);
                ymin = std::min(ymin, val.second);
                ymax = std::max(ymax, val.second);
            }
        }
    }

    plot.xrange(xmin - 1, xmax + 1);
    plot.yrange(ymin - 1, ymax + 1);

    for (int i = 0; i < points.size(); i++) {
        std::vector<double> pointX, pointY;
        for (int j = 0; j < points[i].size(); j++) {
            pointX.push_back(points[i][j].first);
            pointY.push_back(points[i][j].second);
        }
        plot.drawPoints(pointX, pointY).pointType(0).lineWidth(3).label(legend[legendIndexOffset + i]);
    }

    Figure fig = {{plot}};
    Canvas canvas = {{fig}};
    canvas.fontSize(5);
    canvas.show();
}

//==utility==//

template<typename T1, typename T2>
std::vector<std::pair<T1, T2>> vec2vecpair(std::vector<T1> vec1, std::vector<T2> vec2) {
    if (vec1.size() != vec2.size()) throw std::runtime_error("vec1 size doesnt meet vec2 size");
    
    std::vector<std::pair<T1, T2>> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = {vec1[i], vec2[i]};
    }

    return result;
}

//==file==//

inline bool file_exists(std::string filename) {
    if (FILE* file = fopen(filename.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

inline void file_mode() {
    std::string filename = prompt<std::string>("Название файла (input.txt):", "input.txt");
    
    if (!file_exists(filename)) exit_error("no such file");   
    
    std::cout << "Используется файл ввода \"" << filename << "\"" << std::endl;
    file_in.open(filename.c_str());
    READ_MODE = 2;
}

inline void cli_mode() {
    READ_MODE = 1;
}