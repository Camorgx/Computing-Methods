#include "Romberg.h"

#include <cmath>
#include <format>
#include <iostream>
#include <vector>

double ax(double t) {
    return std::sin(t) / (std::sqrt(t) + 1);
}

double ay(double t) {
    return std::log(t + 1) / (t + 1);
}

const double eps = 1e-6;
int M;

double vx(double t) {
    return Romberg(0, t, eps, M, ax);
}

double vy(double t) {
    return Romberg(0, t, eps, M, ay);
}

std::string to_string(const std::vector<double>& vec) {
    std::string res = "[ ";
    for (const double item : vec)
        res += std::format("{0:.6f} ", item);
    res += "]";
    return res;
}

int main() {
    const int M_cases[] = { 4, 8, 12, 16, 20 };
    const double t_cases[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    for (const int m : M_cases) {
        std::cout << std::format("M = {0}\n", M = m);
        failed_to_meet_eps_cnt = 0;
        total_called_cnt = 0;
        std::vector<double> point_x, point_y;
        for (const double t : t_cases) {
            double x = Romberg(0, t, eps, M, vx);
            double y = Romberg(0, t, eps, M, vy);
            point_x.emplace_back(x);
            point_y.emplace_back(y);
        }
        std::cout << std::format("x = {0}\n", to_string(point_x));
        std::cout << std::format("y = {0}\n", to_string(point_y));
        std::cout << std::format("percentage not meet eps = {0:.2f}%\n", 
            static_cast<double>(failed_to_meet_eps_cnt) / total_called_cnt);
        std::cout << std::endl;
    }
    return 0;
}
