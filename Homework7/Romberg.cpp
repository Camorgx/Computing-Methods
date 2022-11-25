#include "Romberg.h"

#include <vector>

int failed_to_meet_eps_cnt;
int total_called_cnt;

double Romberg(double a, double b, double eps, int M, double (*f)(double x)) {
    ++total_called_cnt;
    std::vector<std::vector<double>> r;
    double h = b - a;
    r.emplace_back();
    r.emplace_back(2);
    r[1][1] = (f(a) + f(b)) * h / 2;
    int k = 2;
    for (; k <= M; ++k) {
        r.emplace_back(k + 1);
        double tmp = 0;
        int up = 1 << (k - 2);
        for (int i = 1; i <= up; ++i)
            tmp += f(a + (2 * i - 1) * h / (1 << (k - 1)));
        r[k][1] = (r[k - 1][1] + h / (1 << (k - 2)) * tmp) / 2;
        for (int j = 2; j <= k; ++j) {
            r[k][j] = r[k][j - 1]
                + (r[k][j - 1] - r[k - 1][j - 1]) / ((1 << (2 * (j - 1))) - 1);
        }
        if (abs(r[k][k] - r[k - 1][k - 1]) < eps) break;
    }
    int res_index;
    if (k > M) {
        res_index = k - 1;
        ++failed_to_meet_eps_cnt;
    }
    else res_index = k;
    return r[res_index][res_index];
}
