#include <iostream>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <cmath>

namespace Homework1 {
    using std::cin, std::cout, std::endl;
    using std::fixed, std::setprecision;
    using std::ifstream;
    using std::sqrt, std::abs, std::sin, std::cos;
    using std::numbers::pi;

    const double eps = 1e-8;

    void solve(double p, double x, double y) {
        double l = pi / 2, r = pi;
        do {
            double theta = (l + r) / 2;
            double mirror_x = p * cos(2 * theta);
            double mirror_y = p * sin(2 * theta);
            double t_x = cos(theta);
            double t_y = sin(theta);
            double d_x = 
                (y - mirror_y) * (mirror_x - t_x) / (mirror_y - t_y) 
                + mirror_x;
            if (abs(d_x - x) <= eps) {
                r = theta;
                break;
            }
            if (d_x > x) l = theta;
            else r = theta;
        } while (abs(r - l) >= eps);
        double theta = r;
        double t_x = cos(theta);
        double t_y = sin(theta);
        double k_pb = t_y / (t_x - p);
        double db = sqrt((y - t_y) * (y - t_y) + (x - t_x) * (x - t_x));
        double r_x = t_x + db / sqrt(1 + k_pb * k_pb);
        double r_y = t_y + db * k_pb / sqrt(1 + k_pb * k_pb);
        cout << fixed << setprecision(6) << "T = (" << t_x << ", " << t_y << ") "
            << "R = (" << r_x << ", " << r_y << ")" << endl;
    }

    void test() {
        ifstream fin("Homework1.txt");
        int cnt;
        fin >> cnt;
        for (int i = 0; i < cnt; ++i) {
            double p, x, y;
            fin >> p >> x >> y;
            solve(p, x, y);
        }
        fin.close();
    }
}
