#include <cmath>
#include <vector>
#include <iostream>

using std::cout, std::endl;

const double eps = 1e-5;
std::vector<std::vector<double>> X, Y;
std::vector<double> egvalue;

using matrix = std::vector<std::vector<double>>;
using vector = std::vector<double>;

// m: number of rows
// n: number of colomns
matrix init_matrix(size_t m, size_t n) {
    matrix res;
    for (int i = 0; i < m; ++i)
        res.emplace_back(n);
    return res;
}

void Doolittle_decompose(const matrix& a, matrix& l, matrix& u) {
    size_t n = a.size();
    l = init_matrix(n, n);
    u = init_matrix(n, n);
    for (int i = 0; i < n; ++i)
        l[i][i] = 1;
    for (int k = 0; k < n; ++k) {
        for (int j = k; j < n; ++j) {
            u[k][j] = a[k][j];
            for (int r = 0; r < k; ++r)
                u[k][j] -= l[k][r] * u[r][j];
        }
        for (int i = k + 1; i < n; ++i) {
            l[i][k] = a[i][k];
            for (int r = 0; r < k; ++r)
                l[i][k] -= l[i][r] * u[r][k];
            l[i][k] /= u[k][k];
        }
    }
}

vector Doolittle_solve(
    const vector& b, const matrix& l, const matrix& u) {
    size_t n = b.size();
    vector y(n), x(n);
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j)
            y[i] -= l[i][j] * y[j];
    }
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= u[i][j] * x[j];
        x[i] /= u[i][i];
    }
    return x;
}

void display_matrix(const matrix& mat) {
    size_t m = mat.size();
    size_t n = mat[0].size();
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            std::cout << mat[i][j] << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

double get_norm(const matrix& mat) {
    double ans = 0;
    size_t m = mat.size();
    size_t n = mat[0].size();
    for (int i = 0; i < m; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j)
            sum += std::abs(mat[i][j]);
        ans = std::max(ans, sum);
    }
    return ans;
}

double get_norm(const vector& vec) {
    double norm = 0;
    for (double x : vec)
        norm = std::max(norm, abs(x));
    return norm;
}

int anti_exponentiation(const matrix& a, double& eigenvalue, vector& eigenvector) {
    X.clear();
    Y.clear();
    egvalue.clear();
    size_t n = a.size();
    vector x0(n, 1);
    X.emplace_back(x0);
    matrix l, u;
    Doolittle_decompose(a, l, u);
    int cnt = 0;
    while (true) {
        ++cnt;
        double eigevalue0 = get_norm(x0);
        vector y(n);
        for (int i = 0; i < n; ++i)
            y[i] = x0[i] / eigevalue0;
        Y.emplace_back(y);
        vector x1 = Doolittle_solve(y, l, u);
        X.emplace_back(x1);
        double eigevalue1 = get_norm(x1);
        egvalue.emplace_back(1 / eigevalue1);
        if (std::abs(1 / eigevalue1 - 1 / eigevalue0) < eps) {
            eigenvector = std::move(y);
            eigenvalue = 1 / eigevalue1;
            return cnt;
        }
        x0 = x1;
    }
}

void display_vector(const vector& vec) {
    cout << '(';
    for (int i = 0; i < vec.size() - 1; ++i)
        cout << vec[i] << ", ";
    cout << vec[vec.size() - 1] << ')' << endl;
}

int main() {
    matrix A1 = {
        { 1.0 / 9, 1.0 / 8, 1.0 / 7, 1.0 / 6, 1.0 / 5 },
        { 1.0 / 8, 1.0 / 7, 1.0 / 6, 1.0 / 5, 1.0 / 4 },
        { 1.0 / 7, 1.0 / 6, 1.0 / 5, 1.0 / 4, 1.0 / 3 },
        { 1.0 / 6, 1.0 / 5, 1.0 / 4, 1.0 / 3, 1.0 / 2 },
        { 1.0 / 5, 1.0 / 4, 1.0 / 3, 1.0 / 2, 1.0 / 1 }
    };
    matrix A2 = {
        { 4, -1, 1, 3 },
        { 16, -2, -2, 5 },
        { 16, -3, -1, 7 },
        { 6, -4, 2, 9 }
    };
    matrix test[] = { A1, A2 };
    
    for (const matrix& test : test) {
        double eigenvalue;
        vector eigenvector;
        int cnt = anti_exponentiation(test, eigenvalue, eigenvector);
        cout << "cnt = " << cnt << endl << endl;
        for (int i = 0; i < X.size(); ++i) {
            cout << "X" << i << " = ";
            display_vector(X[i]);
        }
        cout << endl;
        for (int i = 0; i < Y.size(); ++i) {
            cout << "Y" << i << " = ";
            display_vector(Y[i]);
        }
        cout << endl;
        for (int i = 0; i < egvalue.size(); ++i)
            cout << "lambda[" << i << "] = " << egvalue[i] << endl;
        cout << endl;
        cout << "eigenvalue = " << eigenvalue << endl;
        cout << "eigenvector = ";
        display_vector(eigenvector);
        cout << endl;
    }
}
