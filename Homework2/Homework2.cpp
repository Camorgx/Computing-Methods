#include <iostream>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <vector>
#include <cmath>

const double Eps = 1e-6;
const double a = 0.5;
const int n = 100;
const double h = 1.0 / n;

using std::move;
using std::swap;
using std::abs, std::max, std::exp;
using std::cout, std::endl;
using std::ifstream;

using vector = std::vector<double>;
using matrix = std::vector<std::vector<double>>;

std::vector<std::vector<double>> A;
std::vector<double> b;

vector Gauss(const matrix& origin_A, const vector& origin_b) {
	matrix A(origin_A);
	vector b(origin_b);
	int n = static_cast<int>(A.size());
	for (int i = 0; i < n; ++i) {
		int k = i;
		for (int j = i + 1; j < n; ++j)
			if (abs(A[k][i]) < abs(A[j][i])) k = j;
		swap(A[k], A[i]);
		swap(b[k], b[i]);
		for (int j = i + 1; j < n; ++j) {
			A[j][i] /= A[i][i];
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[j][i] * A[i][k];
			b[j] -= A[j][i] * b[i];
		}
	}
	for (int i = n - 1; i >= 0; --i) {
		for (int j = i + 1; j < n; ++j)
			b[i] -= A[i][j] * b[j];
		b[i] /= A[i][i];
	}
	return b;
}

double get_norm(const vector& vec) {
	double norm = 0;
	for (double x : vec)
		norm = max(norm, abs(x));
	return norm;
}

vector Gauss_Seidel(const matrix& origin_A, const vector& origin_b, double eps) {
	matrix r;
	matrix A(origin_A);
	vector b(origin_b);
	int n = static_cast<int>(b.size());
	for (int i = 0; i < n; ++i) {
		vector tmp;
		for (int j = 0; j < n; ++j)
			tmp.emplace_back(-A[i][j] / A[i][i]);
		r.emplace_back(move(tmp));
		r[i][i] = 0;
		b[i] /= A[i][i];
	}
	vector x0(b), delta(n);
	do {
		vector x1(n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < i; ++j)
				x1[i] += r[i][j] * x1[j];
			for (int j = i + 1; j < n; ++j)
				x1[i] += r[i][j] * x0[j];
			x1[i] += b[i];
		}
		for (int i = 0; i < n; ++i)
			delta[i] = x1[i] - x0[i];
		x0 = x1;
	} while (get_norm(delta) > eps);
	return x0;
}

vector get_precise(double eps) {
	vector ans;
	for (int i = 1; i <= n - 1; ++i) {
		double x = i * h;
		ans.emplace_back((1 - a) / (1 - exp(-1 / eps))
			* (1 - exp(-x / eps)) + a * x);
	}
	return ans;
}

void refresh_matrix(double eps) {
	// 刷新系数矩阵 A
	A.clear();
	vector tmp1(n - 1);
	tmp1[0] = -(2 * eps + h);
	tmp1[1] = eps + h;
	A.emplace_back(move(tmp1));
	for (size_t i = 1; i < (size_t)n - 2; ++i) {
		vector tmp(n - 1);
		tmp[i - 1] = eps;
		tmp[i] = -(2 * eps + h);
		tmp[i + 1] = eps + h;
		A.emplace_back(move(tmp));
	}
	vector tmp2(n - 1);
	tmp2[(size_t)n - 3] = eps;
	tmp2[(size_t)n - 2] = -(2 * eps + h);
	A.emplace_back(move(tmp2));

	// 刷新矩阵 b
	b.clear();
	b.emplace_back(0);
	for (int i = 1; i < n - 2; ++i)
		b.emplace_back(a * h * h);
	b.emplace_back(a * h * h - (eps + h));
}

void display_vector(const vector& vec) {
	cout << '(';
	for (int i = 0; i < vec.size() - 1; ++i)
		cout << vec[i] << ", ";
	cout << vec[vec.size() - 1] << ')' << endl;
}

int main() {
	int testCnt;
	ifstream fin("Homework2.txt");
	fin >> testCnt;
	cout.precision(4);
	for (int i = 0; i < testCnt; ++i) {
		double eps;
		fin >> eps;
		cout << "eps = " << eps << endl;
		auto y_precise = get_precise(eps);
		cout << "精确解：\ny = ";
		display_vector(y_precise);
		refresh_matrix(eps);
		auto y_g = Gauss(A, b);
		cout << "列主元法：\ny = ";
		display_vector(y_g);
		vector delta;
		for (int i = 0; i < n - 1; ++i)
			delta.emplace_back(y_g[i] - y_precise[i]);
		cout << "误差范数为 " << get_norm(delta) << ".\n";
		auto y_gs = Gauss_Seidel(A, b, Eps);
		cout << "Gauss-Seidel 消元法：\ny = ";
		display_vector(y_gs);
		for (int i = 0; i < n - 1; ++i)
			delta[i] = y_gs[i] - y_precise[i];
		cout << "误差范数为 " << get_norm(delta) << ".\n";
		cout << endl;
	}
	fin.close();
	return 0;
}
