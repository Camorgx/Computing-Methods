#include <iostream>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <vector>
#include <cmath>

#include "Homework2.h"

using std::move;
using std::vector;
using std::swap;
using std::abs, std::max, std::exp;
using std::cout, std::endl;
using std::fixed, std::setprecision;
using std::ifstream;

vector<double> Homework2::Gauss(double eps) {
	for (int i = 0; i < n - 1; ++i) {
		int k = i;
		for (int j = i + 1; j < n - 1; ++j)
			if (abs(A[k][i]) < abs(A[j][i])) k = j;
		swap(A[k], A[i]);
		swap(b[k], b[i]);
		for (int j = i + 1; j < n - 1; ++j) {
			A[j][i] /= A[i][i];
			for (int k = i + 1; k < n - 1; ++k)
				A[j][k] -= A[j][i] * A[i][k];
			b[j] -= A[j][i] * b[i];
		}
	}
	for (int i = n - 2; i >= 0; --i) {
		for (int j = i + 1; j < n - 1; ++j)
			b[i] -= A[i][j] * b[j];
		b[i] /= A[i][i];
	}
	return b;
}

vector<double> Homework2::Gauss_Seidel(double eps) {
	vector<vector<double>> r;
	for (int i = 0; i < n - 1; ++i) {
		vector<double> tmp;
		for (int j = 0; j < n - 1; ++j)
			tmp.emplace_back(-A[i][j] / A[i][i]);
		r.emplace_back(move(tmp));
		r[i][i] = 0;
		b[i] /= A[i][i];
	}
	vector<double> x0(b), delta(n - 1);
	do {
		vector<double> x1(n - 1);
		for (int i = 0; i < n - 1; ++i) {
			for (int j = 0; j < i; ++j)
				x1[i] += r[i][j] * x1[j];
			for (int j = i + 1; j < n - 1; ++j)
				x1[i] += r[i][j] * x0[j];
			x1[i] += b[i];
		}
		for (int i = 0; i < n - 1; ++i)
			delta[i] = x1[i] - x0[i];
		x0 = x1;
	} while (get_norm(delta) > Eps);
	return x0;
}

vector<double> Homework2::get_precise(double eps) {
	vector<double> ans;
	for (int i = 1; i <= n - 1; ++i) {
		double x = i * h;
		ans.emplace_back((1 - a) / (1 - exp(-1 / eps))
			* (1 - exp(-x / eps)) + a * x);
	}
	return ans;
}

double Homework2::get_norm(const vector<double>& vec) {
	double norm = 0;
	for (double x : vec)
		norm = max(norm, abs(x));
	return norm;
}

void Homework2::refresh_matrix(double eps) {
	// 刷新系数矩阵 A
	A.clear();
	vector<double> tmp1(n - 1);
	tmp1[0] = -(2 * eps + h);
	tmp1[1] = eps + h;
	A.emplace_back(move(tmp1));
	for (size_t i = 1; i < (size_t)n - 2; ++i) {
		vector<double> tmp(n - 1);
		tmp[i - 1] = eps;
		tmp[i] = -(2 * eps + h);
		tmp[i + 1] = eps + h;
		A.emplace_back(move(tmp));
	}
	vector<double> tmp2(n - 1);
	tmp2[(size_t)n - 3] = eps;
	tmp2[(size_t)n - 2] = -(2 * eps + h);
	A.emplace_back(tmp2);

	// 刷新矩阵 B
	b.clear();
	b.emplace_back(0);
	for (int i = 1; i < n - 2; ++i)
		b.emplace_back(a * h * h);
	b.emplace_back(a * h * h - (eps + h));
}

void Homework2::display_vector(const std::vector<double>& vec) {
	cout << '(';
	for (int i = 0; i < vec.size() - 1; ++i)
		cout << vec[i] << ", ";
	cout << vec[vec.size() - 1] << ')' << endl;
}

void Homework2::run_test() {
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
		auto y_g = Gauss(eps);
		cout << "列主元法：\ny = ";
		display_vector(y_g);
		refresh_matrix(eps);
		vector<double> delta;
		for (int i = 0; i < n - 1; ++i)
			delta.emplace_back(y_g[i] - y_precise[i]);
		cout << "误差范数为 " << get_norm(delta) << ".\n";
		auto y_gs = Gauss_Seidel(eps);
		cout << "Gauss-Seidel 消元法：\ny = ";
		display_vector(y_gs);
		for (int i = 0; i < n - 1; ++i)
			delta[i] = y_gs[i] - y_precise[i];
		cout << "误差范数为 " << get_norm(delta) << ".\n";
		cout << endl;
	}
	fin.close();
}
