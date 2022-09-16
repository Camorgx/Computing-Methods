#include <iostream>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <vector>
#include <cmath>

#include "Homework2.h"

using std::move;
using std::vector;
using std::swap, std::max;
using std::abs, std::max;
using std::cout, std::endl;
using std::fixed, std::setprecision;
using std::ifstream;

vector<double> Homework2::Gauss(double eps) {
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

vector<double> Homework2::Gauss_Seidel(double eps) {
	vector<vector<double>> r;
	for (int i = 0; i < n; ++i) {
		vector<double> tmp;
		for (int j = 0; j < n; ++j)
			tmp.emplace_back(-A[i][j] / A[i][i]);
		r.emplace_back(move(tmp));
		r[i][i] = 0;
		b[i] /= A[i][i];
	}
	vector<double> x0(b), delta(n);
	do {
		vector<double> x1(n);
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
	} while (get_norm(delta) > Eps);
	return x0;
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
	vector<double> tmp1(n);
	tmp1[0] = -(2 * eps + 1.0 / n);
	tmp1[1] = eps + 1.0 / n;
	A.emplace_back(move(tmp1));
	for (size_t i = 1; i < (size_t)n - 1; ++i) {
		vector<double> tmp(n);
		tmp[i - 1] = eps;
		tmp[i] = -(2 * eps + 1.0 / n);
		tmp[i + 1] = eps + 1.0 / n;
		A.emplace_back(move(tmp));
	}
	vector<double> tmp2(n);
	tmp2[(size_t)n - 2] = eps;
	tmp2[(size_t)n - 1] = -(2 * eps + 1.0 / n);
	A.emplace_back(tmp2);

	// 刷新矩阵 B
	b.clear();
	for (int i = 0; i < n; ++i)
		b.emplace_back(a / (n * n));
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
		refresh_matrix(eps);
		auto x_g = Gauss(eps);
		cout << "eps = " << eps << endl;
		cout << "列主元法：\nx = ";
		display_vector(x_g);
		refresh_matrix(eps);
		auto x_gs = Gauss_Seidel(eps);
		cout << "Gauss-Seidel 消元法：\nx = ";
		display_vector(x_gs);
		for (int i = 0; i < n; ++i)
			x_gs[i] -= x_g[i];
		cout << "误差：\n";
		display_vector(x_gs);
		cout << "误差范数：\n";
		cout << get_norm(x_gs) << endl;
		cout << endl;
	}
}
