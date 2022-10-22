#include "Matrix.h"

#include <cmath>
#include <exception>
#include <iostream>

Vector Gauss(const Matrix& origin_A, const Vector& origin_b) {
	Matrix A(origin_A);
	Vector b(origin_b);
	if (A.column_size() != A.row_size())
		throw std::invalid_argument("Gauss Algorithm requires square matrix.");
	int n = static_cast<int>(A.column_size());
	for (int i = 0; i < n; ++i) {
		int k = i;
		for (int j = i + 1; j < n; ++j)
			if (abs(A[k][i]) < abs(A[j][i])) k = j;
		std::swap(A[k], A[i]);
		std::swap(b[k], b[i]);
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

double item_square_sum(Matrix& A) {
	double res = 0;
	int n = static_cast<int>(A.column_size());
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i != j)
				res += A[i][j] * A[i][j];
	return res;
}

std::vector<double> Jacobi(const Matrix& origin_A, const double eps) {
	Matrix a(origin_A);
	if (a.column_size() != a.row_size())
		throw std::invalid_argument("Jacobi Algorithm requires square matrix.");
	int n = static_cast<int>(a.column_size());
	while (item_square_sum(a) > eps) {
		// 选取对角线按模最大元素
		double max_mod = 0;
		int p = 0, q = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				double mod = std::abs(a[i][j]);
				//std::cout << std::format("mod = {0}", mod) << std::endl;
				if (i != j && mod > max_mod) {
					max_mod = mod;
					p = i;
					q = j;
				}
			}
		}
		// 确定旋转角度
		double s = (a[q][q] - a[p][p]) / (2 * a[p][q]);
		double t = 0;
		if (std::abs(s) < eps) t = 1;
		else {
			double t1 = -s - std::sqrt(s * s + 1);
			double t2 = -s + std::sqrt(s * s + 1);
			t = std::min(std::abs(t1), std::abs(t2));
		}
		double c = 1 / std::sqrt(1 + t * t);
		double d = t / std::sqrt(1 + t * t);
		// 迭代计算
		Matrix b(a);
		std::cout << "b = \n";
		std::cout << b.to_string() << std::endl;
		for (int i = 0; i < n; ++i) {
			if (i != p && i != q) {
				b[i][p] = b[p][i] = c * a[p][i] - d * a[q][i];
				b[i][q] = b[q][i] = c * a[q][i] + d * a[p][i];
			}
		}
		b[p][p] = a[p][p] -  t * a[p][q];
		b[q][q] = a[q][q] + t * a[p][q];
		b[p][q] = b[q][p] = 0;
		a = b;
		std::cout << "a = \n";
		std::cout << a.to_string() << std::endl;
	}
	std::vector<double> res;
	for (int i = 0; i < n; ++i)
		res.push_back(a[i][i]);
	return res;
}

int main() {
	Matrix m1 = { { 1, 2, 3 }, { 2, 2, 0 }, { 3, 0, 3 } };
	auto res1 = Jacobi(m1, 1e-6);
	for (double x : res1)
		std::cout << x << ' ';
	std::cout << std::endl;
    return 0;
}
