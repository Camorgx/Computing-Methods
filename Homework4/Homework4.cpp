#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <format>
#include <iostream>

Vector get_eigenvector(const Matrix& origin_A, double eigenvalue) {
	Matrix A(origin_A);
	int n = static_cast<int>(A.column_size());
	Matrix I = Matrix::eye(n);
	A -= eigenvalue * I;
	for (int i = 0; i < n; ++i) {
		int k = i;
		for (int j = i + 1; j < n; ++j)
			if (std::abs(A[k][i]) < std::abs(A[j][i])) 
				k = j;
		if (k != i) std::swap(A[k], A[i]);
		if (A[i][i] == 0) continue;
		A[i] *= 1 / A[i][i];
		for (int j = i + 1; j < n; ++j) {
			A[j] -= A[j][i] * A[i];
		}
	}
	Vector b(n);
	b[static_cast<size_t>(n) - 1] = 1;
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
	double sum = 0;
	while ((sum = item_square_sum(a)) > eps) {
		std::cout << std::format("The sum of squares of"
			"non-diagonal elements: {0:.4f}", sum) << std::endl;
		// 选取非对角线按模最大元素
		double max_mod = 0;
		int p = 0, q = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				double mod = std::abs(a[i][j]);
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
			if (std::abs(t1) > std::abs(t2))
				t = t2;
			else t = t1;
		}
		double c = 1 / std::sqrt(1 + t * t);
		double d = t / std::sqrt(1 + t * t);
		// 迭代计算
		Matrix b(a);
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
	}
	std::vector<double> res;
	for (int i = 0; i < n; ++i)
		res.push_back(a[i][i]);
	return res;
}

void SVD(const Matrix& origin_a, Matrix& u, Matrix& v, Matrix& omega, const double eps) {
	Matrix a(origin_a); // m * n
	Matrix at = a.transpose(); // n * m
	Matrix aat = a * at; // m * m
	Matrix ata = at * a; // n * n
	size_t m = a.row_size(), n = a.column_size();
	const Matrix I1 = Matrix::eye(m);
	const Matrix I2 = Matrix::eye(n);
	const Vector zero1(m);
	const Vector zero2(n);
	auto lambda = Jacobi(aat, eps);
	std::sort(lambda.begin(), lambda.end(), std::greater<double>());
	int append_count = static_cast<int>(m - n);
	append_count = std::abs(append_count);
	for (size_t i = 0; i < append_count; ++i)
		lambda.push_back(0);
	u = Matrix(m, m);
	for (int i = 0; i < m; ++i)
		u[i] = get_eigenvector(aat, lambda[i]).unitization();
	omega = Matrix(m, n);
	n = std::min(m, n);
	for (size_t i = 0; i < n; ++i)
		omega[i][i] = std::sqrt(lambda[i]);
	v = Matrix(n, n);
	for (int i = 0; i < n; ++i) {
		if (lambda[i] == 0)
			v[i] = (at * u[i]).unitization();
		else v[i] = 1 / omega[i][i] * (at * u[i]);
	}
	u = u.transpose();
	v = v.transpose();
}

// 返回数组的第一维表示基，第二维是对应的坐标
auto PCA(std::vector<Vector>& data, double eps, int k) {
	std::vector<std::vector<double>> res;
	Vector sum(data[0].length());
	size_t data_size = data.size();
	for (const auto& dat : data)
		sum += dat;
	sum *= 1.0 / data_size;
	for (auto& dat : data)
		dat -= sum;
	Matrix x(data);
	x = x.transpose();
	Matrix covariance = 1.0 / x.column_size() * x * x.transpose();
	std::cout << "Covariance matrix: " << covariance.to_string() << std::endl;
	auto eigenvalues = Jacobi(covariance, eps);
	std::sort(eigenvalues.begin(), eigenvalues.end(), std::greater<double>());
	std::vector<Vector> e;
	for (int i = 0; i < k; ++i) {
		e.emplace_back(get_eigenvector(covariance, eigenvalues[i]).unitization());
	}
	for (int i = 0; i < k; ++i) {
		res.emplace_back(data_size);
		for (size_t j = 0; j < data_size; ++j) {
			res[i][j] = data[j] * e[i];
		}
	}
	return res;
}

int main() {
	const double eps = 1e-6;

	Matrix A = Matrix::random_initialize(4, 3);
	std::cout << "Random matrix A generated: " << A.to_string() << std::endl;
	auto eigen = Jacobi(A * A.transpose(), eps);
	std::cout << "Eigenvalues of A: ";
	for (double value : eigen)
		std::cout << std::format("{0:.4f} ", value);
	std::cout << std::endl;
	Matrix u, omega, v;
	SVD(A, u, v, omega, eps);
	std::cout << "u = " << u.to_string() << std::endl;
	std::cout << "vt = " << v.transpose().to_string() << std::endl;
	std::cout << "omega = " << omega.to_string() << std::endl;
	std::cout << std::endl;

	std::FILE* fin;
	errno_t err = fopen_s(&fin, "iris.txt", "r");
	if (!fin) {
		std::cout << "Failed to open file iris.txt.\n";
		return 0;
	}
	std::vector<Vector> data;
	for (int j = 0; j < 150; ++j) {
		double dat[4];
		int id;
		fscanf_s(fin, "%lf,%lf,%lf,%lf,%d",
			dat, dat + 1, dat + 2, dat + 3, &id);
		data.emplace_back(Vector { dat[0], dat[1], dat[2], dat[3] });
	}
	std::fclose(fin);
	auto res = PCA(data, eps, 2);
	std::string output_x = "x = [";
	std::string output_y = "y = [";
	for (int i = 0; i < 149; ++i) {
		output_x += std::format("{0:.4f}, ", res[0][i]);
		output_y += std::format("{0:.4f}, ", res[1][i]);
	}
	output_x += std::format("{0:.4f}]", res[0][149]);
	output_y += std::format("{0:.4f}]", res[1][149]);
	std::cout << output_x << std::endl;
	std::cout << output_y << std::endl;
    return 0;
}
