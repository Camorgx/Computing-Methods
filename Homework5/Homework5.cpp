#include <exception>
#include <fstream>
#include <format>
#include <iostream>
#include <vector>

using vector = std::vector<double>;

// 使用自然边界条件
void generate_euqations(const vector& x, const vector& y,
	vector& mu, vector& lambda, vector& d) {
	size_t n = x.size() - 1;
	vector h(n);
	mu.resize(n);
	lambda.resize(n);
	d.resize(n);
	for (size_t i = 0; i < n; ++i)
		h[i] = x[i + 1] - x[i];
	for (size_t i = 1; i < n; ++i) {
		lambda[i] = h[i] / (h[i] + h[i - 1]);
		mu[i] = 1 - lambda[i];
	}
	for (size_t i = 1; i < n; ++i) {
		d[i] = (6 / (h[i] + h[i - 1]))
			* ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
	}
}

// 追赶法解三对角方程组
// a[1 .. n], b[1 .. n - 1], c[2 .. n] 
vector Thomas_solve(const vector& a, const vector& b, const vector& c, const vector& f) {
	size_t n = a.size() - 1;
	vector x(n + 1), u(n + 1), v(n + 1), y(n + 1);
	for (size_t k = 1; k <= n; ++k) {
		u[k] = a[k] - c[k] * v[k - 1];
		v[k] = b[k] / u[k];
		y[k] = (f[k] - c[k] * y[k - 1]) / u[k];
	}
	v[n] = 0;
	for (size_t k = n; k >= 1; --k)
		x[k - 1] = y[k] - v[k] * x[k];
	x.pop_back();
	return x;
}

// 使用 M 方法进行三次样条插值
vector interpolation(const vector& x, const vector& y) {
	size_t n = x.size() - 1;
	vector M(n + 1);
	vector mu, lambda, d;
	generate_euqations(x, y, mu, lambda, d);
	mu[1] = 0;
	vector a(n, 2);
	return Thomas_solve(a, lambda, mu, d);
}

int main() {
	const int n = 20;
	vector x(n + 1), y(n + 1);
	std::ifstream fin("point.txt");
	for (int i = 0; i <= n; ++i)
		fin >> x[i] >> y[i];
	fin.close();
	auto res1 = interpolation(x, y);
	x[9] = 10;
	auto res2 = interpolation(x, y);
	for (int i = 0; i < n - 1; ++i) {
		std::cout <<
			std::format("origin[{0}] = {1:.4f}, later[{0}] = {2:.4f}\n", 
				i + 1, res1[i], res2[i]);
	}
	std::cout << std::endl;
	return 0;
}
