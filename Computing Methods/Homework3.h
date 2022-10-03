#pragma once

#include "ITestable.h"

#include <vector>

class Homework3 : public ITestable {
	const double eps = 1e-5;
	std::vector<std::vector<double>> X, Y;
public:
	using matrix = std::vector<std::vector<double>>;
	using vector = std::vector<double>;
	static void Doolittle_decompose(const matrix& a, matrix& l, matrix& u);
	static matrix init_matrix(size_t m, size_t n);
	static vector Doolittle_solve(const vector& b, const matrix& l, const matrix& u);
	int anti_exponentiation(const matrix& a, double& eigenvalue, vector& eigenvector);
	// get the infinite norm of matrix m
	static double get_norm(const matrix& mat);
	static void display_matrix(const matrix& mat);
	void run_test();
};
