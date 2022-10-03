#pragma once

#include "ITestable.h"

#include <vector>

class Homework2 : public ITestable {
	const double Eps = 1e-6;
	const double a = 0.5;
	const int n = 100;
	const double h = 1.0 / n;

	std::vector<std::vector<double>> A;
	std::vector<double> b;
	
	std::vector<double> get_precise(double eps);
	void refresh_matrix(double eps);
public:
	using vector = std::vector<double>;
	using matrix = std::vector<std::vector<double>>;
	static void display_vector(const vector& vec);
	// get the infinity norm of vec
	static double get_norm(const vector& vec);
	static vector Gauss(const matrix& A, const vector& b);
	static vector Gauss_Seidel(const matrix& origin_A, 
		const vector& origin_b, double eps);
	void run_test();
};
