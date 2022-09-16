#pragma once

#include "ITestable.h"

#include <vector>

class Homework2 : public ITestable {
	const double Eps = 1e-6;
	const double a = 0.5;
	const int n = 100;

	std::vector<std::vector<double>> A;
	std::vector<double> b;

	std::vector<double> Gauss(double eps);
	std::vector<double> Gauss_Seidel(double eps);
	void refresh_matrix(double eps);
	void display_vector(const std::vector<double>& vec);
	double get_norm(const std::vector<double>& vec);
public:
	void run_test();
};
