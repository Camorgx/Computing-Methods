#pragma once

#include "ITestable.h"

class Homework1 : public ITestable {
	const double eps = 1e-8;
	void solve(double p, double x, double y);
public:
	void run_test();
};
