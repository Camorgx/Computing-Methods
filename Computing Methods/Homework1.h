#pragma once

#include "ITestable.h"

class Homework1 : public ITestable {
	void solve(double p, double x, double y);
public:
	void run_test();
};
