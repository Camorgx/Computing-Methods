#pragma once

#include "ITestable.h"

class Homework2 : public ITestable {
	void Gauss(double eps);
	void Gauss_Seidel(double eps);
public:
	void run_test() override;
};
