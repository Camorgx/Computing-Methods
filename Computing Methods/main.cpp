#include "Homeworks.h"
#include <cstdlib>

int main() {
	ITestable* test = new Homework3;
	test->run_test();
	delete test;
	return 0;
}