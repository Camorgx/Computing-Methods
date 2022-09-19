#include "Homeworks.h"

int main() {
	ITestable* test = new Homework2;
	test->run_test();
	delete test;
	return 0;
}