#include "Fourier.h"

#include <format>
#include <iostream>
#include <numbers>
#include <random>

using std::numbers::pi;

double f1(double x) {
	return 0.7 * std::sin(4 * pi * x) + std::sin(10 * pi * x);
}

double f2(double x) {
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> dis(0, 1);
	return f1(x) + 0.3 * dis(gen);
}

int main() {
	double (*func[])(double) = { f1, f1, f2 };
	int n[] = { 1 << 4, 1 << 7, 1 << 7 };
	std::string name[] = { "f1", "f1", "f2" };
	for (int i = 0; i < 3; ++i) {
		std::cout << std::format("f = {0}, n = {1}\n", name[i], n[i]);
		auto f = discretization(func[i], n[i]);
		std::cout << "After discretization:\n";
		for (const auto& complex : f)
			std::cout << std::format("{0:.4f} ", complex.real());
		auto g = FFT(f);
		std::cout << "\nReal of g:\n";
		for (const auto& complex : g) 
			std::cout << std::format("{0:.4f} ", complex.real());
		std::cout << "\nImage of g:\n";
		for (const auto& complex : g)
			std::cout << std::format("{0:.4f} ", complex.imag());
		std::cout << "\nMod of g:\n";
		for (const auto& complex : g)
			std::cout << std::format("{0:.4f} ", std::abs(complex));
		auto ifft = IFFT(g);
		std::cout << "\nAfter IFFT:\n";
		for (const auto& complex : ifft)
			std::cout << std::format("{0:.4f} ", complex.real());
		if (i == 2) {
			auto ifft_strip = IFFT(strip_frequency_domain(ifft, 0.25));
			std::cout << "\nIFFT after stripping:\n";
			for (const auto& complex : ifft_strip)
				std::cout << std::format("{0:.4f} ", complex.real());
		}
		std::cout << '\n' << std::endl;
	}
	return 0;
}
