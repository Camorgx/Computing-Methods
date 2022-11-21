#include "Fourier.h"

#include <numbers>

complex_vector discretization(double (*f)(double x), int n) {
	complex_vector res;
	for (int i = 0; i < n; ++i)
		res.emplace_back(f(static_cast<double>(i) / n));
	return res;
}

using std::numbers::pi;
using complex = std::complex<double>;
using namespace std::complex_literals;

complex_vector FFT(const complex_vector& f) {
	int n = static_cast<int>(f.size());
	complex_vector y(n);
	if (n == 1) return f;
	complex wn = std::exp(-2 * pi * 1i / static_cast<double>(n));
	complex w = 1;
	complex_vector f0, f1;
	for (int i = 0; i < n; i += 2)
		f0.emplace_back(f[i]);
	for (int i = 1; i < n; i += 2)
		f1.emplace_back(f[i]);
	auto y0 = FFT(f0);
	auto y1 = FFT(f1);
	int half_n = n >> 1;
	for (int k = 0; k < half_n; ++k) {
		y[k] = (y0[k] + w * y1[k]) / 2.0;
		y[static_cast<size_t>(k) + half_n] = (y0[k] - w * y1[k]) / 2.0;
		w *= wn;
	}
	return y;
}

complex_vector IFFT(const complex_vector& f) {
	int n = static_cast<int>(f.size());
	complex_vector y(n);
	if (n == 1) return f;
	complex wn = std::exp(2 * pi * 1i / static_cast<double>(n));
	complex w = 1;
	complex_vector f0, f1;
	for (int i = 0; i < n; i += 2)
		f0.emplace_back(f[i]);
	for (int i = 1; i < n; i += 2)
		f1.emplace_back(f[i]);
	auto y0 = IFFT(f0);
	auto y1 = IFFT(f1);
	int half_n = n >> 1;
	for (int k = 0; k < half_n; ++k) {
		y[k] = y0[k] + w * y1[k];
		y[static_cast<size_t>(k) + half_n] = y0[k] - w * y1[k];
		w *= wn;
	}
	return y;
}

complex_vector strip_frequency_domain(const complex_vector& f, double proportion) {
	complex_vector res(f);
	for (int i = static_cast<int>(f.size() * proportion); i < f.size(); ++i)
		res[i] = 0;
	return res;
}
