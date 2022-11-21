#pragma once

#include <complex>
#include <vector>

using complex_vector = std::vector<std::complex<double>>;

complex_vector discretization(double (*f)(double x), int n);

complex_vector FFT(const complex_vector& f);

complex_vector IFFT(const complex_vector& f);

complex_vector strip_frequency_domain(const complex_vector& f, double proportion);
