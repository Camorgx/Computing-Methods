#include "Vector.h"

#include <cmath>
#include <format>
#include <exception>

double& Vector::operator[](size_t index) {
	return dat.at(index);
}

const double& Vector::operator[](size_t index) const {
	return dat.at(index);
}

Vector Vector::operator+(const Vector& another) const {
	Vector res(*this);
	res += another;
	return res;
}

Vector Vector::operator-(const Vector& another) const {
	Vector res(*this);
	res -= another;
	return res;
}

Vector Vector::operator-() const {
	Vector res;
	res.len = len;
	for (const double& x : dat)
		res.dat.push_back(-x);
	return res;
}

double Vector::operator*(const Vector& another) const {
	if (len != another.len)
		throw std::runtime_error(std::format(
			"Invalid vector length: {0} and {1}.", len, another.len));
	double res = 0;
	for (int i = 0; i < len; ++i)
		res += dat[i] * another.dat[i];
	return res;
}

Vector Vector::operator*(double num) const {
	Vector res(*this);
	res *= num;
	return res;
}

Vector& Vector::operator*=(double num) {
	for (double& x : dat)
		x *= num;
	return *this;
}

Vector& Vector::operator+=(const Vector& another) {
	if (len != another.len)
		throw std::runtime_error(std::format(
			"Invalid vector length: {0} and {1}.", len, another.len));
	for (int i = 0; i < len; ++i)
		dat[i] += another.dat[i];
	return *this;
}

Vector& Vector::operator-=(const Vector& another) {
	if (len != another.len)
		throw std::runtime_error(std::format(
			"Invalid vector length: {0} and {1}.", len, another.len));
	for (int i = 0; i < len; ++i)
		dat[i] -= another.dat[i];
	return *this;
}

std::string Vector::to_string() const {
	if (len == 0) return "()";
	std::string res = "(";
	for (int i = 0; i < len - 1; ++i)
		res += std::format("{0:.4f}, ", dat[i]);
	res += std::format("{0:.4f})", dat.back());
	return res;
}

double Vector::mod() const {
	double res = 0;
	for (double x : dat)
		res += x * x;
	return std::sqrt(res);
}

Vector& Vector::unitization() {
	double square = mod();
	if (square != 0)
		*this *= 1 / mod();
	return *this;
}

Vector operator*(double num, const Vector& a) {
	Vector res(a);
	res *= num;
	return res;
}
