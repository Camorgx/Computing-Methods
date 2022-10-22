#include "Vector.h"
#include <format>
#include <exception>

double& Vector::operator[](size_t index) {
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

const Vector& Vector::operator+=(const Vector& another) {
	if (len != another.len)
		throw std::runtime_error(std::format(
			"Invalid vector length: {0} and {1}.", len, another.len));
	for (int i = 0; i < len; ++i)
		dat[i] += another.dat[i];
	return *this;
}

const Vector& Vector::operator-=(const Vector& another) {
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
