#pragma once

#include <string>
#include <vector>
#include <ostream>

class Vector {
	std::vector<double> dat;
	size_t len = 0;
public:
	Vector() = default;
	Vector(size_t size) : dat(size), len(size) { }
	Vector(const std::vector<double>& vec) : dat(vec), len(vec.size()) { }
	Vector(const std::initializer_list<double>& list) : dat(list), len(list.size()) { }
	
	double& operator[](size_t index) const;
	Vector operator+(const Vector& another) const;
	Vector operator-(const Vector& another) const;
	Vector operator-() const;
	double operator*(const Vector& another) const;
	Vector operator*(double num) const;
	Vector& operator*=(double num);
	Vector& operator+=(const Vector& another);
	Vector& operator-=(const Vector& another);
	friend Vector operator*(double num, const Vector& a);
	
	std::string to_string() const;
	size_t length() const { return len; }
	const auto begin() const { return dat.begin(); }
	const auto end() const { return dat.end(); }

	double mod() const;
	Vector& unitization();
};
