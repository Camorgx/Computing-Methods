#pragma once

#include "Vector.h"

class Matrix {
	std::vector<Vector> dat;
	size_t rowCnt = 0, columnCnt = 0;
public:
	Matrix() = default;
	Matrix(size_t rowSize, size_t columnSize);
	Matrix(const std::initializer_list<Vector>& list);
	Matrix(const std::initializer_list<std::initializer_list<double>>& list);
	Matrix(const std::vector<Vector>& data);

	Vector& operator[](size_t index) const;
	Matrix operator-() const;
	Matrix operator+(const Matrix& another) const;
	Matrix operator-(const Matrix& another) const;
	Matrix operator*(const Matrix& another) const;
	Vector operator*(const Vector& another) const;
	Matrix operator*(double num) const;
	Matrix& operator+=(const Matrix& another);
	Matrix& operator-=(const Matrix& another);
	Matrix& operator*=(Matrix& another);
	Matrix& operator*=(double num);
	friend Matrix operator*(double num, const Matrix& a);

	size_t row_size() const { return rowCnt; }
	size_t column_size() const { return columnCnt; }
	Matrix transpose() const;
	std::string to_string() const;
	const auto begin() const { return dat.begin(); }
	const auto end() const { return dat.end(); }

	static Matrix eye(size_t size);
};
