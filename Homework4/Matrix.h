#pragma once

#include "Vector.h"

class Matrix {
	std::vector<Vector> dat;
	size_t rowCnt = 0, columnCnt = 0;
public:
	Matrix(size_t rowSize, size_t columnSize);
	Matrix(std::initializer_list<Vector> list);
	Matrix(std::initializer_list<std::initializer_list<double>> list);

	Vector& operator[](size_t index);
	Matrix operator-() const;
	Matrix operator+(const Matrix& another) const;
	Matrix operator-(const Matrix& another) const;
	Matrix operator*(Matrix& another);
	const Matrix& operator+=(const Matrix& another);
	const Matrix& operator-=(const Matrix& another);
	const Matrix& operator*=(Matrix& another);

	size_t row_size() const { return rowCnt; }
	size_t column_size() const { return columnCnt; }
	Matrix transpose();
	std::string to_string() const;
	const auto begin() const { return dat.begin(); }
	const auto end() const { return dat.end(); }
};
