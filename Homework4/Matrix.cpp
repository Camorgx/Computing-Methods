#include "Matrix.h"

#include <exception>
#include <random>

Matrix::Matrix(size_t row, size_t column) {
	rowCnt = row; 
	columnCnt = column;
	for (size_t i = 0; i < row; ++i)
		dat.emplace_back(column);
}

Matrix::Matrix(const std::initializer_list<Vector>& list) {
	rowCnt = list.size();
	bool first = true;
	size_t prev = 0;
	for (const auto& l : list) {
		if (first) {
			prev = l.length();
			columnCnt = l.length();
			first = false;
		}
		else if (l.length() != prev)
			throw std::invalid_argument("Invalid initializer for matrix.");
		dat.emplace_back(l);
		prev = l.length();
	}
}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& list) {
	rowCnt = list.size();
	bool first = true;
	size_t prev = 0;
	for (const auto& l : list) {
		if (first) {
			prev = l.size();
			columnCnt = prev;
			first = false;
		}
		else if (l.size() != prev)
			throw std::invalid_argument("Invalid initializer for matrix.");
		dat.emplace_back(l);
		prev = l.size();
	}
}

Matrix::Matrix(const std::vector<Vector>& data) {
	rowCnt = data.size();
	columnCnt = data[0].length();
	for (size_t i = 0; i < rowCnt; ++i) {
		if (data[i].length() != columnCnt)
			throw std::invalid_argument("Invalid initializer for matrix.");
		dat.emplace_back(data[i]);
	}
}

Vector& Matrix::operator[](size_t index) const {
	return const_cast<Vector&>(dat.at(index));
}

Matrix Matrix::operator-() const {
	Matrix res(*this);
	for (auto& x : res.dat)
		x = -x;
	return res;
}

Matrix Matrix::operator+(const Matrix& another) const {
	Matrix res(*this);
	res += another;
	return res;
}

Matrix Matrix::operator-(const Matrix& another) const {
	Matrix res(*this);
	res -= another;
	return res;
}

Matrix Matrix::operator*(const Matrix& another) const {
	if (columnCnt != another.rowCnt)
		throw std::runtime_error("Invliad size for multiply operation.");
	Matrix res(rowCnt, another.columnCnt);
	for (size_t i = 0; i < rowCnt; ++i) {
		for (size_t j = 0; j < another.columnCnt; ++j) {
			for (size_t k = 0; k < columnCnt; ++k) {
				res[i][j] += dat[i][k] * another.dat[k][j];
			}
		}
	}
	return res;
}

Vector Matrix::operator*(const Vector& another) const {
	size_t n = another.length();
	if (columnCnt != n)
		throw std::runtime_error("Invliad size for multiply operation.");
	Vector res(n);
	for (size_t i = 0; i < rowCnt; ++i)
		res[i] = dat[i] * another;
	return res;
}

Matrix Matrix::operator*(double num) const {
	Matrix res(*this);
	res *= num;
	return res;
}

Matrix& Matrix::operator+=(const Matrix& another) {
	if (!(rowCnt == another.rowCnt && columnCnt == another.columnCnt))
		throw std::runtime_error("Invliad size for add operation.");
	for (size_t i = 0; i < rowCnt; ++i)
		dat[i] += another.dat[i];
	return *this;
}

Matrix& Matrix::operator-=(const Matrix& another) {
	if (!(rowCnt == another.rowCnt && columnCnt == another.columnCnt))
		throw std::runtime_error("Invliad size for minus operation.");
	for (size_t i = 0; i < rowCnt; ++i)
		dat[i] -= another.dat[i];
	return *this;
}

Matrix& Matrix::operator*=(Matrix& another) {
	*this = *this * another;
	return *this;
}

Matrix& Matrix::operator*=(double num) {
	for (Vector& row : dat)
		row *= num;
	return *this;
}


Matrix Matrix::transpose() const {
	Matrix res(columnCnt, rowCnt);
	for (size_t i = 0; i < rowCnt; ++i) {
		for (size_t j = 0; j < columnCnt; ++j) {
			res[j][i] = dat[i][j];
		}
	}
	return res;
}

std::string Matrix::to_string() const {
	if (rowCnt == 0) return "[]";
	std::string res = "[\n";
	for (size_t i = 0; i < rowCnt; ++i) {
		std::string row = dat[i].to_string();
		row = row.substr(1, row.length() - 2);
		res += row + ";\n";
	}
	res += "]";
	return res;
}

Matrix Matrix::eye(size_t size) {
	Matrix res(size, size);
	for (size_t i = 0; i < size; ++i)
		res[i][i] = 1;
	return res;
}

Matrix Matrix::random_initialize(size_t row, size_t column) {
	Matrix res(row, column);
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> dis(0, 1);
	for (size_t i = 0; i < row; ++i)
		for (size_t j = 0; j < column; ++j)
			res[i][j] = dis(gen);
	return res;
}

Matrix operator*(double num, const Matrix& a) {
	Matrix res(a);
	res *= num;
	return res;
}
