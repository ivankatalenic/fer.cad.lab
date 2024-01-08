#pragma once

#include "matrix.hh"

#include <string>
#include <istream>
#include <ostream>
#include <functional>

namespace MatrixUtils {
	void print(const Matrix& m, std::ostream& o);
	void printToConsole(const Matrix& m);
	void printToFile(const Matrix& m, const std::string& filename);

	Matrix load(std::istream& i);
	Matrix loadFromFile(const std::string& filename);

	// Deprecated
	Matrix abs(const Matrix& m);
	// Deprecated
	bool   is_less(const Matrix& a, const Matrix& b);
	// Deprecated
	bool   is_positive(const Matrix& m);

	bool unary_predicate(const Matrix& m, std::function<bool(double)> predicate);
	bool binary_predicate(const Matrix& a, const Matrix& b, std::function<bool(double, double)> predicate);

	Matrix unary_op(const Matrix& m, std::function<double(double)> op);
	Matrix binary_op(const Matrix& a, const Matrix& b, std::function<double(double, double)> op);

	Matrix identity(const int size);
}
