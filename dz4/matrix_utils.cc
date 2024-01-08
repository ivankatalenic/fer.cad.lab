#include "matrix_utils.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

void MatrixUtils::print(const Matrix& m, std::ostream& o) {
	const int rows{m.getRows()};
	const int columns{m.getColumns()};

	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < columns; ++j) {
			o << m[i][j];

			if (j < columns - 1) {
				o << '\t';
			}
		}

		o << '\n';
	}
}

void MatrixUtils::printToConsole(const Matrix& m) {
	print(m, std::cout);
}

void MatrixUtils::printToFile(const Matrix& m, const std::string& filename) {
	std::ofstream file_stream(filename);
	if (!file_stream) {
		throw std::runtime_error("can't write to a file " + filename);
	}

	print(m, file_stream);
}

Matrix MatrixUtils::load(std::istream& i) {
	std::vector<double> values;

	// Determine the number of columns (the number of doubles in the first line)
	constexpr int LINE_SIZE_MAX{512};
	char line[LINE_SIZE_MAX];

	i.getline(line, LINE_SIZE_MAX);

	std::string str(line);
	std::istringstream string_stream(str);
	int columns{0};
	for (double v; string_stream >> v;) {
		++columns;

		values.push_back(v);
	}

	// Extract the remaining values
	for (double v; i >> v;) {
		values.push_back(v);
	}

	int rows{static_cast<int>(values.size()) / columns};

	Matrix result(rows, columns);

	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < columns; ++j) {
			result[i][j] = values[i * columns + j];
		}
	}

	return result;
}

Matrix MatrixUtils::loadFromFile(const std::string& filename) {
	std::ifstream file_stream(filename);
	if (!file_stream) {
		throw std::invalid_argument("can't open the matrix file " + filename);
	}

	return load(file_stream);
}

Matrix MatrixUtils::abs(const Matrix& m) {
	Matrix ret(m);
	for (int i{0}; i < m.getRows(); ++i) {
		for (int j{0}; j < m.getColumns(); ++j) {
			ret[i][j] = std::abs(ret[i][j]);
		}
	}

	return ret;
}

bool MatrixUtils::is_less(const Matrix& a, const Matrix& b) {
	if (a.getRows() != b.getRows()) {
		throw std::invalid_argument("the matrices don't have the same number of rows");
	}
	if (a.getColumns() != b.getColumns()) {
		throw std::invalid_argument("the matrices don't have the same number of columns");
	}

	for (int i{0}; i < a.getRows(); ++i) {
		for (int j{0}; j < a.getColumns(); ++j) {
			if (a[i][j] >= b[i][j]) {
				return false;
			}
		}
	}

	return true;
}

bool MatrixUtils::is_positive(const Matrix& m) {
	for (int i{0}; i < m.getRows(); ++i) {
		for (int j{0}; j < m.getColumns(); ++j) {
			if (m[i][j] <= 0.0) {
				return false;
			}
		}
	}
	return true;
}

bool MatrixUtils::unary_predicate(const Matrix& m, std::function<bool(double)> predicate) {
	if (predicate == nullptr) {
		throw std::invalid_argument("the predicate is null");
	}

	for (int i{0}; i < m.getRows(); ++i) {
		for (int j{0}; j < m.getColumns(); ++j) {
			if (!predicate(m[i][j])) {
				return false;
			}
		}
	}

	return true;
}

bool MatrixUtils::binary_predicate(const Matrix& a, const Matrix& b, std::function<bool(double, double)> predicate) {
	if (a.getRows() != b.getRows()) {
		throw std::invalid_argument("the matrices don't have the same number of rows");
	}
	if (a.getColumns() != b.getColumns()) {
		throw std::invalid_argument("the matrices don't have the same number of columns");
	}
	if (predicate == nullptr) {
		throw std::invalid_argument("the predicate is null");
	}

	for (int i{0}; i < a.getRows(); ++i) {
		for (int j{0}; j < a.getColumns(); ++j) {
			if (!predicate(a[i][j], b[i][j])) {
				return false;
			}
		}
	}

	return true;
}

Matrix MatrixUtils::unary_op(const Matrix& m, std::function<double(double)> op) {
	if (op == nullptr) {
		throw std::invalid_argument("the operation is null");
	}

	Matrix ret(m.getRows(), m.getColumns());

	for (int i{0}; i < m.getRows(); ++i) {
		for (int j{0}; j < m.getColumns(); ++j) {
			ret[i][j] = op(m[i][j]);
		}
	}

	return ret;
}

Matrix MatrixUtils::binary_op(const Matrix& a, const Matrix& b, std::function<double(double, double)> op) {
	if (a.getRows() != b.getRows()) {
		throw std::invalid_argument("the matrices don't have the same number of rows");
	}
	if (a.getColumns() != b.getColumns()) {
		throw std::invalid_argument("the matrices don't have the same number of columns");
	}
	if (op == nullptr) {
		throw std::invalid_argument("the operation is null");
	}

	Matrix ret(a.getRows(), a.getColumns());

	for (int i{0}; i < a.getRows(); ++i) {
		for (int j{0}; j < a.getColumns(); ++j) {
			ret[i][j] = op(a[i][j], b[i][j]);
		}
	}

	return ret;
}
