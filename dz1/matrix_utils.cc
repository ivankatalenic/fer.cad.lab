#include "matrix_utils.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

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
