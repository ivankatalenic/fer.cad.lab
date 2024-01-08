#pragma once

#include "matrix.hh"

#include <string>
#include <istream>
#include <ostream>

namespace MatrixUtils {
	void print(const Matrix& m, std::ostream& o);
	void printToConsole(const Matrix& m);
	void printToFile(const Matrix& m, const std::string& filename);

	Matrix load(std::istream& i);
	Matrix loadFromFile(const std::string& filename);
}
