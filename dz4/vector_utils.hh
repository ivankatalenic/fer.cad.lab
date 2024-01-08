#pragma once

#include "matrix.hh"

#include <string>

namespace VectorUtils {
	Matrix identity_vector(int position, int size);
	Matrix uniform_vector(int size, double val);
	double length(const Matrix& m);
	
	std::string format(const Matrix& v);
}
