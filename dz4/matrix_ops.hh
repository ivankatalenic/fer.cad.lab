#pragma once

#include "matrix.hh"

#include <tuple>

namespace MatrixOps {
	Matrix subForward(const Matrix& L, const Matrix& b);
	Matrix subBackward(const Matrix& U, const Matrix& y);

	Matrix                     decomposeLU(const Matrix& m);
	std::tuple<Matrix, Matrix> decomposeLUP(const Matrix& m);

	Matrix inverse(const Matrix& m);

	double determinant(const Matrix& m);
}
