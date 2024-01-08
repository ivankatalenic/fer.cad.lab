#pragma once

#include "matrix.hh"

#include <cmath>
#include <stdexcept>

namespace fun {
	namespace first {
		double f(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return 100 * pow(x2 - pow(x1, 2), 2) + pow(1 - x1, 2);
		}
		Matrix gradient(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return Matrix{
				{2 * (200 * pow(x1, 3) - 200 * x1 * x2 + x1 - 1)},
				{200 * (x2 - pow(x1, 2))}
			};
		}
		Matrix hessian(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};

			return Matrix{
				{-400 * (x2 - pow(x1, 2)) + 800 * pow(x1, 2) + 2, -400 * x1},
				{-400 * x1, 200}
			};
		}
	}
	namespace second {
		double f(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return pow(x1 - 4, 2) + 4 * pow(x2 - 2, 2);
		}
		Matrix gradient(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return Matrix{
				{2 * (x1 - 4)},
				{8 * (x2 - 2)}
			};
		}
		Matrix hessian(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}

			return Matrix{
				{2, 0},
				{0, 8}
			};
		}
	}
	namespace third {
		double f(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return pow(x1 - 2, 2) + pow(x2 + 3, 2);
		}
		Matrix gradient(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x2{v[1][0]};
			return Matrix{
				{2 * (x1 - 2)},
				{2 * (x2 + 3)}
			};
		}
		Matrix hessian(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}

			return Matrix{
				{2, 0},
				{0, 2}
			};
		}
	}
	namespace fourth {
		double f(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x1_2{pow(x1, 2)};
			const double x1_4{pow(x1, 4)};
			const double x2{v[1][0]};
			return 0.25 * x1_4 - x1_2 + 2 * x1 + pow(x2 - 1, 2);
		}
		Matrix gradient(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x1_3{pow(x1, 3)};
			const double x2{v[1][0]};
			return Matrix{
				{x1_3 - 2 * x1 + 2},
				{2 * (x2 - 1)     }
			};
		}
		Matrix hessian(const Matrix& v) {
			if (v.getColumns() != 1) {
				throw std::invalid_argument("the input is not a vector");
			}
			if (v.getRows() != 2) {
				throw std::invalid_argument("the input is not a 2-vector");
			}
			const double x1{v[0][0]};
			const double x1_2{pow(x1, 2)};
			return Matrix{
				{3 * x1_2 - 2, 0},
				{0           , 2}
			};
		}
	}
}
