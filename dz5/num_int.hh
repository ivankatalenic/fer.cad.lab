#pragma once

#include "matrix.hh"

#include <functional>

namespace Predictor {
	Matrix euler(
		const Matrix& x,
		const Matrix& A,
		const Matrix& B,
		double t,
		double period,
		bool is_time_dependent
	);
}

namespace Corrector {
	Matrix euler_inverted(
		const Matrix& x,
		const Matrix& x_predicted,
		const Matrix& A,
		const Matrix& B,
		double t,
		double period,
		bool is_time_dependent
	);
	Matrix trapeze(
		const Matrix& x,
		const Matrix& x_predicted,
		const Matrix& A,
		const Matrix& B,
		double t,
		double period,
		bool is_time_dependent
	);
}

namespace NumInt {
	Matrix euler(
		const Matrix& A,
		const Matrix& B,
		const Matrix& initial,
		const double period,
		const double t_max,
		const bool is_time_dependent,
		std::function<void(double t, const Matrix& x)> callback = nullptr
	);
	Matrix euler_inverted(
		const Matrix& A,
		const Matrix& B,
		const Matrix& initial,
		const double period,
		const double t_max,
		const bool is_time_dependent,
		std::function<void(double t, const Matrix& x)> callback = nullptr
	);
	Matrix trapeze(
		const Matrix& A,
		const Matrix& B,
		const Matrix& initial,
		const double period,
		const double t_max,
		const bool is_time_dependent,
		std::function<void(double t, const Matrix& x)> callback = nullptr
	);
	Matrix runge_kutta(
		const Matrix& A,
		const Matrix& B,
		const Matrix& initial,
		const double period,
		const double t_max,
		const bool is_time_dependent,
		std::function<void(double t, const Matrix& x)> callback = nullptr
	);
	Matrix predictor_corrector(
		std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, double, double, bool)> predictor,
		std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, const Matrix&, double, double, bool)> corrector,
		const Matrix& A,
		const Matrix& B,
		const Matrix& initial,
		const double period,
		const double t_max,
		const bool is_time_dependent,
		const int iter_max,
		std::function<void(double t, const Matrix& x)> callback = nullptr
	);
}
