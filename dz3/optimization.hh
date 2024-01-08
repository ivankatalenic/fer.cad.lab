#pragma once

#include "matrix.hh"

#include <functional>

namespace Optimization {

	struct Interval {
		double l;
		double r;
	};

	Interval unimodalIntervalSearch(
		std::function<double(double)> f,
		double                        start,
		double                        step
	);
	Interval goldenSectionSearch(
		std::function<double(double)> f,
		Interval                      interval,
		double                        precision = 1e-6
	);
	double goldenSearch(
		std::function<double(double)> f,
		double                        start,
		double                        step,
		double                        precision = 1e-6
	);

	Matrix axisSearch(
		std::function<double(const Matrix&)> f,
		const Matrix&                        start_v,
		double                               step,
		const Matrix&                        precision_v
	);
	Matrix hookeJeeves(
		std::function<double(const Matrix&)> f,
		const Matrix&                        start_v,
		Matrix                               delta_v,
		const Matrix&                        precision_v
	);
	Matrix nelderMead(
		std::function<double(const Matrix&)> f,
		const Matrix&                        start_v, 
		const Matrix&                        step_v,
		double                               epsilon = 1e-6,
		double                               alpha   = 1.0,
		double                               beta    = 0.5,
		double                               gamma   = 2.0,
		double                               sigma   = 0.5
	);

	Matrix gradientDescent(
		std::function<double(const Matrix&)> f,
		std::function<Matrix(const Matrix&)> gradient,
		const Matrix&                        start_v,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix gradientDescentGolden(
		std::function<double(const Matrix&)> f,
		std::function<Matrix(const Matrix&)> gradient,
		const Matrix&                        start_v,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix newtonRaphson(
		std::function<double(const Matrix&)> f,
		std::function<Matrix(const Matrix&)> gradient,
		std::function<Matrix(const Matrix&)> hessian,
		const Matrix&                        start_v,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix newtonRaphsonGolden(
		std::function<double(const Matrix&)> f,
		std::function<Matrix(const Matrix&)> gradient,
		std::function<Matrix(const Matrix&)> hessian,
		const Matrix&                        start_v,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix gaussNewton(
		std::function<Matrix(const Matrix&)> G,
		std::function<Matrix(const Matrix&)> J,
		const Matrix&                        start,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix gaussNewtonGolden(
		std::function<Matrix(const Matrix&)> G,
		std::function<Matrix(const Matrix&)> J,
		const Matrix&                        start,
		const double                         epsilon,
		const int                            iter_max
	);
	Matrix boxMethod(
		std::function<double(const Matrix&)>                     f,
		const Matrix&                                            bound_lower,
		const Matrix&                                            bound_upper,
		const std::vector<std::function<double(const Matrix&)>>& constraints,
		Matrix                                                   start,
		const double                                             epsilon,
		const double                                             alpha,
		const int                                                iter_max
	);
	Matrix transMethod(
		std::function<double(const Matrix&)>                     f,
		const std::vector<std::function<double(const Matrix&)>>& constr_ineq,
		const std::vector<std::function<double(const Matrix&)>>& constr_eq,
		Matrix                                                   start,
		const double                                             epsilon,
		double                                                   t,
		const int                                                iter_max
	);
}
