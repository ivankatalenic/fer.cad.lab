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
}
