#include "optimization.hh"

#include "float_utils.hh"
#include "vector_utils.hh"
#include "matrix_utils.hh"

#include <stdexcept>
#include <cmath>
#include <vector>

#ifdef PRINT
#include <iostream>
#endif

using namespace Optimization;

Interval Optimization::unimodalIntervalSearch(
	std::function<double(double)> f,
	double                        start,
	double                        step
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is nul");
	}
	if (step < 0.0 || FloatUtils::is_equal(step, 0.0, 1e-6)) {
		throw std::invalid_argument("the step must be positive");
	}

	Interval interval{start - step, start + step};

	double m{start};
	double fl, fm, fr;
	int    step_factor{1};

	fm = f(start);
	fl = f(interval.l);
	fr = f(interval.r);

	if (fm < fr && fm < fl) {
		return interval;
	} else if (fm > fr) {
		do {
			interval.l = m;
			m          = interval.r;
			fm         = fr;
			interval.r = start + step * (step_factor *= 2);
			fr         = f(interval.r);
		} while (fm > fr);
	} else {
		do {
			interval.r = m;
			m          = interval.l;
			fm         = fl;
			interval.l = start - step * (step_factor *= 2);
			fl         = f(interval.l);
		} while (fm > fl);
	}

	return interval;
}

Interval Optimization::goldenSectionSearch(
	std::function<double(double)> f,
	Interval                      interval,
	double                        precision
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is nil");
	}
	if (interval.l > interval.r) {
		throw std::invalid_argument("the provided interval is invalid (l > r)");
	}
	if (precision < 0) {
		throw std::invalid_argument("the precision is negative");
	}

	const double K{0.5 * (std::sqrt(5.0) - 1.0)};

	double c{interval.r - K * (interval.r - interval.l)};
	double d{interval.l + K * (interval.r - interval.l)};
	double fc{f(c)};
	double fd{f(d)};

	while ((interval.r - interval.l) > precision) {
#ifdef PRINT
		const double fl{f(interval.l)};
		const double fr{f(interval.r)};

		std::cout << "a="   << interval.l;
		std::cout << "\tc=" << c;
		std::cout << "\td=" << d;
		std::cout << "\tb=" << interval.r << std::endl;

		std::cout << "f(a)="   << fl;
		std::cout << "\tf(c)=" << fc;
		std::cout << "\tf(d)=" << fd;
		std::cout << "\tf(b)=" << fr << std::endl;

		std::cout << "--------------------------------" << std::endl;
#endif
		if (fc < fd) {
			interval.r = d;
			d          = c;
			c          = interval.r - K * (interval.r - interval.l);
			fd         = fc;
			fc         = f(c);
		} else {
			interval.l = c;
			c          = d;
			d          = interval.l + K * (interval.r - interval.l);
			fc         = fd;
			fd         = f(d);
		}
	}

	return interval;
}

double Optimization::goldenSearch(
	std::function<double(double)> f,
	double                        start,
	double                        step,
	double                        precision
) {
	Interval unimodal{Optimization::unimodalIntervalSearch(f, start, step)};
	Interval res{Optimization::goldenSectionSearch(f, unimodal, precision)};

	return (res.l + res.r) / 2.0;
}

Matrix Optimization::axisSearch(
	std::function<double(const Matrix&)> f,
	const Matrix&                        start_v,
	double                               step,
	const Matrix&                        precision_v
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (precision_v.getColumns() != 1) {
		throw std::invalid_argument("the precision matrix vector doesn't have exactly one column");
	}
	if (precision_v.getRows() != start_v.getRows()) {
		throw std::invalid_argument("the precision vector doesn't have the size equal to the start vector");
	}

	const int n{start_v.getRows()};
	Matrix    res_v(start_v);
	Matrix    base_v(res_v);
	do {
		base_v = res_v;
		for (int i{0}; i < n; ++i) {
			const Matrix direction_v(VectorUtils::identity_vector(i, n));

			std::function<double(double)> simple_func([&](double l) -> double {
				return f(res_v + l * direction_v);
			});
			const double lambda{Optimization::goldenSearch(simple_func, 0.0, step, precision_v[i][0])};

			res_v = res_v + lambda * direction_v;
		}
	} while (MatrixUtils::is_less(precision_v, MatrixUtils::abs(res_v - base_v)));

	return res_v;
}

static Matrix explore(
	std::function<double(const Matrix&)> f,
	const Matrix&                        start,
	const Matrix&                        delta
) {
	Matrix x(start);
	for (int i{0}; i < start.getRows(); ++i) {
		const double p{f(x)};
		x[i][0] = x[i][0] + delta[i][0];
		double n{f(x)};

		if (n > p) {
			x[i][0] = x[i][0] - 2 * delta[i][0];
			n = f(x);
			if (n > p) {
				x[i][0] = x[i][0] + delta[i][0];
			}
		}
	}
	return x;
}

Matrix Optimization::hookeJeeves(
	std::function<double(const Matrix&)> f,
	const Matrix&                        start_v,
	Matrix                               delta_v,
	const Matrix&                        precision_v
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (precision_v.getColumns() != 1) {
		throw std::invalid_argument("the precision matrix vector doesn't have exactly one column");
	}
	if (precision_v.getRows() != start_v.getRows()) {
		throw std::invalid_argument("the precision vector doesn't have the size equal to the start vector");
	}
	if (delta_v.getColumns() != 1) {
		throw std::invalid_argument("the delta matrix vector doesn't have exactly one column");
	}
	if (delta_v.getRows() != start_v.getRows()) {
		throw std::invalid_argument("the delta vector doesn't have the size equal to the start vector");
	}
	if (!MatrixUtils::is_positive(delta_v)) {
		throw std::invalid_argument("the delta vector has negative values");
	}
	if (!MatrixUtils::is_positive(precision_v)) {
		throw std::invalid_argument("the precision vector has negative values");
	}

	Matrix xp(start_v);
	Matrix xb(start_v);
	Matrix xn(start_v);

	do {
		xn = explore(f, xp, delta_v);
		const double fxn{f(xn)};
		const double fxb{f(xb)};
#ifdef PRINT
		const double fxp{f(xp)};

		std::cout << "Xb="   << VectorUtils::format(xb);
		std::cout << "\tXp=" << VectorUtils::format(xp);
		std::cout << "\tXn=" << VectorUtils::format(xn) << std::endl;

		std::cout << "f(Xb)="   << fxb;
		std::cout << "\tf(Xp)=" << fxp;
		std::cout << "\tf(Xn)=" << fxn << std::endl;

		std::cout << "--------------------------------" << std::endl;
#endif

		if (fxn < fxb) {
			xp = 2 * xn - xb;
			xb = xn;
		} else {
			delta_v = delta_v * 0.5;
			xp = xb;
		}
	} while (MatrixUtils::is_less(precision_v, delta_v));

	return xb;
}

struct PointValue {
	Matrix point;
	double value;
	bool   value_valid;
};

class Simplex {
public:
	Simplex(std::function<double(const Matrix&)> f): func{f} {};
	void append(Matrix m) {
		point_values.push_back({m, 0.0, false});
	}
	Matrix get(int i) {
		return point_values[i].point;
	}
	void set(int i, Matrix m) {
		point_values[i] = {m, 0.0, false};
	}
	void set(int i, Matrix m, double value) {
		point_values[i] = {m, value, true};
	}
	double value(int i) {
		PointValue& pv{point_values[i]};
		if (!pv.value_valid) {
			pv.value       = func(pv.point);
			pv.value_valid = true;
		}
		return pv.value;
	}
	int size() {
		return point_values.size();
	}
private:
	std::function<double(const Matrix&)> func;
	std::vector<PointValue>              point_values; 
};

static int find_highest(Simplex& simplex) {
	int best{0};
	for (int i{1}; i < simplex.size(); ++i) {
		if (simplex.value(i) > simplex.value(best)) {
			best = i;
		}
	}
	return best;
}

static int find_lowest(Simplex& simplex) {
	int worst{0};
	for (int i{1}; i < simplex.size(); ++i) {
		if (simplex.value(i) < simplex.value(worst)) {
			worst = i;
		}
	}
	return worst;
}

static std::vector<Matrix> filter_worst(Simplex& simplex) {
	std::vector<Matrix> ret;
	const int worst{find_highest(simplex)};
	for (int i{0}; i < simplex.size(); ++i) {
		if (i == worst) {
			continue;
		}

		ret.push_back(simplex.get(i));
	}

	return ret;
}

static Matrix centroid(const std::vector<Matrix>& list) {
	Matrix ret(list[0]);

	for (int i{1}; i < static_cast<int>(list.size()); ++i) {
		ret = ret + list[i];
	}

	return (1.0 / list.size()) * ret;
}

static Matrix reflect(const Matrix& centroid, const Matrix& high, double alpha) {
	return (1.0 + alpha) * centroid - alpha * high;
}

static Matrix expand(const Matrix& centroid, const Matrix& reflection, double gamma) {
	return (1.0 - gamma) * centroid + gamma * reflection;
}

static Matrix contract(const Matrix& centroid, const Matrix& high, double beta) {
	return (1.0 - beta) * centroid + beta * high;
}

static void shrink(Simplex& simplex, const Matrix& point, double sigma) {
	for (int i{0}; i < simplex.size(); ++i) {
		Matrix new_point(point + sigma * (simplex.get(i) - point));
		simplex.set(i, new_point);
	}
}

static double value_deviation(Simplex& simplex) {
	const int n{simplex.size()};
	double    avg_sum{0.0};
	for (int i{0}; i < n; ++i) {
		avg_sum += simplex.value(i);
	}
	const double avg{avg_sum / n};

	double dev{0.0};
	for (int i{0}; i < n; ++i) {
		const double d{simplex.value(i) - avg};
		dev += d * d;
	}

	return std::sqrt(dev / (n + 1));
}

static double pos_deviation(Simplex& simplex) {
	const int n{simplex.size()};
	Matrix centroid(simplex.get(0));
	for (int i{1}; i < n; ++i) {
		centroid = centroid + simplex.get(i);
	}
	centroid = (1.0 / n) * centroid;

	double dev{0.0};
	for (int i{0}; i < n; ++i) {
		const double d{VectorUtils::length(simplex.get(i) - centroid)};
		dev += d * d;
	}

	return std::sqrt(dev / (n + 1));
}

Matrix Optimization::nelderMead(
	std::function<double(const Matrix&)> f,
	const Matrix&                        start_v,
	const Matrix&                        step_v,
	double                               epsilon,
	double                               alpha,
	double                               beta,
	double                               gamma,
	double                               sigma
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (step_v.getColumns() != 1) {
		throw std::invalid_argument("the step matrix vector doesn't have exactly one column");
	}
	if (step_v.getRows() != start_v.getRows()) {
		throw std::invalid_argument("the step vector doesn't have the size equal to the start vector");
	}
	if (!MatrixUtils::is_positive(step_v)) {
		throw std::invalid_argument("the step vector has negative values");
	}

	const int n{start_v.getRows()};
	Simplex simplex(f);
	simplex.append(start_v);
	for (int i{0}; i < n; ++i) {
		simplex.append(step_v[i][0] * VectorUtils::identity_vector(i, n) + start_v);
	}

	Matrix xc(centroid(filter_worst(simplex)));
	do {
		const int h{find_highest(simplex)};
		const int l{find_lowest(simplex)};
		xc = centroid(filter_worst(simplex));
#ifdef PRINT
		std::cout << "Xc=" << VectorUtils::format(xc);
		std::cout << "\tf(Xc)=" << f(xc) << std::endl;
#endif
		
		const Matrix xr(reflect(xc, simplex.get(h), alpha));
		const double fxr{f(xr)};

		if (fxr < simplex.value(l)) {
			const Matrix xe(expand(xc, xr, gamma));
			const double fxe{f(xe)};

			if (fxe < simplex.value(l)) {
				simplex.set(h, xe, fxe);
			} else {
				simplex.set(h, xr, fxr);
			}
		} else {
			// Check if the reflected point is worse than all simplex points except the worst simplex point
			bool is_worse{true};
			for (int j{0}; j <= n; ++j) {
				if (j == h) {
					continue;
				}

				if (fxr <= simplex.value(j)) {
					is_worse = false;
					break;
				}
			}

			if (is_worse) {
				if (fxr < simplex.value(h)) {
					simplex.set(h, xr, fxr);
				}

				const Matrix xk(contract(xc, simplex.get(h), beta));
				const double fxk{f(xk)};

				if (fxk < simplex.value(h)) {
					simplex.set(h, xk, fxk);
				} else {
					shrink(simplex, simplex.get(l), sigma);
				}
			} else {
				simplex.set(h, xr, fxr);
			}
		}
	} while (value_deviation(simplex) > epsilon || pos_deviation(simplex) > epsilon);

	return xc;
}
