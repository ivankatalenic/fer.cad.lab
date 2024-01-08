#include "optimization.hh"

#include "float_utils.hh"
#include "vector_utils.hh"
#include "matrix_utils.hh"
#include "matrix_ops.hh"

#include <stdexcept>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include <limits>

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
		throw std::invalid_argument("the provided invalid is invalid (l > r)");
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

static int find_second_highest(Simplex& simplex) {
	int best{0};
	for (int i{1}; i < simplex.size(); ++i) {
		if (simplex.value(i) > simplex.value(best)) {
			best = i;
		}
	}
	int second_best{0};
	for (int i{1}; i < simplex.size(); ++i) {
		if (i == best) {
			continue;
		}
		if (simplex.value(i) > simplex.value(second_best)) {
			second_best = i;
		}
	}
	return second_best;
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

static Matrix centroid(Simplex& simplex) {
	Matrix ret(simplex.get(0));

	for (int i{1}; i < static_cast<int>(simplex.size()); ++i) {
		ret = ret + simplex.get(i);
	}

	return (1.0 / simplex.size()) * ret;
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

Matrix Optimization::gradientDescent(
	std::function<double(const Matrix&)> f,
	std::function<Matrix(const Matrix&)> gradient,
	const Matrix&                        start_v,
	const double                         epsilon,
	const int                            iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (gradient == nullptr) {
		throw std::invalid_argument("the provided gradient is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix pos(start_v);
	for (int i{0}; i < iter_max; ++i) {
#ifdef PRINT
		std::cout << VectorUtils::format(pos) << ": " << f(pos) << std::endl;
#endif

		const auto grad{gradient(pos)};
		if (grad.getLength() < epsilon) {
			return pos;
		}

		pos = pos - grad;
	}
	std::cout << "reached max iterations " << iter_max << " while in gradient descent" << std::endl;
	return pos;
}

Matrix Optimization::gradientDescentGolden(
	std::function<double(const Matrix&)> f,
	std::function<Matrix(const Matrix&)> gradient,
	const Matrix&                        start_v,
	const double                         epsilon,
	const int                            iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (gradient == nullptr) {
		throw std::invalid_argument("the provided gradient is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix pos(start_v);
	for (int i{0}; i < iter_max; ++i) {
#ifdef PRINT
		std::cout << VectorUtils::format(pos) << ": " << f(pos) << std::endl;
#endif

		const auto grad{gradient(pos)};

		const double lambda{goldenSearch(
			[&](double val) -> double {
				return f(pos + val * grad);
			},
			0,
			1.0
		)};
		if (std::abs(lambda * grad.getLength()) < epsilon) {
			return pos;
		}

		pos = pos + lambda * grad;
	}
	std::cout << "reached max iterations " << iter_max << " while in gradient descent" << std::endl;
	return pos;
}

Matrix Optimization::newtonRaphson(
	std::function<double(const Matrix&)> f,
	std::function<Matrix(const Matrix&)> gradient,
	std::function<Matrix(const Matrix&)> hessian,
	const Matrix&                        start_v,
	const double                         epsilon,
	const int                            iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (gradient == nullptr) {
		throw std::invalid_argument("the provided gradient is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix pos(start_v);
	for (int i{0}; i < iter_max; ++i) {
#ifdef PRINT
		std::cout << VectorUtils::format(pos) << ": " << f(pos) << std::endl;
#endif

		const auto grad{gradient(pos)};
		const auto hess{hessian(pos)};
		const auto step{MatrixOps::inverse(hess) * (-1 * grad)};
		if (step.getLength() < epsilon) {
			return pos;
		}

		pos = pos + step;
	}
	std::cout << "reached max iterations " << iter_max << " while in Newton-Raphson method" << std::endl;
	return pos;
}

Matrix Optimization::newtonRaphsonGolden(
	std::function<double(const Matrix&)> f,
	std::function<Matrix(const Matrix&)> gradient,
	std::function<Matrix(const Matrix&)> hessian,
	const Matrix&                        start_v,
	const double                         epsilon,
	const int                            iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (gradient == nullptr) {
		throw std::invalid_argument("the provided gradient is null");
	}
	if (start_v.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix pos(start_v);
	for (int i{0}; i < iter_max; ++i) {
#ifdef PRINT
		std::cout << VectorUtils::format(pos) << ": " << f(pos) << std::endl;
#endif

		const auto grad{gradient(pos)};
		const auto hess{hessian(pos)};
		const auto step{MatrixOps::inverse(hess) * (-1 * grad)};

		const double lambda{goldenSearch(
			[&](double val) -> double {
				return f(pos + val * step);
			},
			0,
			1.0
		)};
		if (std::abs(lambda * step.getLength()) < epsilon) {
			return pos;
		}

		pos = pos + lambda * step;
	}
	std::cout << "reached max iterations " << iter_max << " while in Newton-Raphson method" << std::endl;
	return pos;
}

Matrix Optimization::gaussNewton(
	std::function<Matrix(const Matrix&)> G,
	std::function<Matrix(const Matrix&)> J,
	const Matrix&                        start,
	const double                         epsilon,
	const int                            iter_max
) {
	if (G == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (J == nullptr) {
		throw std::invalid_argument("the provided Jacobian function is null");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix x{start};
	for (int i{0}; i < iter_max; i++) {
		const Matrix Gx{G(x)};
		const Matrix Jx{J(x)};
		const Matrix A{(~Jx) * Jx};
		const Matrix g{(~Jx) * Gx};
		
		const Matrix step{MatrixOps::inverse(A) * (-1 * g)};

		if (step.getLength() < epsilon) {
			return x;
		}

		x = x + step;
	}
	std::cout << "reached max iterations " << iter_max << " while in Gauss-Newton method" << std::endl;
	return x;
}

Matrix Optimization::gaussNewtonGolden(
	std::function<Matrix(const Matrix&)> G,
	std::function<Matrix(const Matrix&)> J,
	const Matrix&                        start,
	const double                         epsilon,
	const int                            iter_max
) {
	if (G == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (J == nullptr) {
		throw std::invalid_argument("the provided Jacobian function is null");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	Matrix x{start};
	for (int i{0}; i < iter_max; i++) {
		const Matrix Gx{G(x)};
		const Matrix Jx{J(x)};
		const Matrix A{(~Jx) * Jx};
		const Matrix g{(~Jx) * Gx};
		
		const Matrix step{MatrixOps::inverse(A) * (-1 * g)};

		const double lambda{goldenSearch(
			[&](double val) -> double {
				const Matrix new_x{x + val * step};
				const Matrix new_Gx{G(new_x)};
				const int rows{new_Gx.getRows()};
				double sum{0};
				for (int i{0}; i < rows; i++) {
					sum += new_Gx[i][0] * new_Gx[i][0];
				}
				return sum;
			},
			0,
			1.0
		)};
		if (std::abs(lambda * step.getLength()) < epsilon) {
			return x;
		}

		x = x + lambda * step;
	}
	std::cout << "reached max iterations " << iter_max << " while in Gauss-Newton method" << std::endl;
	return x;
}

static bool is_inside_bounds(const Matrix& lower, const Matrix& upper, const Matrix& point) {
	if (point.getColumns() != 1 || lower.getColumns() != 1 || upper.getColumns() != 1) {
		throw std::invalid_argument("provided bounds or the point isn't a vector");
	}
	if (point.getRows() != lower.getRows()) {
		throw std::invalid_argument("the dimension of lower bound is not the same as the one from the point");
	}
	if (point.getRows() != upper.getRows()) {
		throw std::invalid_argument("the dimension of upper bound is not the same as the one from the point");
	}

	for (int i{0}; i < point.getRows(); ++i) {
		const double p{point[i][0]};
		const double l{lower[i][0]};
		const double u{upper[i][0]};
		if (p < l || p > u) {
			return false;
		}
	}
	return true;
}

static bool are_constraints_satisfied(
	const std::vector<std::function<double(const Matrix&)>>& constraints,
	const Matrix& point
) {
	for (auto c : constraints) {
		if (c(point) < 0) {
			return false;
		}
	}
	return true;
}

static double random() {
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
	return dist(gen);
}

static Matrix inside_point(
	const std::vector<std::function<double(const Matrix&)>>& constraints,
	const Matrix&                                            start,
	const double                                             epsilon
) {
	Matrix pos(start);
	Matrix old(start);
	do {
		old = pos;
		pos = hookeJeeves(
			[&](const Matrix& v) -> double {
				double ret{0.0};
				for (auto c : constraints) {
					if (c(v) < 0) {
						ret = ret - c(v);
					}
				}
				return ret;
			},
			start,
			VectorUtils::uniform_vector(start.getRows(), 1.0),
			VectorUtils::uniform_vector(start.getRows(), 1e-6)
		);
	} while ((pos - old).getLength() > epsilon);

	return pos;
}

Matrix Optimization::boxMethod(
	std::function<double(const Matrix&)>                     f,
	const Matrix&                                            bound_lower,
	const Matrix&                                            bound_upper,
	const std::vector<std::function<double(const Matrix&)>>& constraints,
	Matrix                                                   start,
	const double                                             epsilon,
	const double                                             alpha,
	const int                                                iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (start.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	const int n{start.getRows()};

	if (!is_inside_bounds(bound_lower, bound_upper, start)
		|| !are_constraints_satisfied(constraints, start)) {
		
		auto constr(constraints);
		for (int i{0}; i < n; ++i) {
			constr.push_back([=](const Matrix& m) -> double {
				return m[i][0] - bound_lower[i][0];
			});
			constr.push_back([=](const Matrix& m) -> double {
				return bound_upper[i][0] - m[i][0];
			});
		}
		start = inside_point(constraints, start, 1e-3);
	}

	Simplex simplex(f);
	Matrix cent(start);

	for (int t{0}; t < 2 * n; ++t) {
		Matrix p(n, 1);
		for (int i{0}; i < n; ++i) {
			const double r{random()};
			p[i][0] = r * bound_lower[i][0] + (1.0 - r) * bound_upper[i][0];
		}
		while (!are_constraints_satisfied(constraints, p)) {
			p = 0.5 * (p + cent);
		}
		simplex.append(p);
		cent = centroid(simplex);
	}

	for (int i{0}; i < iter_max; ++i) {
		const int worst{find_highest(simplex)};
		const int second_worst{find_second_highest(simplex)};

		cent = centroid(filter_worst(simplex));
		Matrix xr{reflect(cent, simplex.get(worst), alpha)};

		for (int j{0}; j < n; ++j) {
			if (xr[j][0] < bound_lower[j][0]) {
				xr[j][0] = bound_lower[j][0];
			} else if (xr[j][0] > bound_upper[j][0]) {
				xr[j][0] = bound_upper[j][0];
			}
		}

		while (!are_constraints_satisfied(constraints, xr)) {
			xr = 0.5 * (xr + cent);
		}

		if (f(xr) > simplex.value(second_worst)) {
			xr = 0.5 * (xr + cent);
		}

		simplex.set(worst, xr);

		if (value_deviation(simplex) < epsilon) {
			break;
		}
	}

	return simplex.get(find_lowest(simplex));
}

static double trans_func(
	std::function<double(const Matrix&)>                     f,
	const std::vector<std::function<double(const Matrix&)>>& constr_ineq,
	const std::vector<std::function<double(const Matrix&)>>& constr_eq,
	Matrix                                                   x,
	double                                                   t
) {
	double sum{f(x)};
	for (auto g : constr_ineq) {
		const double val{g(x)};
		if (val < 0) {
			return std::numeric_limits<double>::max();
		}
		sum += -1.0/t * log(val);
	}
	for (auto h : constr_eq) {
		const double val{h(x)};
		sum += t * val * val;
	}
	return sum;
}

Matrix Optimization::transMethod(
	std::function<double(const Matrix&)>                     f,
	const std::vector<std::function<double(const Matrix&)>>& constr_ineq,
	const std::vector<std::function<double(const Matrix&)>>& constr_eq,
	Matrix                                                   start,
	const double                                             epsilon,
	double                                                   t,
	const int                                                iter_max
) {
	if (f == nullptr) {
		throw std::invalid_argument("the provided function is null");
	}
	if (start.getColumns() != 1) {
		throw std::invalid_argument("the start point matrix doesn't have exactly one column");
	}
	if (epsilon < 0.0) {
		throw std::invalid_argument("the provided epsilon is negative");
	}
	if (iter_max < 0) {
		throw std::invalid_argument("the max number of iterations is negative");
	}

	if (!are_constraints_satisfied(constr_ineq, start)) {
		start = inside_point(constr_ineq, start, 1e-3);
	}

	Matrix pos(start);
	Matrix old(start);
	for (int i{0}; i < iter_max; ++i) {
		old = pos;

		pos = hookeJeeves(
			[&](const Matrix& x) -> double {
				return trans_func(f, constr_ineq, constr_eq, x, t);
			},
			pos,
			VectorUtils::uniform_vector(start.getRows(), 1.0),
			VectorUtils::uniform_vector(start.getRows(), 1e-6)
		);

		if ((pos - old).getLength() < epsilon) {
			return pos;
		}

		t *= 10;
	}

	std::cout << "reached max iterations " << iter_max << " while optimizing with punishment and barrier method" << std::endl;
	return pos;
}
