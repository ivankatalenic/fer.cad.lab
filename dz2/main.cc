#include "matrix.hh"
#include "matrix_utils.hh"
#include "matrix_ops.hh"
#include "optimization.hh"
#include "vector_utils.hh"
#include "float_utils.hh"

#include <iostream>
#include <functional>
#include <cmath>
#include <random>

void zad1() {
	std::cout << __func__ << std::endl;

	int evals;

	std::function<double(double)> f([&](double x) -> double {
		evals++;

		const double y{x - 3};
		return y * y;
	});
	std::function<double(const Matrix&)> fm([&](const Matrix& m) -> double {
		evals++;

		const double y{m[0][0] - 3};
		return y * y;
	});

	const double start{10.0};
	const Matrix startm{{10.0}};

	double min;
	Matrix minm(1, 1);

	double step{1.0};
	Matrix stepm{{1.0}};
	
	double precision{1e-6};
	Matrix precisionm{{1e-6}};

	evals = 0;
	min = Optimization::goldenSearch(f, start, step, precision);
	std::cout << "golden section search:" << std::endl;
	std::cout << "min: "   << min << std::endl;
	std::cout << "evals: " << evals << std::endl << std::endl;

	evals = 0;
	minm = Optimization::axisSearch(fm, startm, step, precisionm);
	std::cout << "axis search:" << std::endl;
	std::cout << "min: "   << VectorUtils::format(minm) << std::endl;
	std::cout << "evals: " << evals << std::endl << std::endl;

	evals = 0;
	minm = Optimization::nelderMead(fm, startm, stepm, precision);
	std::cout << "nelder mead search:" << std::endl;
	std::cout << "min: "   << VectorUtils::format(minm) << std::endl;
	std::cout << "evals: " << evals << std::endl << std::endl;

	evals = 0;
	minm = Optimization::hookeJeeves(fm, startm, stepm, precisionm);
	std::cout << "hooke jeeves search:" << std::endl;
	std::cout << "min: "   << VectorUtils::format(minm) << std::endl;
	std::cout << "evals: " << evals << std::endl << std::endl;
	
	std::cout << std::endl;
}

void zad2() {
	std::cout << __func__ << std::endl;

	int eval1, eval2, eval3, eval4;
	
	std::function<double(const Matrix&)> f1([&](const Matrix& m) -> double {
		eval1++;

		const double a{m[1][0] - m[0][0] * m[0][0]};
		const double b{1 - m[0][0]};
		return 100 * a * a + b * b;
	});
	std::function<double(const Matrix&)> f2([&](const Matrix& m) -> double {
		eval2++;

		const double a{m[0][0] - 4};
		const double b{m[1][0] - 2};
		return a * a + 4 * b * b;
	});
	std::function<double(const Matrix&)> f3([&](const Matrix& m) -> double {
		eval3++;

		double sum{0.0};
		for (int i{1}; i <= 5; ++i) {
			const double y{m[i - 1][0] - i};
			sum += y * y;
		}
		return sum;
	});
	std::function<double(const Matrix&)> f4([&](const Matrix& m) -> double {
		eval4++;

		const double a{(m[0][0] - m[1][0]) * (m[0][0] + m[1][0])};
		const double b{m[0][0] * m[0][0] + m[1][0] * m[1][0]};
		return std::abs(a) + std::sqrt(b);
	});

	const Matrix start1{{-1.9}, {1.0}}, start2{{0.1}, {0.3}}, start3{{0}, {0}, {0}, {0}, {0}}, start4{{5.1}, {1.1}};

	Matrix min1(2, 1), min2(2, 1), min3(5, 1), min4(2,1);

	double step{1.0};
	Matrix step1{{1.0}, {1.0}}, step2{{1.0}, {1.0}}, step3{{1.0}, {1.0}, {1.0}, {1.0}, {1.0}}, step4{{1.0}, {1.0}};
	Matrix delta1{{1.0}, {1.0}}, delta2{{1.0}, {1.0}}, delta3{{1.0}, {1.0}, {1.0}, {1.0}, {1.0}}, delta4{{1.0}, {1.0}};
	
	const Matrix prec1{{1e-6}, {1e-6}}, prec2{{1e-6}, {1e-6}}, prec3{{1e-6}, {1e-6}, {1e-6}, {1e-6}, {1e-6}}, prec4{{1e-6}, {1e-6}};
	const double epsilon{1e-6};

	std::cout << "method\tf1\tf2\tf3\tf4" << std::endl;
	
	std::cout << "------------------------------------" << std::endl;
	eval1 = eval2 = eval3 = eval4 = 0;
	min1 = Optimization::axisSearch(f1, start1, step, prec1);
	min2 = Optimization::axisSearch(f2, start2, step, prec2);
	min3 = Optimization::axisSearch(f3, start3, step, prec3);
	min4 = Optimization::axisSearch(f4, start4, step, prec4);
	std::cout << "axis:\t" 
		<< VectorUtils::format(min1) << "\t"
		<< VectorUtils::format(min2) << "\t"
		<< VectorUtils::format(min3) << "\t"
		<< VectorUtils::format(min4) << std::endl;
	std::cout << "evals:\t"
		<< eval1 << "\t"
		<< eval2 << "\t"
		<< eval3 << "\t"
		<< eval3 << std::endl;

	std::cout << "------------------------------------" << std::endl;
	eval1 = eval2 = eval3 = eval4 = 0;
	min1 = Optimization::nelderMead(f1, start1, step1, epsilon);
	min2 = Optimization::nelderMead(f2, start2, step2, epsilon);
	min3 = Optimization::nelderMead(f3, start3, step3, epsilon);
	min4 = Optimization::nelderMead(f4, start4, step4, epsilon);
	std::cout << "nelder-mead:\t"
		<< VectorUtils::format(min1) << "\t"
		<< VectorUtils::format(min2) << "\t"
		<< VectorUtils::format(min3) << "\t"
		<< VectorUtils::format(min4) << std::endl;
	std::cout << "evals:\t"
		<< eval1 << "\t"
		<< eval2 << "\t"
		<< eval3 << "\t"
		<< eval3 << std::endl;

	std::cout << "------------------------------------" << std::endl;
	eval1 = eval2 = eval3 = eval4 = 0;
	min1 = Optimization::hookeJeeves(f1, start1, delta1, prec1);
	min2 = Optimization::hookeJeeves(f2, start2, delta2, prec2);
	min3 = Optimization::hookeJeeves(f3, start3, delta3, prec3);
	min4 = Optimization::hookeJeeves(f4, start4, delta4, prec4);
	std::cout << "hooke-jeeves:\t" 
		<< VectorUtils::format(min1) << "\t"
		<< VectorUtils::format(min2) << "\t"
		<< VectorUtils::format(min3) << "\t"
		<< VectorUtils::format(min4) << std::endl;
	std::cout << "evals:\t"
		<< eval1 << "\t"
		<< eval2 << "\t"
		<< eval3 << "\t"
		<< eval3 << std::endl;
	
	std::cout << std::endl;
}

void zad3() {
	int eval;

	std::function<double(const Matrix&)> f4([&](const Matrix& m) -> double {
		eval++;

		const double a{(m[0][0] - m[1][0]) * (m[0][0] + m[1][0])};
		const double b{m[0][0] * m[0][0] + m[1][0] * m[1][0]};
		return std::abs(a) + std::sqrt(b);
	});

	const Matrix start{{5}, {5}};
	Matrix min(2, 1);
	Matrix step{{1.0}, {1.0}};
	Matrix delta{{1.0}, {1.0}};
	const Matrix prec{{1e-6}, {1e-6}};
	const double epsilon{1e-6};

	std::cout << "methods\tf4" << std::endl;

	std::cout << "------------------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f4, start, step, epsilon);
	std::cout << "nelder-mead: "<< VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl;

	std::cout << "------------------------------------" << std::endl;
	eval =  0;
	min = Optimization::hookeJeeves(f4, start, delta, prec);
	std::cout << "hooke-jeeves: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl;
}

void zad4() {
	int eval;

	std::function<double(const Matrix&)> f1([&](const Matrix& m) -> double {
		eval++;

		const double a{m[1][0] - m[0][0] * m[0][0]};
		const double b{1 - m[0][0]};
		return 100 * a * a + b * b;
	});

	Matrix start{{0.5}, {0.5}};
	Matrix step{{1.0}, {1.0}};
	const double epsilon{1e-6};
	Matrix min(2, 1);
	
	step = {{1.0}, {1.0}};
	std::cout << "step = 1, start = (0.5, 0.5) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{2.0}, {2.0}};
	std::cout << "step = 2, start = (0.5, 0.5) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{4.0}, {4.0}};
	std::cout << "step = 4, start = (0.5, 0.5) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{8.0}, {8.0}};
	std::cout << "step = 8, start = (0.5, 0.5) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{16.0}, {16.0}};
	std::cout << "step = 16, start = (0.5, 0.5) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	start = {{20}, {20}};

	step = {{1.0}, {1.0}};
	std::cout << "step = 1, start = (20, 20) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{2.0}, {2.0}};
	std::cout << "step = 2, start = (20, 20) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{4.0}, {4.0}};
	std::cout << "step = 4, start = (20, 20) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{8.0}, {8.0}};
	std::cout << "step = 8, start = (20, 20) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;

	step = {{16.0}, {16.0}};
	std::cout << "step = 16, start = (20, 20) --------------------------" << std::endl;
	eval = 0;
	min = Optimization::nelderMead(f1, start, step, epsilon);
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl << std::endl;
}

void zad5() {
	int eval{0};
	std::function<double(const Matrix&)> f6([&](const Matrix& m) -> double {
		eval++;

		const double square_sum{m[0][0] * m[0][0] + m[1][0] * m[1][0]};
		const double sin{std::sin(std::sqrt(square_sum))};
		const double square_sin{sin * sin};
		const double up{square_sin - 0.5};
		const double big{1 + 0.001 * square_sum};
		const double down{big * big};

		return 0.5 + up / down;
	});

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-50, 50);

	const Matrix start{{dis(gen)}, {dis(gen)}};
	const Matrix step{{1.0}, {1.0}};
	const Matrix delta{{1.0}, {1.0}};
	const double epsilon{1e-6};

	std::cout << "start: " << VectorUtils::format(start) << std::endl;
	Matrix min(Optimization::nelderMead(f6, start, step, epsilon));
	std::cout << "x_min: " << VectorUtils::format(min) << std::endl;
	std::cout << "eval: " << eval << std::endl;
	std::cout << "is min found: " << (f6(min) < 1e-4) << std::endl;
}

int main() {
	std::cout << std::boolalpha;

	// zad1();
	// zad2();
	// zad3();
	// zad4();
	zad5();

	return 0;
}
