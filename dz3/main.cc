#include "matrix.hh"
#include "matrix_utils.hh"
#include "matrix_ops.hh"
#include "optimization.hh"
#include "vector_utils.hh"
#include "float_utils.hh"
#include "functions.hh"

#include <iostream>
#include <functional>
#include <cmath>
#include <random>

void zad1() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f;
	int evals_grad;

	std::function<double(const Matrix&)> f([&](const Matrix& x) -> double {
		evals_f++;
		return fun::third::f(x);
	});
	std::function<Matrix(const Matrix&)> gradient([&](const Matrix& x) -> Matrix {
		evals_grad++;
		return fun::third::gradient(x);
	});

	const int iter_max{10000000};
	const Matrix start{{0}, {0}};
	Matrix min(1, 1);

	evals_f = 0;
	evals_grad = 0;
	min = Optimization::gradientDescent(f, gradient, start, 1e-6, iter_max);
	std::cout << "regular gradient descent:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_grad << std::endl << std::endl;

	evals_f = 0;
	evals_grad = 0;
	min = Optimization::gradientDescentGolden(f, gradient, start, 1e-6, iter_max);
	std::cout << "golden gradient descent:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_grad << std::endl << std::endl;
	
	std::cout << std::endl;
}

void zad2() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f1;
	int evals_f1_grad;
	int evals_f1_hess;
	int evals_f2;
	int evals_f2_grad;
	int evals_f2_hess;

	std::function<double(const Matrix&)> f1([&](const Matrix& x) -> double {
		evals_f1++;
		return fun::first::f(x);
	});
	std::function<Matrix(const Matrix&)> f1_grad([&](const Matrix& x) -> Matrix {
		evals_f1_grad++;
		return fun::first::gradient(x);
	});
	std::function<Matrix(const Matrix&)> f1_hess([&](const Matrix& x) -> Matrix {
		evals_f1_hess++;
		return fun::first::hessian(x);
	});

	std::function<double(const Matrix&)> f2([&](const Matrix& x) -> double {
		evals_f2++;
		return fun::second::f(x);
	});
	std::function<Matrix(const Matrix&)> f2_grad([&](const Matrix& x) -> Matrix {
		evals_f2_grad++;
		return fun::second::gradient(x);
	});
	std::function<Matrix(const Matrix&)> f2_hess([&](const Matrix& x) -> Matrix {
		evals_f2_hess++;
		return fun::second::hessian(x);
	});

	const int    iter_max{100000};
	const Matrix start1{{-1.9}, {2}};
	const Matrix start2{{0.1}, {0.3}};
	Matrix       min(2, 1);
	const double epsilon{1e-6};

	std::cout << ">>>>>>>>>>>>>> f1" << std::endl;

	evals_f1      = 0;
	evals_f1_grad = 0;
	evals_f1_hess = 0;
	min = Optimization::gradientDescentGolden(f1, f1_grad, start1, epsilon, iter_max);
	std::cout << "golden gradient descent:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start1) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f1 << std::endl;
	std::cout << "gradient evals: " << evals_f1_grad << std::endl;
	std::cout << "hessian evals: " << evals_f1_hess << std::endl << std::endl;

	evals_f1      = 0;
	evals_f1_grad = 0;
	evals_f1_hess = 0;
	min = Optimization::newtonRaphsonGolden (f1, f1_grad, f1_hess, start1, epsilon, iter_max);
	std::cout << "Newton-Raphson golden search method:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start1) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f1 << std::endl;
	std::cout << "gradient evals: " << evals_f1_grad << std::endl;
	std::cout << "hessian evals: " << evals_f1_hess << std::endl << std::endl;

	std::cout << ">>>>>>>>>>>>> f2" << std::endl;

	evals_f2      = 0;
	evals_f2_grad = 0;
	evals_f2_hess = 0;
	min = Optimization::gradientDescentGolden(f2, f2_grad, start2, epsilon, iter_max);
	std::cout << "golden gradient descent:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start2) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f2 << std::endl;
	std::cout << "gradient evals: " << evals_f2_grad << std::endl;
	std::cout << "hessian evals: " << evals_f2_hess << std::endl << std::endl;

	evals_f2      = 0;
	evals_f2_grad = 0;
	evals_f2_hess = 0;
	min = Optimization::newtonRaphsonGolden(f2, f2_grad, f2_hess, start2, epsilon, iter_max);
	std::cout << "Newton-Raphson golden search method:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start2) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f2 << std::endl;
	std::cout << "gradient evals: " << evals_f2_grad << std::endl;
	std::cout << "hessian evals: " << evals_f2_hess << std::endl << std::endl;
	
	std::cout << std::endl;
}

void zad3() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f;
	int evals_f_grad;
	int evals_f_hess;

	std::function<double(const Matrix&)> f([&](const Matrix& x) -> double {
		evals_f++;
		return fun::fourth::f(x);
	});
	std::function<Matrix(const Matrix&)> f_grad([&](const Matrix& x) -> Matrix {
		evals_f_grad++;
		return fun::fourth::gradient(x);
	});
	std::function<Matrix(const Matrix&)> f_hess([&](const Matrix& x) -> Matrix {
		evals_f_hess++;
		return fun::fourth::hessian(x);
	});

	const int    iter_max{100000};
	const Matrix start1{{3}, {3}};
	const Matrix start2{{1}, {2}};
	const double epsilon{1e-6};
	Matrix min(1, 1);

	evals_f      = 0;
	evals_f_grad = 0;
	evals_f_hess = 0;
	min = Optimization::newtonRaphson(f, f_grad, f_hess, start1, epsilon, iter_max);
	std::cout << "Newton-Raphson:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start1) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;
	std::cout << "hessian evals: " << evals_f_hess << std::endl << std::endl;

	evals_f      = 0;
	evals_f_grad = 0;
	evals_f_hess = 0;
	min = Optimization::newtonRaphson(f, f_grad, f_hess, start2, epsilon, iter_max);
	std::cout << "Newton-Raphson:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start2) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;
	std::cout << "hessian evals: " << evals_f_hess << std::endl << std::endl;

	evals_f      = 0;
	evals_f_grad = 0;
	evals_f_hess = 0;
	min = Optimization::newtonRaphsonGolden(f, f_grad, f_hess, start1, epsilon, iter_max);
	std::cout << "Newton-Raphson Golden Search:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start1) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;
	std::cout << "hessian evals: " << evals_f_hess << std::endl << std::endl;

	evals_f      = 0;
	evals_f_grad = 0;
	evals_f_hess = 0;
	min = Optimization::newtonRaphsonGolden(f, f_grad, f_hess, start2, epsilon, iter_max);
	std::cout << "Newton-Raphson Golden Search:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start2) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;
	std::cout << "hessian evals: " << evals_f_hess << std::endl << std::endl;
	
	std::cout << std::endl;
}

void zad4() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f;
	int evals_f_grad;

	std::function<Matrix(const Matrix&)> f([&](const Matrix& x) -> Matrix {
		evals_f++;
		const double x1{x[0][0]};
		const double x2{x[1][0]};
		return Matrix{{10 * (x2 - pow(x1, 2))}, {1 - x1}};
	});
	std::function<Matrix(const Matrix&)> f_grad([&](const Matrix& x) -> Matrix {
		evals_f_grad++;
		const double x1{x[0][0]};
		return Matrix{
			{-20 * x1, 10},
			{-1      , 0 }
		};
	});

	const int    iter_max{100000};
	const Matrix start{{-1.9}, {2}};
	const double epsilon{1e-6};
	Matrix min(1, 1);

	evals_f      = 0;
	evals_f_grad = 0;
	min = Optimization::gaussNewton(f, f_grad, start, epsilon, iter_max);
	std::cout << "Gauss-Newton:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;

	std::cout << std::endl;
}

void zad5() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f;
	int evals_f_grad;

	std::function<Matrix(const Matrix&)> f([&](const Matrix& x) -> Matrix {
		evals_f++;
		const double x1{x[0][0]};
		const double x2{x[1][0]};
		return Matrix{
			{pow(x1, 2) + pow(x2, 2) - 1},
			{x2 - pow(x1, 2)            }
		};
	});
	std::function<Matrix(const Matrix&)> f_grad([&](const Matrix& x) -> Matrix {
		evals_f_grad++;
		const double x1{x[0][0]};
		const double x2{x[1][0]};
		return Matrix{
			{2 * x1 , 2 * x2},
			{-2 * x1, 1     }
		};
	});

	const int    iter_max{100000};
	const double epsilon{1e-6};
	Matrix min(1, 1);

	const Matrix starts[]{Matrix{{-2}, {2}}, Matrix{{2}, {2}}, Matrix{{2}, {-2}}};

	for (int i{0}; i < 3; i++) {
		evals_f      = 0;
		evals_f_grad = 0;
		const Matrix& start{starts[i]};
		try {
			min = Optimization::gaussNewton(f, f_grad, start, epsilon, iter_max);
		} catch (const std::exception& e) {
			std::cout << "caught an exception: " << e.what() << std::endl;
		}
		std::cout << "Gauss-Newton:" << std::endl;
		std::cout << "start: " << VectorUtils::format(start) << std::endl;
		std::cout << "min: " << VectorUtils::format(min) << std::endl;
		std::cout << "function evals: " << evals_f << std::endl;
		std::cout << "gradient evals: " << evals_f_grad << std::endl << std::endl;
	}

	std::cout << std::endl;
}

void zad6() {
	std::cout << __func__ << std::endl << std::endl;

	int evals_f;
	int evals_f_grad;

	const Matrix pairs{{1, 3}, {2, 4}, {3, 4}, {5, 5}, {6, 6}, {7, 8}};

	std::function<Matrix(const Matrix&)> f([&](const Matrix& x) -> Matrix {
		evals_f++;
		const double x1{x[0][0]};
		const double x2{x[1][0]};
		const double x3{x[2][0]};
		Matrix ret(pairs.getRows(), 1);
		for (int i{0}; i < pairs.getRows(); i++) {
			const double t{pairs[i][0]};
			const double y{pairs[i][1]};
			ret[i][0] = x1 * exp(x2 * t) + x3 - y;
		}
		return ret;
	});
	std::function<Matrix(const Matrix&)> f_grad([&](const Matrix& x) -> Matrix {
		evals_f_grad++;
		const double x1{x[0][0]};
		const double x2{x[1][0]};
		Matrix ret(pairs.getRows(), 3);
		for (int i{0}; i < pairs.getRows(); i++) {
			const double t{pairs[i][0]};
			ret[i][0] = exp(x2 * t);
			ret[i][1] = x1 * exp(x2 * t) * t;
			ret[i][2] = 1;
		}
		return ret;
	});

	const int    iter_max{100000};
	const double epsilon{1e-6};
	Matrix min(1, 1);
	const Matrix start{{1}, {1}, {1}};

	min = Optimization::gaussNewtonGolden(f, f_grad, start, epsilon, iter_max);
	std::cout << "Gauss-Newton Golden Search:" << std::endl;
	std::cout << "start: " << VectorUtils::format(start) << std::endl;
	std::cout << "min: " << VectorUtils::format(min) << std::endl;
	std::cout << "function evals: " << evals_f << std::endl;
	std::cout << "gradient evals: " << evals_f_grad << std::endl;

	std::cout << std::endl;
}

int main(int argc, char* argv[]) {
	std::cout << std::boolalpha;

	if (argc != 2) {
		throw std::invalid_argument("there should be one argument to the program: the assignment id (1..6)");
	}
	const int id{atoi(argv[1])};

	switch (id) {
		case 1: zad1(); break;
		case 2: zad2(); break;
		case 3: zad3(); break;
		case 4: zad4(); break;
		case 5: zad5(); break;
		case 6: zad6(); break;
		default: throw std::invalid_argument("invalid assignment id: should be 1..6");
	}

	return 0;
}
