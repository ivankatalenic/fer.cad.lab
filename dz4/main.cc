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

	const Matrix bound_lower{{-100}, {-100}};
	const Matrix bound_upper{{100}, {100}};
	std::function<double(const Matrix&)> constr1{[](const Matrix& v) -> double {
		const double x1{v[0][0]}, x2{v[1][0]};
		return x2 - x1;
	}};
	std::function<double(const Matrix&)> constr2{[](const Matrix& v) -> double {
		const double x1{v[0][0]};
		return 2 - x1;
	}};
	const double epsilon{1e-6};
	const double alpha{1.3};
	const int iter_max{100000};
	Matrix min(1, 1);

	min = Optimization::boxMethod(
		fun::first::f,
		bound_lower,
		bound_upper,
		{constr1, constr2},
		fun::first::start,
		epsilon,
		alpha,
		iter_max
	);
	std::cout << "Function 1 optimum: " << VectorUtils::format(min) << std::endl;

	min = Optimization::boxMethod(
		fun::second::f,
		bound_lower,
		bound_upper,
		{constr1, constr2},
		fun::second::start,
		epsilon,
		alpha,
		iter_max
	);
	std::cout << "Function 2 optimum: " << VectorUtils::format(min) << std::endl;
	
	std::cout << std::endl;
}

void zad2() {
	std::cout << __func__ << std::endl << std::endl;

	std::function<double(const Matrix&)> constr1{[](const Matrix& v) -> double {
		const double x1{v[0][0]}, x2{v[1][0]};
		return x2 - x1;
	}};
	std::function<double(const Matrix&)> constr2{[](const Matrix& v) -> double {
		const double x1{v[0][0]};
		return 2 - x1;
	}};
	const double epsilon{1e-6};
	const int iter_max{100000};
	int t = 1;
	Matrix min(1, 1);

	Matrix start{fun::first::start};
	std::cout << "Start: " << VectorUtils::format(start) << std::endl;
	min = Optimization::transMethod(fun::first::f, {constr1, constr2}, {}, start, epsilon, t, iter_max);
	std::cout << "Function 1 optimum: " << VectorUtils::format(min) << std::endl;

	std::cout << '\n';

	start = {{0.1}, {2}};
	std::cout << "Start: " << VectorUtils::format(start) << std::endl;
	min = Optimization::transMethod(fun::first::f, {constr1, constr2}, {}, start, epsilon, t, iter_max);
	std::cout << "Function 1 optimum: " << VectorUtils::format(min) << std::endl;

	std::cout << '\n';

	start = fun::second::start;
	std::cout << "Start: " << VectorUtils::format(start) << std::endl;
	min = Optimization::transMethod(fun::second::f, {constr1, constr2}, {}, start, epsilon, t, iter_max);
	std::cout << "Function 2 optimum: " << VectorUtils::format(min) << std::endl;
	
	std::cout << std::endl;
}

void zad3() {
	std::cout << __func__ << std::endl << std::endl;

	std::function<double(const Matrix&)> constr1{[](const Matrix& v) -> double {
		const double x1{v[0][0]}, x2{v[1][0]};
		return 3 - x1 - x2;
	}};
	std::function<double(const Matrix&)> constr2{[](const Matrix& v) -> double {
		const double x1{v[0][0]}, x2{v[1][0]};
		return 3 + 1.5 * x1 - x2;
	}};
	std::function<double(const Matrix&)> constr3{[](const Matrix& v) -> double {
		const double x2{v[1][0]};
		return x2 - 1;
	}};
	const double epsilon{1e-6};
	const int iter_max{100000};
	Matrix min(1, 1);

	int t = 1;
	Matrix start{fun::fourth::start};
	std::cout << "Start: " << VectorUtils::format(start) << std::endl;
	min = Optimization::transMethod(fun::fourth::f, {constr1, constr2}, {constr3}, start, epsilon, t, iter_max);
	std::cout << "Optimum: " << VectorUtils::format(min) << std::endl;

	t = 1;
	start = {{5}, {5}};
	std::cout << "Start: " << VectorUtils::format(start) << std::endl;
	min = Optimization::transMethod(fun::fourth::f, {constr1, constr2}, {constr3}, start, epsilon, t, iter_max);
	std::cout << "Optimum: " << VectorUtils::format(min) << std::endl;
	
	std::cout << std::endl;
}

int main(int argc, char* argv[]) {
	std::cout << std::boolalpha;

	if (argc != 2) {
		throw std::invalid_argument("there should be one argument to the program: the assignment id (1..3)");
	}
	const int id{atoi(argv[1])};

	switch (id) {
		case 1: zad1(); break;
		case 2: zad2(); break;
		case 3: zad3(); break;
		default: throw std::invalid_argument("invalid assignment id: should be 1..3");
	}

	return 0;
}
