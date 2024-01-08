#include "matrix.hh"
#include "matrix_utils.hh"
#include "matrix_ops.hh"

#include <iostream>

void zad1() {
	std::cout << __func__ << std::endl;

	const double a{25e9};

	double m{1e-30};

	const double a1{a / m * m};
	std::cout << "a == (a / m * m): " << (a == a1) << std::endl;
	
	std::cout << std::endl;
}

void zad2() {
	std::cout << __func__ << std::endl;

	Matrix A{{3, 9, 6}, {4, 12, 12}, {1, -1, 1}};
	Matrix b{{12}, {12}, {1}};

	Matrix inv(MatrixOps::inverse(A));
	MatrixUtils::printToConsole(inv * b);
	
	std::cout << std::endl;
}

void zad3() {
	std::cout << __func__ << std::endl;

	Matrix A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

	const Matrix lu{MatrixOps::decomposeLU(A)};
	const auto [lup, perm] = MatrixOps::decomposeLUP(A);

	std::cout << "lu:" << std::endl;
	MatrixUtils::printToConsole(lu);
	std::cout << std::endl;
	std::cout << "lup:" << std::endl;	
	MatrixUtils::printToConsole(lup);
	std::cout << std::endl;
	std::cout << "perm:" << std::endl;
	MatrixUtils::printToConsole(perm);
	std::cout << std::endl;
}

void zad4() {
	std::cout << __func__ << std::endl;

	Matrix A{{1e-6, 3e6, 2e6}, {1e6, 2e6, 3e6}, {2e6, 1e6, 2e6}};
	Matrix b{{12e6 + 1e-6}, {14e6}, {1e7}};
	
	Matrix lu{MatrixOps::decomposeLU(A)};
	auto [lup, perm] = MatrixOps::decomposeLUP(A);

	Matrix y_lu{MatrixOps::subForward(lu, b)};
	Matrix y_lup{MatrixOps::subForward(lup, perm * b)};

	Matrix x_lu{MatrixOps::subBackward(lu, y_lu)};
	Matrix x_lup{MatrixOps::subBackward(lup, y_lup)};

	std::cout << "lu:" << std::endl;
	MatrixUtils::printToConsole(x_lu);
	std::cout << std::endl;
	std::cout << "lup:" << std::endl;
	MatrixUtils::printToConsole(x_lup);
	std::cout << std::endl;
}

void zad5() {
	std::cout << __func__ << std::endl;

	Matrix A{{0, 1, 2}, {2, 0, 3}, {3, 5, 1}};
	Matrix b{{6}, {9}, {3}};
	
	auto [lup, perm] = MatrixOps::decomposeLUP(A);

	Matrix y_lup{MatrixOps::subForward(lup, perm * b)};

	Matrix x_lup{MatrixOps::subBackward(lup, y_lup)};

	std::cout << "solution:" << std::endl;
	MatrixUtils::printToConsole(x_lup);
	std::cout << std::endl;
}

void zad6() {
	std::cout << __func__ << std::endl;

	Matrix A{{4e9, 1e9, 3e9}, {4, 2, 7}, {3e-10, 5e-10, 2e-10}};
	Matrix b{{9e9}, {15}, {15e-10}};
	
	auto [lup, perm] = MatrixOps::decomposeLUP(A);

	Matrix y_lup{MatrixOps::subForward(lup, perm * b)};

	Matrix x_lup{MatrixOps::subBackward(lup, y_lup)};

	std::cout << "lup:" << std::endl;
	MatrixUtils::printToConsole(x_lup);
	std::cout << std::endl;
}

void zad7() {
	std::cout << __func__ << std::endl;

	Matrix A{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

	std::cout << "solution:" << std::endl;
	MatrixUtils::printToConsole(MatrixOps::inverse(A));
	std::cout << std::endl;
}

void zad8() {
	std::cout << __func__ << std::endl;

	Matrix A{{4, -5, -2}, {5, -6, -2}, {-8, 9, 3}};

	std::cout << "solution:" << std::endl;
	MatrixUtils::printToConsole(MatrixOps::inverse(A));
	std::cout << std::endl;
}

void zad9() {
	std::cout << __func__ << std::endl;

	Matrix A{{4, -5, -2}, {5, -6, -2}, {-8, 9, 3}};

	std::cout << "det: " << MatrixOps::determinant(A) << std::endl;
	std::cout << std::endl;
}

void zad10() {
	std::cout << __func__ << std::endl;

	Matrix A{{3, 9, 6}, {4, 12, 12}, {1, -1, 1}};

	std::cout << "det: " << MatrixOps::determinant(A) << std::endl;
	std::cout << std::endl;
}

int main() {
	std::cout << std::boolalpha;

	zad1();
	zad2();
	zad3();
	zad4();
	zad5();
	zad6();
	// zad7(); // singular
	zad8();
	zad9();
	zad10();

	return 0;
}
