#include "num_int.hh"

#include "vector_utils.hh"
#include "matrix_utils.hh"
#include "exporter_file.hh"

#include <iostream>
#include <cmath>
#include <memory>

void print_res(const std::string& str, const Matrix& v) {
	std::cout << str;
	std::cout << ": [";
	for (int i{0}; i < v.getRows(); ++i) {
		std::cout << v[i][0] << " ";
	}
	std::cout << "]" << std::endl;
}

void zad1() {
	std::cout << __func__ << std::endl;

	const Matrix A{{0, 1}, {-1, 0}};
	const Matrix B{{0, 0}, {0, 0}};
	const Matrix initial{{1}, {1}};
	const double period{0.01};
	const double t_max{10};

	std::shared_ptr<ExporterFile> file;

	Matrix error(2, 1);
	auto callback{[&](double t, const Matrix& v) {
		Matrix valid{{std::cos(t) + std::sin(t)}, {std::cos(t) - std::sin(t)}};
		error = error + MatrixUtils::abs(v - valid);

		file->write_point(t, v);
	}};

	file = std::make_shared<ExporterFile>("output/zad1_euler.txt");
	const auto r1{NumInt::euler(A, B, initial, period, t_max, false, callback)};
	print_res("euler", r1);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

	file = std::make_shared<ExporterFile>("output/zad1_euler_inverted.txt");
	const auto r2{NumInt::euler_inverted(A, B, initial, period, t_max, false, callback)};
	print_res("euler inverted", r2);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

	file = std::make_shared<ExporterFile>("output/zad1_trapeze.txt");
	const auto r3{NumInt::trapeze(A, B, initial, period, t_max, false, callback)};
	print_res("trapeze", r3);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

	file = std::make_shared<ExporterFile>("output/zad1_runge_kutta.txt");
	const auto r4{NumInt::runge_kutta(A, B, initial, period, t_max, false, callback)};
	print_res("runge-kutta", r4);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

	file = std::make_shared<ExporterFile>("output/zad1_predictor_corrector_euler_euler_inverted.txt");
	const auto r5{NumInt::predictor_corrector(Predictor::euler, Corrector::euler_inverted, A, B, initial, period, t_max, false, 2, callback)};
	print_res("predictor euler, corrector euler inverted", r5);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

	file = std::make_shared<ExporterFile>("output/zad1_predictor_corrector_euler_trapeze.txt");
	const auto r6{NumInt::predictor_corrector(Predictor::euler, Corrector::trapeze, A, B, initial, period, t_max, false, 1, callback)};
	print_res("predictor euler, corrector trapeze", r6);
	print_res("error", error);
	error = Matrix(2, 1);
	std::cout << std::endl;

}

void zad2() {
	std::cout << __func__ << std::endl;

	std::shared_ptr<ExporterFile> file;
	auto callback{[&](double t, const Matrix& v) {
		file->write_point(t, v);
	}};

	const Matrix A{{0, 1}, {-200, -102}};
	const Matrix B{{0, 0}, {0, 0}};
	const Matrix initial{{1}, {-2}};
	const double period{0.1};
	const double t_max{1};

	file = std::make_shared<ExporterFile>("output/zad2_euler.txt");
	const auto r1{NumInt::euler(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad2_euler_inverted.txt");
	const auto r2{NumInt::euler_inverted(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad2_trapeze.txt");
	const auto r3{NumInt::trapeze(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad2_runge_kutta.txt");
	const auto r4{NumInt::runge_kutta(A, B, initial, 0.01, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad2_predictor_corrector_euler_euler_inverted.txt");
	const auto r5{NumInt::predictor_corrector(Predictor::euler, Corrector::euler_inverted, A, B, initial, 0.01, t_max, false, 2, callback)};

	file = std::make_shared<ExporterFile>("output/zad2_predictor_corrector_euler_trapeze.txt");
	const auto r6{NumInt::predictor_corrector(Predictor::euler, Corrector::trapeze, A, B, initial, period, t_max, false, 1, callback)};


	print_res("euler", r1);
	print_res("euler inverted", r2);
	print_res("trapeze", r3);
	print_res("runge-kutta", r4);
	print_res("predictor euler, corrector euler inverted", r5);
	print_res("predictor euler, corrector trapeze", r6);
}

void zad3() {
	std::cout << __func__ << std::endl;

	std::shared_ptr<ExporterFile> file;
	auto callback{[&](double t, const Matrix& v) {
		file->write_point(t, v);
	}};

	const Matrix A{{0, -2}, {1, -3}};
	const Matrix B{{2, 0}, {0, 3}};
	const Matrix initial{{1}, {3}};
	const double period{0.01};
	const double t_max{10};

	file = std::make_shared<ExporterFile>("output/zad3_euler.txt");
	const auto r1{NumInt::euler(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad3_euler_inverted.txt");
	const auto r2{NumInt::euler_inverted(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad3_trapeze.txt");
	const auto r3{NumInt::trapeze(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad3_runge_kutta.txt");
	const auto r4{NumInt::runge_kutta(A, B, initial, period, t_max, false, callback)};

	file = std::make_shared<ExporterFile>("output/zad3_predictor_corrector_euler_euler_inverted.txt");
	const auto r5{NumInt::predictor_corrector(Predictor::euler, Corrector::euler_inverted, A, B, initial, period, t_max, false, 2, callback)};

	file = std::make_shared<ExporterFile>("output/zad3_predictor_corrector_euler_trapeze.txt");
	const auto r6{NumInt::predictor_corrector(Predictor::euler, Corrector::trapeze, A, B, initial, period, t_max, false, 1, callback)};


	print_res("euler", r1);
	print_res("euler inverted", r2);
	print_res("trapeze", r3);
	print_res("runge-kutta", r4);
	print_res("predictor euler, corrector euler inverted", r5);
	print_res("predictor euler, corrector trapeze", r6);
}

void zad4() {
	std::cout << __func__ << std::endl;

	std::shared_ptr<ExporterFile> file;
	auto callback{[&](double t, const Matrix& v) {
		file->write_point(t, v);
	}};

	const Matrix A{{1, -5}, {1, -7}};
	const Matrix B{{5, 0}, {0, 3}};
	const Matrix initial{{-1}, {3}};
	const double period{0.01};
	const double t_max{1};

	file = std::make_shared<ExporterFile>("output/zad4_euler.txt");
	const auto r1{NumInt::euler(A, B, initial, period, t_max, true, callback)};

	file = std::make_shared<ExporterFile>("output/zad4_euler_inverted.txt");
	const auto r2{NumInt::euler_inverted(A, B, initial, period, t_max, true, callback)};

	file = std::make_shared<ExporterFile>("output/zad4_trapeze.txt");
	const auto r3{NumInt::trapeze(A, B, initial, period, t_max, true, callback)};

	file = std::make_shared<ExporterFile>("output/zad4_runge_kutta.txt");
	const auto r4{NumInt::runge_kutta(A, B, initial, period, t_max, true, callback)};

	file = std::make_shared<ExporterFile>("output/zad4_predictor_corrector_euler_euler_inverted.txt");
	const auto r5{NumInt::predictor_corrector(Predictor::euler, Corrector::euler_inverted, A, B, initial, period, t_max, true, 2, callback)};

	file = std::make_shared<ExporterFile>("output/zad4_predictor_corrector_euler_trapeze.txt");
	const auto r6{NumInt::predictor_corrector(Predictor::euler, Corrector::trapeze, A, B, initial, period, t_max, true, 1, callback)};


	print_res("euler", r1);
	print_res("euler inverted", r2);
	print_res("trapeze", r3);
	print_res("runge-kutta", r4);
	print_res("predictor euler, corrector euler inverted", r5);
	print_res("predictor euler, corrector trapeze", r6);
}

int main() {
	zad1();
	// zad2();
	// zad3();
	// zad4();

	return 0;
}
