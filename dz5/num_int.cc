#include "num_int.hh"

#include "float_utils.hh"
#include "vector_utils.hh"
#include "matrix_utils.hh"
#include "matrix_ops.hh"

#include <stdexcept>

Matrix NumInt::euler(
	const Matrix& A,
	const Matrix& B,
	const Matrix& initial,
	const double period,
	const double t_max,
	const bool is_time_dependent,
	std::function<void(double t, const Matrix& x)> callback
) {
	if (period < 0) {
		throw std::invalid_argument("the period can't be negative");
	}
	if (t_max < 0) {
		throw std::invalid_argument("the integration limit can't be negative");
	}

	if (callback != nullptr)
		callback(0, initial);

	double t{period};
	Matrix x(initial);
	Matrix time_matrix{VectorUtils::uniform_vector(B.getColumns(), 1)};

	while (t < t_max || FloatUtils::is_equal(t, t_max, 1e-6)) {
		if (is_time_dependent) {
			time_matrix = VectorUtils::uniform_vector(B.getColumns(), t - period);
		}

		const Matrix derivative{A * x + B * time_matrix};
		x = x + period * derivative;

		if (callback != nullptr)
			callback(t, x);

		t += period;
	}
	if (!FloatUtils::is_equal(t - period, t_max, 1e-6)) {
		const double last_period{t_max - (t - period)};
		t = t_max;

		if (is_time_dependent) {
			time_matrix = VectorUtils::uniform_vector(B.getColumns(), t - last_period);
		}

		const Matrix derivative{A * x + B * time_matrix};
		x = x + last_period * derivative;

		if (callback != nullptr)
			callback(t, x);

		t += last_period;
	}

	return x;
}

Matrix NumInt::euler_inverted(
	const Matrix& A,
	const Matrix& B,
	const Matrix& initial,
	const double period,
	const double t_max,
	const bool is_time_dependent,
	std::function<void(double t, const Matrix& x)> callback
) {
	if (period < 0) {
		throw std::invalid_argument("the period can't be negative");
	}
	if (t_max < 0) {
		throw std::invalid_argument("the integration limit can't be negative");
	}

	const Matrix P{MatrixOps::inverse(MatrixUtils::identity(A.getRows()) - A * period)};
	const Matrix Q{P * period * B};

	if (callback != nullptr)
		callback(0, initial);

	double t{period};
	Matrix x(initial);
	Matrix time_matrix{VectorUtils::uniform_vector(B.getColumns(), 1)};
	while (t < t_max || FloatUtils::is_equal(t, t_max, 1e-6)) {
		if (is_time_dependent) {
			time_matrix = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		x = P * x + Q * time_matrix;

		if (callback != nullptr)
			callback(t, x);

		t += period;
	}
	if (!FloatUtils::is_equal(t - period, t_max, 1e-6)) {
		const double last_period{t_max - (t - period)};
		t = t_max;

		if (is_time_dependent) {
			time_matrix = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		x = P * x + Q * time_matrix;

		if (callback != nullptr)
			callback(t, x);

		t += last_period;
	}

	return x;
}

Matrix NumInt::trapeze(
	const Matrix& A,
	const Matrix& B,
	const Matrix& initial,
	const double period,
	const double t_max,
	const bool is_time_dependent,
	std::function<void(double t, const Matrix& x)> callback
) {
	if (period < 0) {
		throw std::invalid_argument("the period can't be negative");
	}
	if (t_max < 0) {
		throw std::invalid_argument("the integration limit can't be negative");
	}

	const Matrix V{MatrixOps::inverse(MatrixUtils::identity(A.getRows()) - A * 0.5 * period)};
	const Matrix R{V * (MatrixUtils::identity(A.getRows()) + A * 0.5 * period)};
	const Matrix S{V * 0.5 * period * B};

	if (callback != nullptr)
		callback(0, initial);

	double t{period};
	Matrix x(initial);
	Matrix time_matrix1{VectorUtils::uniform_vector(B.getColumns(), 1)};
	Matrix time_matrix2{VectorUtils::uniform_vector(B.getColumns(), 1)};
	while (t < t_max || FloatUtils::is_equal(t, t_max, 1e-6)) {
		if (is_time_dependent) {
			time_matrix1 = VectorUtils::uniform_vector(B.getColumns(), t - period);
			time_matrix2 = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		x = R * x + S * (time_matrix1 + time_matrix2);

		if (callback != nullptr)
			callback(t, x);

		t += period;
	}
	if (!FloatUtils::is_equal(t - period, t_max, 1e-6)) {
		const double last_period{t_max - (t - period)};
		t = t_max;

		if (is_time_dependent) {
			time_matrix1 = VectorUtils::uniform_vector(B.getColumns(), t - last_period);
			time_matrix2 = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		x = R * x + S * (time_matrix1 + time_matrix2);

		if (callback != nullptr)
			callback(t, x);

		t += last_period;
	}

	return x;
}

Matrix NumInt::runge_kutta(
	const Matrix& A,
	const Matrix& B,
	const Matrix& initial,
	const double period,
	const double t_max,
	const bool is_time_dependent,
	std::function<void(double t, const Matrix& x)> callback
) {
	if (period < 0) {
		throw std::invalid_argument("the period can't be negative");
	}
	if (t_max < 0) {
		throw std::invalid_argument("the integration limit can't be negative");
	}

	const Matrix V{MatrixOps::inverse(MatrixUtils::identity(A.getRows()) - A * 0.5 * period)};
	const Matrix R{V * (MatrixUtils::identity(A.getRows()) + A * 0.5 * period)};
	const Matrix S{V * 0.5 * period * B};

	if (callback != nullptr)
		callback(0, initial);

	double t{period};
	Matrix x(initial);
	Matrix time_matrix1{VectorUtils::uniform_vector(B.getColumns(), 1)};
	Matrix time_matrix2{VectorUtils::uniform_vector(B.getColumns(), 1)};
	Matrix time_matrix3{VectorUtils::uniform_vector(B.getColumns(), 1)};
	while (t < t_max || FloatUtils::is_equal(t, t_max, 1e-6)) {
		if (is_time_dependent) {
			time_matrix1 = VectorUtils::uniform_vector(B.getColumns(), t - period);
			time_matrix2 = VectorUtils::uniform_vector(B.getColumns(), t - 0.5 * period);
			time_matrix3 = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		const Matrix m1{A * x + B * time_matrix1};
		const Matrix m2{A * (x + 0.5 * period * m1) + B * time_matrix2};
		const Matrix m3{A * (x + 0.5 * period * m2) + B * time_matrix2};
		const Matrix m4{A * (x + period * m3) + B * time_matrix3};

		x = x + (period / 6) * (m1 + 2*m2 + 2*m3 + m4);

		if (callback != nullptr)
			callback(t, x);

		t += period;
	}
	if (!FloatUtils::is_equal(t - period, t_max, 1e-6)) {
		const double last_period{t_max - (t - period)};
		t = t_max;

		if (is_time_dependent) {
			time_matrix1 = VectorUtils::uniform_vector(B.getColumns(), t - period);
			time_matrix2 = VectorUtils::uniform_vector(B.getColumns(), t - 0.5 * period);
			time_matrix3 = VectorUtils::uniform_vector(B.getColumns(), t);
		}

		const Matrix m1{A * x + B * time_matrix1};
		const Matrix m2{A * (x + 0.5 * period * m1) + B * time_matrix2};
		const Matrix m3{A * (x + 0.5 * period * m2) + B * time_matrix2};
		const Matrix m4{A * (x + period * m3) + B * time_matrix3};

		x = x + (period / 6) * (m1 + 2*m2 + 2*m3 + m4);

		if (callback != nullptr)
			callback(t, x);

		t += last_period;
	}

	return x;
}

Matrix NumInt::predictor_corrector(
	std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, double, double, bool)> predictor,
	std::function<Matrix(const Matrix&, const Matrix&, const Matrix&, const Matrix&, double, double, bool)> corrector,
	const Matrix& A,
	const Matrix& B,
	const Matrix& initial,
	const double period,
	const double t_max,
	const bool is_time_dependent,
	const int iter_max,
	std::function<void(double t, const Matrix& x)> callback
) {
	if (predictor == nullptr) {
		throw std::invalid_argument("the predictor is not specified");
	}
	if (corrector == nullptr) {
		throw std::invalid_argument("the corrector is not specified");
	}
	if (period < 0) {
		throw std::invalid_argument("the period can't be negative");
	}
	if (t_max < 0) {
		throw std::invalid_argument("the integration limit can't be negative");
	}
	if (iter_max < 1) {
		throw std::invalid_argument("the number of iterations should be greater than 1");
	}

	if (callback != nullptr)
		callback(0, initial);
	
	Matrix x(initial);
	double t{period};
	while (t < t_max || FloatUtils::is_equal(t, t_max, 1e-6)) {
		Matrix x_predicted{predictor(x, A, B, t, period, is_time_dependent)};
		for (int i{0}; i < iter_max; ++i) {
			x_predicted = corrector(x, x_predicted, A, B, t, period, is_time_dependent);
		}
		x = x_predicted;

		if (callback != nullptr)
			callback(t, x);

		t += period;
	}
	if (!FloatUtils::is_equal(t - period, t_max, 1e-6)) {
		const double last_period{t_max - (t - period)};
		t = t_max;

		Matrix x_predicted{predictor(x, A, B, t, last_period, is_time_dependent)};
		for (int i{0}; i < iter_max; ++i) {
			x_predicted = corrector(x, x_predicted, A, B, t, last_period, is_time_dependent);
		}
		x = x_predicted;

		if (callback != nullptr)
			callback(t, x);

		t += last_period;
	}

	return x;
}

Matrix Predictor::euler(
	const Matrix& x,
	const Matrix& A,
	const Matrix& B,
	double t,
	double period,
	bool is_time_dependent
) {
	double val{t - period};
	if (!is_time_dependent) {
		val = 1;
	}
	const Matrix time_matrix{VectorUtils::uniform_vector(B.getColumns(), val)};

	const Matrix derivative{A * x + B * time_matrix};
	return x + period * derivative;
}

Matrix Corrector::euler_inverted(
	const Matrix& x,
	const Matrix& x_predicted,
	const Matrix& A,
	const Matrix& B,
	double t,
	double period,
	bool is_time_dependent
) {
	double val{t};
	if (!is_time_dependent) {
		val = 1;
	}
	const Matrix time_matrix{VectorUtils::uniform_vector(B.getColumns(), val)};
	
	const Matrix derivative{A * x_predicted + B * time_matrix};
	return x + period * derivative;
}

Matrix Corrector::trapeze(
	const Matrix& x,
	const Matrix& x_predicted,
	const Matrix& A,
	const Matrix& B,
	double t,
	double period,
	bool is_time_dependent
) {
	Matrix time_matrix1{VectorUtils::uniform_vector(B.getColumns(), 1)};
	Matrix time_matrix2{VectorUtils::uniform_vector(B.getColumns(), 1)};
	if (is_time_dependent) {
		time_matrix1 = VectorUtils::uniform_vector(B.getColumns(), t - period);
		time_matrix2 = VectorUtils::uniform_vector(B.getColumns(), t);
	}

	const Matrix derv1{A * x + B * time_matrix1};
	const Matrix derv2{A * x_predicted + B * time_matrix2};
	return x + 0.5 * period * (derv1 + derv2);
}
