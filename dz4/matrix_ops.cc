#include "matrix_ops.hh"

#include <stdexcept>
#include <cmath>
#include <tuple>
#include <optional>
#include <utility>
#include <vector>

static Matrix getIdentity(int rows) {
	Matrix result(rows, rows);

	for (int i{0}; i < rows; ++i) {
		result[i][i] = 1.0;
	}

	return result;
}

static bool isEqual(double a, double b, double epsilon) {
	if (std::abs(a - b) < epsilon) {
		return true;
	}
	return false;
}

static void switchRows(Matrix& m, int first_row, int second_row) {
	if (first_row < 0 || first_row > m.getRows() - 1) {
		throw std::out_of_range("the index for first row is out of range");
	}
	if (second_row < 0 || second_row > m.getRows() - 1) {
		throw std::out_of_range("the index for second row is out of range");
	}

	if (first_row == second_row) {
		return;
	}

	for (int i{0}; i < m.getColumns(); ++i) {
		std::swap(m[first_row][i], m[second_row][i]);
	}
}

static std::tuple<Matrix, std::optional<std::tuple<Matrix, int>>> decompose(const Matrix& m, bool is_pivoting) {
	if (m.getColumns() != m.getRows()) {
		throw std::invalid_argument("the matrix doesn't have the row count equal to the column count");
	}

	Matrix result(m);
	Matrix perm_mat(getIdentity(m.getRows()));
	int    perm_cnt{0};

	const int n{result.getRows()};
	for (int i{0}; i < n - 1; ++i) {
		if (is_pivoting) {
			int index{i};
			for (int j{i + 1}; j < n; ++j) {
				if (std::abs(result[j][i]) > std::abs(result[index][i])) {
					index = j;
				}
			}

			if (isEqual(std::abs(result[index][i]), 0.0, 1e-5)) {
				throw std::invalid_argument("the provided matrix is singular");
			}

			if (index != i) {
				switchRows(result,   index, i);
				switchRows(perm_mat, index, i);
				perm_cnt++;
			}
		}

		for (int j{i + 1}; j < n; ++j) {
			const double denom{result[i][i]};
			if (isEqual(denom, 0.0, 1e-6)) {
				throw std::domain_error("attempted to divide by zero while decomposing a matrix");
			}
			result[j][i] /= denom;
			
			for (int k{i + 1}; k < n; ++k) {
				result[j][k] -= result[j][i] * result[i][k];
			}
		}
	}

	std::optional<std::tuple<Matrix, int>> optional;
	if (is_pivoting) {
		optional = std::make_optional(std::make_tuple(perm_mat, perm_cnt));
	}

	return std::make_tuple(result, optional);
}

static Matrix getColumn(const Matrix& m, const int column_index) {
	if (column_index < 0 || column_index > m.getColumns() - 1) {
		throw std::out_of_range("the column index is out of range");
	}

	Matrix result(m.getRows(), 1);
	for (int i{0}; i < m.getRows(); i++) {
		result[i][0] = m[i][column_index];
	}

	return result;
}

static Matrix mergeColumns(const std::vector<Matrix>& c) {
	if (c.size() == 0) {
		throw std::invalid_argument("can't create a matrix with zero columns");
	}

	const int columns{static_cast<int>(c.size())};
	const int rows{c[0].getRows()};

	Matrix result(rows, columns);

	for (int j{0}; j < columns; ++j) {
		const Matrix& part{c[j]};

		if (part.getColumns() != 1) {
			throw std::invalid_argument("the matrix part doesn't have the exactly one column");
		}
		if (part.getRows() != rows) {
			throw std::invalid_argument("the matrix parts don't all have the same number of rows");
		}
		
		for (int i{0}; i < rows; ++i) {
			result[i][j] = part[i][0];
		}
	}

	return result;
}

Matrix MatrixOps::subForward(const Matrix& L, const Matrix& b) {
	if (L.getRows() != L.getColumns()) {
		throw std::invalid_argument("the matrix must have the row count equal to the column count");
	}
	if (b.getColumns() != 1 || b.getRows() != L.getColumns()) {
		throw std::invalid_argument("the vector doesn't have the compatible dimension");
	}

	Matrix result(b);

	const int n{b.getRows()};
	for (int i{0}; i < n - 1; ++i) {
		for (int j{i + 1}; j < n; ++j) {
			result[j][0] -= L[j][i] * result[i][0];
		}
	}

	return result;
}

Matrix MatrixOps::subBackward(const Matrix& U, const Matrix& y) {
	if (U.getRows() != U.getColumns()) {
		throw std::invalid_argument("the matrix must have the row count equal to the column count");
	}
	if (y.getColumns() != 1 || y.getRows() != U.getColumns()) {
		throw std::invalid_argument("the vector doesn't have the compatible dimension");
	}

	Matrix result(y);

	const int n{y.getRows()};
	for (int i{n - 1}; i >= 0; --i) {
		result[i][0] /= U[i][i];

		for (int j{0}; j < i; ++j) {
			result[j][0] -= U[j][i] * result[i][0];
		}
	}

	return result;
}

Matrix MatrixOps::decomposeLU(const Matrix& m) {
	auto [result, perm] = decompose(m, false);
	return result;
}

std::tuple<Matrix, Matrix> MatrixOps::decomposeLUP(const Matrix& m) {
	auto [result, opt] = decompose(m, true);
	auto [perm_mat, perm_cnt] = opt.value();
	return std::make_tuple(result, perm_mat);
}

Matrix MatrixOps::inverse(const Matrix& m) {
	if (m.getColumns() != m.getRows()) {
		throw std::invalid_argument("the matrix doesn't have the row count equal to the column count");
	}

	const double det{determinant(m)};
	if (isEqual(det, 0.0, 1e-6)) {
		throw std::domain_error("can't find an inverse of a singluar matrix");
	}

	const auto [lu_mat, perm_mat] = decomposeLUP(m);

	std::vector<Matrix> inverse_columns;
	for (int i{0}; i < m.getColumns(); ++i) {
		const Matrix y{subForward(lu_mat, getColumn(perm_mat, i))};
		const Matrix x{subBackward(lu_mat, y)};

		inverse_columns.push_back(x);
	}

	return mergeColumns(inverse_columns);
}

double MatrixOps::determinant(const Matrix& m) {
	if (m.getRows() != m.getColumns()) {
		throw std::invalid_argument("the matrix must have the row number equal to the column number");
	}

	const auto [lu_mat, opt]        = decompose(m, true);
	const auto [perm_mat, perm_cnt] = opt.value();

	double det{1.0};
	for (int i{0}; i < lu_mat.getRows(); ++i) {
		det *= lu_mat[i][i];
	}

	return (perm_cnt % 2 == 0) ? det : -det;
}
