#include "matrix.hh"

#include <stdexcept>
#include <cmath>

Matrix::Matrix(int r, int c):
	rows{r},
	columns{c},
	values{new double[rows * columns]{0.0}}
{
	if (r < 1 || c < 1) {
		throw std::invalid_argument("the provided matrix dimensions are invalid: they must be 1 or higher");
	}
}

Matrix::Matrix(std::initializer_list<Row> l) {
	const int cnt_row{static_cast<int>(l.size())};
	if (cnt_row == 0) {
		throw std::invalid_argument("can't construct a matrix with zero rows");
	}

	int  cnt_column{0};
	bool cnt_column_set{false};
	for (const Row& row : l) {
		if (row.getSize() == 0) {
			throw std::invalid_argument("can't construct a matrix with zero columns");
		}
		if (!cnt_column_set) {
			cnt_column_set = true;
			cnt_column     = row.getSize();
			continue;
		}
		if (row.getSize() != cnt_column) {
			throw std::invalid_argument("can't construct a matrix with inconsistent row sizes");
		}
	}

	rows    = cnt_row;
	columns = cnt_column;
	values  = new double[rows * columns];

	int i{0};
	for (const Row& row : l) {
		for (int j{0}; j < row.getSize(); j++) {
			values[i * columns + j] = row[j];
		}
		i++;
	}
}

Matrix::Matrix(const Matrix& other):
	rows{other.rows},
	columns{other.columns},
	values{new double[rows * columns]}
{
	for (int i{0}; i < rows * columns; ++i) {
		values[i] = other.values[i];
	}
}

Matrix::Matrix(Matrix&& other):
	rows{other.rows},
	columns{other.columns},
	values{other.values}
{
	other.rows = 0;
	other.columns = 0;
	other.values = nullptr;
}

int Matrix::getRows() const {
	return rows;
}

int Matrix::getColumns() const {
	return columns;
}

double Matrix::getLength() const {
	if (this->getColumns() != 1) {
		throw std::domain_error("the matrix is not a vector with one column");
	}

	double square_sum{0.0};
	for (int i{0}; i < rows; ++i) {
		square_sum += values[i] * values[i];
	}

	return std::sqrt(square_sum);
}

Matrix::~Matrix() {
	delete[] values;
}

Matrix& Matrix::operator=(const Matrix& other) {
	if (this == &other) {
		return *this;
	}

	const int size_new{other.rows * other.columns};
	const int size_old{rows * columns};

	if (size_new != size_old) {
		delete[] values;
		values = new double[size_new];
	}

	for (int i{0}; i < size_new; ++i) {
		values[i] = other.values[i];
	}

	rows = other.rows;
	columns = other.columns;

	return *this;
}

Matrix& Matrix::operator=(Matrix&& other) {
	if (this == &other) {
		return *this;
	}

	delete[] values;

	values = other.values;
	rows = other.rows;
	columns = other.columns;

	other.values = nullptr;
	other.rows = 0;
	other.columns = 0;

	return *this;
}

Matrix Matrix::operator+(const Matrix& other) const {
	if (other.rows != rows) {
		throw std::invalid_argument("matrices have incompatible row sizes");
	}
	if (other.columns != columns) {
		throw std::invalid_argument("matrices have incompatible column sizes");
	}

	Matrix result(rows, columns);
	for (int i{0}; i < rows * columns; ++i) {
		result.values[i] = values[i] + other.values[i];
	}

	return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
	if (other.rows != rows) {
		throw std::invalid_argument("matrices have incompatible row sizes");
	}
	if (other.columns != columns) {
		throw std::invalid_argument("matrices have incompatible column sizes");
	}

	Matrix result(rows, columns);
	for (int i{0}; i < rows * columns; ++i) {
		result.values[i] = values[i] - other.values[i];
	}

	return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
	if (columns != other.rows) {
		throw std::invalid_argument("matrices have incompatible sizes");
	}

	Matrix result(rows, other.columns);
	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < other.columns; ++j) {
			double sum{0.0};
			
			for (int k{0}; k < columns; ++k) {
				sum += (*this)[i][k] * other[k][j];
			}

			result[i][j] = sum;
		}
	}

	return result;
}

Matrix Matrix::operator~() const {
	Matrix result(*this);

	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < columns; ++j) {
			result[j][i] = (*this)[i][j];
		}
	}

	return result;
}

Matrix Matrix::operator*(const double& s) const {
	Matrix result(*this);

	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < columns; ++j) {
			result[i][j] *= s;
		}
	}

	return result;
}

RowView Matrix::operator[](int row_index) {
	if (row_index < 0 || row_index > rows - 1) {
		throw std::out_of_range("the provided row index is out of range");
	}
	return RowView(values + row_index * columns, columns);
}

const RowView Matrix::operator[](int row_index) const {
	if (row_index < 0 || row_index > rows - 1) {
		throw std::out_of_range("the provided row index is out of range");
	}
	return RowView(values + row_index * columns, columns);
}

bool Matrix::operator==(const Matrix& other) const {
	if (this == &other) {
		return true;
	}

	if (other.rows != rows || other.columns != columns) {
		return false;
	}

	for (int i{0}; i < rows * columns; ++i) {
		if (values[i] != other.values[i]) {
			return false;
		}
	}

	return true;
}

Matrix operator*(const double& s, const Matrix& m) {
	return m * s;
}

Row::Row(std::initializer_list<double> l):
	values{new double[l.size()]},
	size{static_cast<int>(l.size())}
{
	int i{0};
	for (double v : l) {
		values[i++] = v;
	}
}

int Row::getSize() const {
	return size;
}

double& Row::operator[](int index) {
	if (index < 0 || index > size - 1) {
		throw std::out_of_range("the provided index is out of range");
	}
	return values[index];
}

double  Row::operator[](int index) const {
	if (index < 0 || index > size - 1) {
		throw std::out_of_range("the provided index is out of range");
	}
	return values[index];
}


Row::~Row() {
	delete[] values;
}

RowView::RowView(double* v, int s): values{v}, size{s} {
	if (s < 1) {
		throw std::invalid_argument("the size value must be 1 or higher");
	}
	if (v == nullptr) {
		throw std::invalid_argument("the values pointer must not be null");
	}
}

double& RowView::operator[](int index) {
	if (index < 0 || index > size - 1) {
		throw std::out_of_range("the provided index is out of range");
	}
	return values[index];
}

double RowView::operator[](int index) const {
	if (index < 0 || index > size - 1) {
		throw std::out_of_range("the provided index is out of range");
	}
	return values[index];
}
