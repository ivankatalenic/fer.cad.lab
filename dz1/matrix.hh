#pragma once

#include <initializer_list>

class Matrix {
public:
	Matrix(int rows, int columns);
	Matrix(std::initializer_list<class Row> l);
	
	Matrix(const Matrix& other);
	Matrix(Matrix&& other);

	int getRows() const;
	int getColumns() const;

	~Matrix();

	Matrix& operator=(const Matrix& other);
	Matrix& operator=(Matrix&& other);

	Matrix operator+(const Matrix& other) const;
	Matrix operator-(const Matrix& other) const;
	Matrix operator*(const Matrix& other) const;
	Matrix operator~() const;
	Matrix operator*(const double& s) const;

	class RowView operator[](int row_index);
	const class RowView operator[](int row_index) const;

	bool operator==(const Matrix& other) const;

private:
	int rows;
	int columns;

	double* values;
};

class Row {
public:
	Row(std::initializer_list<double> l);
	~Row();

	int getSize() const;

	double& operator[](int index);
	double  operator[](int index) const;
private:
	double*   values;
	const int size;
};

class RowView {
public:
	RowView(double* values, int size);

	double& operator[](int index);
	double  operator[](int index) const;
private:
	double*   values;
	const int size;
};

Matrix operator*(const double& s, const Matrix& m);
