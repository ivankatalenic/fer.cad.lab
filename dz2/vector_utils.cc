#include "vector_utils.hh"

#include <stdexcept>
#include <cmath>

Matrix VectorUtils::identity_vector(int pos, int size) {
	if (pos < 0) {
		throw std::invalid_argument("the position can't be negative");
	}
	if (size < 1) {
		throw std::invalid_argument("the size can't be less than one");
	}
	if (pos >= size) {
		throw std::out_of_range("the position is greater than the size");
	}

	Matrix ret_v(size, 1);
	ret_v[pos][0] = 1.0;
	return ret_v;
}

double VectorUtils::length(const Matrix& m) {
	double sum{0.0};
	for (int i{0}; i < m.getRows(); ++i) {
		const double comp{m[i][0]};
		sum += comp * comp;
	}
	return std::sqrt(sum);
}

std::string VectorUtils::format(const Matrix& v) {
	const int n{v.getRows()};

	std::string ret;
	ret += "(";
	for (int i{0}; i < n; ++i) {
		ret += std::to_string(v[i][0]);
		if (i != n - 1) {
			ret += ", ";
		}
	}
	ret += ")";

	return ret;
}
