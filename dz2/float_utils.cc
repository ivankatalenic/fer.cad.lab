#include "float_utils.hh"

#include <cmath>
#include <stdexcept>

bool FloatUtils::is_equal(double a, double b, double precision) {
	if (precision < 0.0) {
		throw std::invalid_argument("the precision can't be negative");
	}

	if (std::abs(a - b) < precision) {
		return true;
	}
	
	return false;
}
