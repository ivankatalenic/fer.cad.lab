#pragma once

#include "matrix.hh"

#include <string>
#include <fstream>

class ExporterFile {
public:
	ExporterFile(const std::string& filename);

	void write_point(double t, const Matrix& m);
private:
	std::ofstream filestream;
};
