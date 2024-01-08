#include "exporter_file.hh"

ExporterFile::ExporterFile(const std::string& filename): filestream(filename) {

}

void ExporterFile::write_point(double t, const Matrix& m) {
	filestream << t << '\t';

	const int rows{m.getRows()};
	const int columns{m.getColumns()};

	for (int i{0}; i < rows; ++i) {
		for (int j{0}; j < columns; ++j) {
			filestream << m[i][j];
			
			const bool is_last_value{i == rows - 1 && j == columns - 1};
			if (!is_last_value) {
				filestream << '\t';
			}
		}
	}

	filestream << '\n';
}
