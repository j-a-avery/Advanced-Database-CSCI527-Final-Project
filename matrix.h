#pragma once
#ifndef MATRIX
#define MATRIX

#include <string>
#include <iostream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <initializer_list>


using
std::string,
std::ostream,
std::stringstream,
std::istringstream,
std::ifstream,
std::cerr,
std::cout,
std::flush,
std::ws,
std::fixed,
std::setw,
std::setprecision,
std::endl;


class matrix_dimension {
	int rows, cols;

public:
	friend std::ostream& operator<<(std::ostream& os, matrix_dimension& d);
	matrix_dimension(int rows, int cols) {
		this->rows = rows;
		this->cols = cols;
	}
};

std::ostream& operator<<(std::ostream& os, matrix_dimension& d) {
	os << '[' << d.rows << ", " << d.cols << ']';
	return os;
}


class matrix {
private:
	int _rows;
	int _cols;
	double** _data;
	int _precision = 3; // Number of decimal places to display

public:
	matrix() {
		_rows = 0;
		_cols = 0;
		_data = nullptr;
	}

	matrix(int rows, int cols) {
		_rows = rows;
		_cols = cols;
		_data = new double*[rows];

		for (int row = 0; row < rows; row++) {
			_data[row] = new double[cols];
		}
	};

	matrix(int rows, int cols, std::initializer_list<double> data) {
		_rows = rows;
		_cols = cols;

		_data = new double*[rows];

		int row = 0;
		int col = 0;

		for (double datum : data) {
			if (col == 0) {
				_data[row] = new double[cols];
			}

			_data[row][col] = datum;

			col++;
			if (col == cols) {
				col = 0;
				row++;
			}
		}
	};

	~matrix() {
		for (int row = 0; row < _rows; row++) {
			delete[] _data[row];
		}

		delete[] _data;
	};

	inline double& at(int row, int col) {
		return _data[row][col];
	};

	inline int rows() const {
		return _rows;
	};

	inline int cols() const {
		return _cols;
	};

	inline matrix_dimension& shape() const {
		auto md = new matrix_dimension(_rows, _cols);
		return *md;
	};

	inline matrix& setprecision(int p) {
		_precision = p;
		return *this;
	}

	inline int precision() const {
		return _precision;
	}

	// Equality operator
	// matrix A == matrix B iff the dimensions and all elements are the same
	bool operator==(matrix& that) {
		if (this->rows() != that.rows()) {
			return false;
		}

		if (this->cols() != that.cols()) {
			return false;
		}

		for (int row = 0; row < this->rows(); row++) {
			for (int col = 0; col < this->cols(); col++) {
				if (this->at(row, col) != that.at(row, col)) {
					return false;
				}
			}
		}

		return true;
	}

	bool operator!=(matrix& that) {
		return !(*this == that);
	}
};	


ostream& operator<<(ostream& os, matrix& m) {
	double max = std::numeric_limits<double>::min();
	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			if (m.at(row, col) > max) {
				max = m.at(row, col);
			}
		}
	}

	std::stringstream s;
	s << (int)max;
	int width = (int)s.str().length() + 4 + m.precision();

	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			os << fixed << setw(width) << setprecision(m.precision()) << m.at(row, col);
		}

		os << endl;
	}

	os << "\t[" << m.rows() << " x " << m.cols() << "] matrix" << endl;

	return os;
}


matrix& load_csv(string path) {
    // First, determine the shape of the data in the file
    int numrows = 0;
    int numcols = 0;
    int maxcols = 0;

    // Try to open the file
    ifstream datafile(path);

    if (!datafile.is_open()) {
        cerr << "There was a problem opening \"" << path << "\"" << endl;
        return *new matrix();
    }

    string line;        // The contents of each line in the file
    double next_item;   // The next item to be read from the file

    cout << "Reading data file " << path << std::endl;

    while (getline(datafile, line)) {
        numrows += 1;
        numcols = 0;

        istringstream buffer{ line };

		while (buffer >> next_item) {
			numcols += 1;

			if ((buffer >> ws).peek() == ',') {
				buffer.ignore();
			}

			if (numcols > maxcols) {
				maxcols = numcols;
			}
		}
    }

    cout << "Read " << numrows << " rows x " << numcols << " columns." << endl;

    // Create a matrix with the file's contents
    cout << "Building matrix. " << flush;
    matrix* M = new matrix(numrows, maxcols);

    // Return to the beginning of the file
    datafile.clear();
    datafile.seekg(0);
    int row = 0;

    // Fead the file's contents into the new matrix
    while (getline(datafile, line)) {
        int col = 0;
        istringstream buffer{ line };

        while (buffer >> next_item) {
            M->at(row, col) = next_item;
            col += 1;

            if ((buffer >> ws).peek() == ',') {
                buffer.ignore();
            }
        }

        row += 1;

        if (row % (numrows / 10) == 0) {
            cout << ". ";
        }
    }

    cout << "Done!" << endl;

    // return the matrix
    return *M;
}

#endif