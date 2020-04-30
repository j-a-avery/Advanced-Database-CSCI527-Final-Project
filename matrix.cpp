//
// Created by james on 4/14/20.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>

#include "matrix.h"


/**********************************************************
 *
 *          Constructors, Destructors, and Matrix Creators
 *
 **********************************************************/

/**********************************************************
 *  Default Constructor
 *
 *  @tparam Number The type of the data. Currently limited to int and double
 *  @return An empty 0-dimensional matrix
 *
 *  This constructor should be avoided; It is included just for completeness.
 **********************************************************/
template <typename Number>
matrix<Number>::matrix() {
    _numrows = 0;
    _numcols = 0;
    _data = new Number*[0];
}


/**********************************************************
 *  Preferred Constructor
 *  Constructs a `rows` x `cols` matrix
 *
 *  @tparam Number The type of the data. Currently limited to int and double
 *  @param cols Width of the matrix
 *  @param rows Height of the matrix
 *  @return A new rows x cols matrix object
 **********************************************************/
template <typename Number>
matrix<Number>::matrix(int rows, int cols) {
    _numrows = rows;
    _numcols = cols;

    _data = new Number*[_numrows];
    for (int i = 0; i < _numrows; i++) {
        _data[i] = new Number[_numcols];
    }
}


/**********************************************************
 * Creates a matrix filled with numeric data from a CSV file
 *
 * @tparam Number The type of the data. Currently limited to int and double
 * @param filename Path to the file
 * @return A reference to a matrix containing the file's contents
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::from_csv(std::string filename) {
    // First, determine the shape of the data in the file
    // It is assumed that the data in the file may be jagged,
    // but the resulting matrix will not be
    int numrows = 0;
    int numcols = 0;
    int maxcols = 0;

    // Try to open the file.
    std::ifstream datafile(filename);

    if (!datafile.is_open()) {
        // Todo: Make this throw an exception instead of returning an empty matrix
        std::cerr << "There was a problem opening " << filename << "." << std::endl;
        return *new matrix<Number>();
    }

    std::string line;   // The contents of each line in the file
    Number next_item;   // The next item to be read from the file

    std::cout << "Reading data file \"" << filename << "\"" << std::endl;

    while (std::getline(datafile, line)) {
        numrows += 1;
        numcols = 0;
        std::istringstream buffer{line};

        while (buffer >> next_item) {
            numcols += 1;

            if ((buffer >> std::ws).peek() == ',') {
                buffer.ignore();
            }

            if (numcols > maxcols) {
                maxcols = numcols;
            }
        }

    }

    std::cout << "Read " << numrows << " rows × " << numcols << " columns." << std::endl;

    // Create a matrix with the file's contents
    std::cout << "Building matrix. ";
    datafile.clear(); // Return to the beginning of the file
    datafile.seekg(0);
    matrix<Number>* M = new matrix<Number>(numrows, numcols); // Create the matrix object
    int row = 0;

    // Read the file's contents into the new matrix
    while (std::getline(datafile, line)) {
        int col = 0;
        std::istringstream buffer{line};

        while (buffer >> next_item) {
            M->at(row, col) = next_item;
            col += 1;

            if ((buffer >> std::ws).peek() == ',') {
                buffer.ignore();
            }
        }
        row += 1;

        if (row % (numrows / 10) == 0) {
            std::cout << ". ";
        }
    }
    std::cout << "Done!" << std::endl;


    // Return the matrix
    return *M;
}


/**********************************************************
 * Creates a matrix filled with a specified value
 *
 * @tparam Number The type of the data in the matrix
 * @param rows Height of the matrix
 * @param cols Width of the matrix
 * @param value Value to fill the matrix with
 * @return A new rows x cols matrix with each element set to value
 **********************************************************/
 template <typename Number>
 matrix<Number>& matrix<Number>::filled_with(int rows, int cols, Number value) {
     matrix<Number>* M = new matrix<Number>(rows, cols);

     for (int i = 0; i < rows; i++) {
         for (int j = 0; j < rows; j++) {
             M->at(i, j) = value;
         }
     }

     return *M;
 }


/**********************************************************
 * Creates a zero matrix
 *
 * @tparam Number The type of the data. Currently limited to int and double
 * @param rows Height of the matrix
 * @param cols Width of the matrix
 * @return A new rows x cols matrix with each element set to 0
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::zero(int rows, int cols) {
    return matrix<Number>::filled_with(rows, cols, (Number)0);
}


/**********************************************************
 * Creates a ones matrix
 *
 * @tparam Number The type of the data. Currently limited to int and double
 * @param rows Height of the matrix
 * @param cols Width of the matrix
 * @return A new rows x cols matrix with each element set to 1
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::one(int rows, int cols) {
    return matrix<Number>::filled_with(rows, cols, (Number)1);
}


/**********************************************************
 * Creates a copy of the matrix
 *
 * @tparam Number
 * @return A copy of the matrix
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::copy() {
    matrix<Number>* _copy = new matrix<Number>(this->rows(), this->cols());

    for (int i = 0; i < this->rows(); i++) {
        for (int j = 0; j < this->cols(); j++) {
            _copy->at(i, j) = this->at(i, j);
        }
    }

    return *_copy;
}


/**********************************************************
 * Casts the matrix to double
 *
 * @tparam Number
 * @return A copy of the matrix with each element cast to double
 **********************************************************/
template <typename Number>
matrix<double>& matrix<Number>::as_double() {
    matrix<double>* _copy = new matrix<double>(this->rows(), this->cols());

    for (int i = 0; i < this->rows(); i++) {
        for (int j = 0; j < this->cols(); j++) {
            _copy->at(i, j) = (double)(this->at(i, j));
        }
    }

    return *_copy;
}


/**********************************************************
 * Casts the matrix to int
 *
 * @tparam Number
 * @return A copy of the matrix with each element cast to int
 **********************************************************/
template <typename Number>
matrix<int>& matrix<Number>::as_int() {
    matrix<int>* _copy = new matrix<int>(this->rows(), this->cols());

    for (int i = 0; i < this->rows(); i++) {
        for (int j = 0; j < this->cols(); j++) {
            _copy->at(i, j) = (int)(this->at(i, j));
        }
    }

    return *_copy;
}



/**********************************************************
 *
 *          Indexing functions
 *
 **********************************************************/


template <typename Number>
Number& matrix<Number>::at(int row, int col) {
    return _data[row][col];
}


/**********************************************************
 *
 * @return A new matrix whose _data points to the _data of the
 *         row_index-th row of the parent matrix
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::row(int row_index) {
    matrix<Number>* r = new matrix<Number>(1, this->cols());
    r->_data[0] = this->_data[row_index];

    return *r;
}


/**********************************************************
 *
 * @tparam Number
 * @return The only value of a 1x1 matrix
 */
template <typename Number>
Number& matrix<Number>::value() {
    if (this->rows() != 1 || this->cols() != 1) {
        throw std::exception(); // Todo: Make this specific
    }

    return this->at(0, 0);
}



/**********************************************************
 *
 *          Linear Algebra Operations
 *
 **********************************************************/


/**********************************************************
 * @return the transpose of the matrix
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::T() {
    matrix<Number>* mT = new matrix<Number>(this->cols(), this->rows());

    for (int row = 0; row < this->rows(); row++) {
        for (int col = 0; col < this->cols(); col++) {
            mT->at(col, row) = this->at(row, col);
        }
    }

    return *mT;
}


/**********************************************************
 * Dot product
 *
 * @param that A matrix
 * @return A reference to a new matrix containing the dot product
 *         of this . that
 *
 * The two matrices should be compatible, i.e.
 * A in F[m x n] and B in F[p x q] then it should be the case that n == p
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::dot(matrix<Number> that) {
    // Ensure that the matrices are compatible.
    // If they are not, then return an empty matrix
    // todo: Add more robust dimension checking. This is a terrible way to do this
    if (this->cols() != that.rows()) {
        return *(new matrix<Number>());
    }

    // A dot B = sum(ab) for a in A and b in B
    matrix<Number>* DP = &matrix<Number>::zero(this->rows(), that.cols());

    for (int i = 0; i < this->rows(); i++) {
        for (int j = 0; j < that.cols(); j++) {
            for (int k = 0; k < this->cols(); k++) {
                DP->at(i, j) += this->at(i, k) * that.at(k, j);
            }
        }
    }

    return *DP;
}


/**********************************************************
 * Hadamard product
 *
 * @param that
 * @return A reference to a new matrix whose elements are the Hadamard product
 *         of the other two products.
 *
 * Hadamard(A, B) = a * b for a in A, for b in B
 * As with subtraction, broadcasts for one-row matrices
 */
template <typename Number>
matrix<Number>& matrix<Number>::Hadamard(matrix<Number>& that) {
    // Ensure that the matrices are compatible
    // Todo: Add better dimension checking
    if (this->cols() != that.cols()) {
        return *(new matrix<Number>());
    }

    if (this->rows() != that.rows() && that.rows() != 1) {
        return *(new matrix<Number>());
    }

    // If `that` has the same dimensions as `this`,
    // then multiply as normal by incrementing its row index.
    // If `that` has only one row,
    // then broadcast it to each row, by not incrementing its row index.
    int that_inc = (that.rows() == 1 ? 0 : 1);

    matrix<Number>* Hprod = new matrix<Number>(this->rows(), this->cols());

    for (int this_row = 0, that_row = 0;
            this_row < this->rows();
            this_row++, that_row += that_inc) {
        for (int col = 0; col < this->cols(); col++) {
            Hprod->at(this_row, col) = this->at(this_row, col) * that.at(that_row, col);
        }
    }

    return *Hprod;
}


/**********************************************************
 * Matrix norm
 *
 * For a 1xn or nx1 matrix, this is equivalent to the Euclidean norm.
 * For an n x m matrix, this is the Frobenius norm
 * (i.e. the Euclidean norm generalized for multiple dimensions).
 * A.norm() = sqrt(sum(a*a)) for a in A
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::norm(DIM axis) {
    // The clearest, most concise and readable way to write this would be a one-liner,
    // return this->square().sum(axis).sqrt();
    // but that caused a really low-level malloc/free error,
    // so rather than figure out how to fix it, I'm taking the easy road
    // and just rewriting the function to calculate the norm manually
    matrix<Number>* norms = new matrix();
    Number sum; // This should be a local variable in each case,
                // but one of the fundamental design philosophies of C++
                // seems to be to make everything as
                // obtuse, inconsistent, and complicated as possible.

    switch (axis) {
        case DIM::ROW:
            // Calculate a matrix whose rows are the L2 norms of each row vector
            norms = new matrix<Number>(this->rows(), 1);

            for (int row = 0; row < this->rows(); row++) {
                sum = 0;
                for (int col = 0; col < this->cols(); col++) {
                    sum += std::pow(this->at(row, col), 2);
                }
                norms->at(row, 0) = std::sqrt(sum);
            }
            break;
        case DIM::COL:
            // Calculate a matrix whose columns are the L2 norms of each column vector
            norms = new matrix<Number>(1, this->cols());

            for (int col = 0; col < this->cols(); col++) {
                sum = 0;
                for (int row = 0; row < this->rows(); row++) {
                    sum += std::pow(this->at(row, col), 2);
                }
                norms->at(0, col) = std::sqrt(sum);
            }
            break;
        case DIM::WHOLE:
            // Calculate the Frobenius norm of the matrix
            sum = 0;
            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    sum += std::pow(this->at(row, col), 2);
                }
            }

            sum = std::sqrt(sum);

            norms = new matrix(1, 1);

            norms->at(0, 0) = sum;
            // ...
            break;
        default:
            // This should never happen
            // error
            break;
    }

    return *norms;
}


/**********************************************************
 * Tilde function, named after a non-standard notation used by
 * Boyd & Vandenberghe, 2018, were x̃ = x - x̄,
 * i.e. x.tilde() == x - x.mean()
 *
 * This is useful in statistical functions, Spearman correlation in particular.
 *
 * @tparam Number
 * @return
 **********************************************************/
template <typename Number>
matrix<double>& matrix<Number>::tilde() {
    return (this->copy().as_double() - this->mean(DIM::WHOLE).value());
}


/**********************************************************
 * Scalar subtraction
 *
 * @param scalar The number to subtract from each element
 * @return A copy of the matrix, with the scalar subtracted from each element
 *
 * A - n = a - n for al a in A
 **********************************************************/
template <typename R>
matrix<R>& matrix<R>::operator-(R scalar) {
    matrix<R>* product = &this->copy();

    for (int i = 0; i < product->rows(); i++) {
        for (int j = 0; j < product->cols(); j++) {
            product->at(i, j) -= scalar;
        }
    }

    return *product;
}


/**********************************************************
 * Matrix subtraction, with broadcasting where applicable
 *
 * @param that The matrix to subtract
 * @return A copy of the matrix, with the elements of `that` subtracted from
 *         the elements of `this`
 *
 * This function does broadcast across its rows,
 * so if a 1 x n matrix is subtracted from a m x n matrix,
 * then the row will be treated as an m x n matrix where each row
 * is the same as the 1-row matrix.
 */
template <typename R>
matrix<R>& matrix<R>::operator-(matrix<R>& that) {
    // The matrices must either be the same dimensions,
    // or `that` must be 1-row with the same number of columns as `this`
    // todo: Add proper dimension checking
    if (this->cols() != that.cols()) {
        return *(new matrix<R>(0, 0)); // still a terrible way to do this
    }

    if (that.rows() != this->rows() && that.rows() != 1) {
        return *(new matrix<R>(0, 0));
    }

    // Create a new matrix to be returned
    matrix<R>* product = new matrix<R>(this->rows(), this->cols());

    // If `that` has the same dimensions as `this`,
    // then subtract as normal by incrementing its row index.
    // If `that` has only one row,
    // then broadcast it to each row, by not incrementing its row index.
    int that_inc = (that.rows() == 0 ? 0 : 1);

    // Perform the subtraction
    for (int this_i = 0, that_i = 0; this_i < this->rows(); this_i++, that_i += that_inc) {
        for (int j = 0; j < this->cols(); j++) {
            product->at(this_i, j) = this->at(this_i, j) - that.at(that_i, j);
        }
    }

    return *product;
}


/**********************************************************
 * Scalar division
 *
 * @param scalar
 * @return A copy of the matrix with each element divided by the scalar
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::operator/(Number scalar) {
    matrix* that = &this->copy();

    for (int row = 0; row < that->rows(); row++) {
        for (int col = 0; col < that->cols(); col++) {
            that->at(row, col) /= scalar;
        }
    }

    return *that;
}


/**********************************************************
 * In-place scalar division
 *
 * @param scalar
 * @return the original matrix with each element divided by the scalar
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::operator/=(Number scalar) {
    for (int row = 0; row < this->rows(); row++) {
        for (int col = 0; col < this->cols(); col++) {
            this->at(row, col) /= scalar;
        }
    }

    return *this;
}


/**********************************************************
 * Matrix equality.
 *
 * A == B = true iff a = b for a in A, b in B
 * @tparam Number
 * @param that
 * @return
 */
template <typename Number>
bool matrix<Number>::operator==(matrix<Number>& that) {
    if (this->rows() != that.rows() || this->cols() != that.cols()) {
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

template <typename Number>
bool matrix<Number>::operator!=(matrix<Number>& that) {
    return !(*this == that);
}


/**********************************************************
 * @return a matrix with the absolute value of each element
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::abs() {
    matrix<Number>* that = &this->copy();

    for (int row = 0; row < that->rows(); row++) {
        for (int col = 0; col < that->cols(); col++) {
            that->at(row, col) = std::abs(this->at(row, col));
        }
    }

    return *that;
}


/**********************************************************
 * @return a matrix with each element squared
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::square() {
    matrix<Number>* that = &this->copy();

    for (int row = 0; row < that->rows(); row++) {
        for (int col = 0; col < that->cols(); col++) {
            that->at(row, col) = std::pow(this->at(row, col), 2);
        }
    }

    return *that;
}


/**********************************************************
 * @return a matrix with the square root of each element
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::sqrt() {
    matrix<Number>* that = &this->copy();

    for (int row = 0; row < that->rows(); row++) {
        for (int col = 0; col < that->cols(); col++) {
            that->at(row, col) = std::sqrt(this->at(row, col));
        }
    }

    return *that;
}


/**********************************************************
 * @return a matrix with the absolute value of each element
 **********************************************************/
template <typename Number>
matrix<Number>& matrix<Number>::pow(double exponent) {
    matrix<Number>* that = &this->copy();

    for (int row = 0; row < that->rows(); row++) {
        for (int col = 0; col < that->cols(); col++) {
            that->at(row, col) = std::pow(this->at(row, col), exponent);
        }
    }

    return *that;
}


template <typename R>
matrix<R>& matrix<R>::sum(DIM axis) {
    matrix<R>* sums = nullptr; // Initialize unnecessarily to get Clang-Tidy off my back.

    switch (axis) {
        case DIM::ROW:
            // sum across rows
            sums = &matrix::zero(this->rows(), 1);
            for (int i = 0; i < this->rows(); i++) {
                for (int j = 0; j < this->cols(); j++) {
                    sums->at(i, 0) += this->at(i, j);
                }
            }
            break;
        case DIM::COL:
            // sum along columns
            sums = &matrix::zero(1, this->cols());
            for (int i = 0; i < this->rows(); i++) {
                for (int j = 0; j < this->cols(); j++) {
                    sums->at(0, j) += this->at(i, j);
                }
            }
            break;
        case DIM::WHOLE:
            // sum the entire matrix
            sums = new matrix<R>(1, 1);
            for (int i = 0; i < this->rows(); i++) {
                for (int j = 0; j < this->cols(); j++) {
                    sums->at(0, 0) += this->at(i, j);
                }
            }
            break;
    }

    return *sums;
}


template <typename Number>
matrix<Number>& matrix<Number>::min(DIM axis) {
    matrix<Number>* mins = nullptr;

    switch(axis) {
        case DIM::ROW:
            // find the min element of each row
            mins = new matrix<Number>(this->rows(), 1);

            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    mins->at(row, 0) = std::numeric_limits<Number>::max();
                    if (this->at(row, col) < mins->at(row, 0)) {
                        mins->at(row, 0) = this->at(row, col);
                    }
                }
            }

            break;

        case DIM::COL:
            // find the minimum element of each column
            mins = new matrix<Number>(1, this->cols());

            for (int col = 0; col < this->cols(); col++) {
                for (int row = 0; row < this->rows(); row++) {
                    mins->at(0, col) = std::numeric_limits<Number>::max();
                    if (this->at(row, col) < mins->at(0, col)) {
                        mins->at(0, col) = this->at(row, col);
                    }
                }
            }

            break;

        case DIM::WHOLE:
            // find the maximum element of the whole matrix
            mins = new matrix<Number>(1, 1);
            mins->at(0, 0) = std::numeric_limits<Number>::max();

            for (int col = 0; col < this->cols(); col++) {
                for (int row = 0; row < this->rows(); row++) {
                    if (this->at(row, col) < mins->value()) {
                        mins->at(0, 0) = this->at(row, col);
                    }
                }
            }

            break;
    }

    return *mins;
}


template <typename Number>
matrix<Number>& matrix<Number>::max(DIM axis) {
    matrix<Number>* maxes = nullptr;

    switch(axis) {
        case DIM::ROW:
            // find the max element of each row
            maxes = new matrix<Number>(this->rows(), 1);

            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    maxes->at(row, 0) = std::numeric_limits<Number>::min();

                    if (this->at(row, col) > maxes->at(row, 0)) {
                        maxes->at(row, 0) = this->at(row, col);
                    }
                }
            }

            break;

        case DIM::COL:
            // find the maximum element of each column
            maxes = new matrix<Number>(1, this->cols());

            for (int col = 0; col < this->cols(); col++) {
                for (int row = 0; row < this->rows(); row++) {
                    maxes ->at(0, col) = std::numeric_limits<Number>::min();

                    if (this->at(row, col) > maxes->at(0, col)) {
                        maxes->at(0, col) = this->at(row, col);
                    }
                }
            }

            break;

        case DIM::WHOLE:
            // find the maximum value of the whole array
            maxes = new matrix<Number>(1, 1);
            maxes->at(0, 0) = std::numeric_limits<Number>::min();

            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    if (this->at(row, col) > maxes->value()) {
                        maxes->at(0, 0) = this->at(row, col);
                    }
                }
            }

            break;
    }

    return *maxes;
}


template <typename Number>
matrix<int>& matrix<Number>::argmax(DIM axis) {
    matrix<int>* argmaxes = nullptr;

    switch (axis) {
        case DIM::ROW:
            // Find the column index of the maximum element of each row
            argmaxes = new matrix<int>(this->rows(), 1);

            for (int row = 0; row < this->rows(); row++) {
                Number maxval = std::numeric_limits<Number>::min();  // The current row max is the smallest possible value of Number
                for (int col = 0; col < this->cols(); col++) {
                    if (this->at(row, col) > maxval) {
                        maxval = this->at(row, col);
                        argmaxes->at(row, 0) = col;
                    }
                }
            }

            break;
        case DIM::COL:
            // Find the row index of the maximum element of each column
            argmaxes = new matrix<int>(1, this->cols());

            for (int col = 0; col < this->cols(); col++) {
                Number maxval = std::numeric_limits<Number>::min(); // The current col max is the smallest possible value of Number
                for (int row = 0; row < this->rows(); row++) {
                    if (this->at(row, col) > maxval) {
                        maxval = this->at(row, col);
                        argmaxes->at(0, col) = row;
                    }
                }
            }

            break;
        case DIM::WHOLE:
            // Find the index of the maximum element of the entire matrix
            argmaxes = new matrix<int>(1, 2);

            Number maxval = std::numeric_limits<Number>::min(); // The smallest possible value of Number

            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    if (this->at(row, col) > maxval) {
                        maxval = this->at(row, col);
                        argmaxes->at(0, 0) = row;
                        argmaxes->at(0, 1) = col;
                    }
                }
            }
            break;
    }

    return *argmaxes;
}


template <typename Number>
matrix<int>& matrix<Number>::argmin(DIM axis) {
    matrix<int>* argmins = nullptr;

    switch (axis) {
        case DIM::ROW:
            // Find the column index of the minimum element of each row
            argmins = new matrix<int>(this->rows(), 1);

            for (int row = 0; row < this->rows(); row++) {
                Number minval = std::numeric_limits<Number>::max();  // The current row min is the largest possible value of Number
                for (int col = 0; col < this->cols(); col++) {
                    if (this->at(row, col) < minval) {
                        minval = this->at(row, col);
                        argmins->at(row, 0) = col;
                    }
                }
            }

            break;
        case DIM::COL:
            // Find the row index of the minimum element of each column
            argmins = new matrix<int>(1, this->cols());

            for (int col = 0; col < this->cols(); col++) {
                Number minval = std::numeric_limits<Number>::max(); // The current col min is the largest possible value of Number
                for (int row = 0; row < this->rows(); row++) {
                    if (this->at(row, col) < minval) {
                        minval = this->at(row, col);
                        argmins->at(0, col) = row;
                    }
                }
            }

            break;
        case DIM::WHOLE:
            // Find the index of the minimum element of the entire matrix
            argmins = new matrix<int>(1, 2);

            Number minval = std::numeric_limits<Number>::max(); // The largest possible value of Number

            for (int row = 0; row < this->rows(); row++) {
                for (int col = 0; col < this->cols(); col++) {
                    if (this->at(row, col) < minval) {
                        minval = this->at(row, col);
                        argmins->at(0, 0) = row;
                        argmins->at(0, 1) = col;
                    }
                }
            }
            break;
    }

    return *argmins;
}


template <typename Number>
matrix<double>& matrix<Number>::mean(DIM axis) {
    matrix<double>* _mean = &this->as_double().sum(axis);

    switch (axis) {
        case DIM::ROW:
            *_mean /= this->cols();
            break;
        case DIM::COL:
            *_mean /= this->rows();
            break;
        case DIM::WHOLE:
            *_mean /= this->rows() * this->cols();
            break;
    }

    return *_mean;
}