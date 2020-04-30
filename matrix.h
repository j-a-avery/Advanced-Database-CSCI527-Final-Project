//
// Created by James on 4/14/20.
// A lightweight 2D C-style array wrapped up with size data
// enhanced with row-wise linear algebra operations
//

#ifndef PROJECT_MATRIX_H
#define PROJECT_MATRIX_H

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <initializer_list>

enum DIM {
    ROW,
    COL,
    WHOLE
};

template <typename Number>
class matrix {
private:
    Number** _data;
    int _numrows;
    int _numcols;

public:
    // Constructors & destructor
    matrix();       // An empty 0x0 placeholder matrix
    matrix(int rows, int cols);

    // Matrix creation functions
    static matrix<Number>& from_csv(std::string filename);
    static matrix<Number>& filled_with(int rows, int cols, Number value);
    static matrix<Number>& zero(int rows, int cols);
    static matrix<Number>& one(int rows, int cols);
    matrix<Number>& copy();
    matrix<double>& as_double();
    matrix<int>& as_int();


    // Index operators
    Number& at(int row, int col);
    matrix<Number>& row(int row_index); // TODO: Make this refer to the row, rather than copy it
    Number& value();    // Accesses the only value of a 1x1 matrix


    // Matrix characteristics
    inline int rows() const { return _numrows; }
    inline int cols() const { return _numcols; }


    // Linear algebra and element-wise operations
    matrix<Number>& T();                        // transpose
    matrix<Number>& dot(matrix<Number>);        // dot product: A.dot(B) = sum(a*b)
    matrix<Number>& Hadamard(matrix<Number>&);  // Hadamard product
    matrix<Number>& norm(DIM axis=DIM::WHOLE);  // Frobenius/Euclidean norm
    matrix<double>& tilde();                    // x.tilde() = x - x.mean()
    matrix<Number>& operator-(Number);          // scalar subtraction
    matrix<Number>& operator-(matrix<Number>&); // element-wise subtraction
    matrix<Number>& operator/(Number);          // scalar division
    matrix<Number>& operator/=(Number);         // in-place scalar division
    bool operator==(matrix<Number>&);                    // matrix equality
    bool operator!=(matrix<Number>&);
    // The following functions apply the operation to each element of the matrix
    matrix<Number>& abs();                      // absolute value
    matrix<Number>& square();                   // square each element
    matrix<Number>& sqrt();                     // square root of each element
    matrix<Number>& pow(double exponent);       // raise each element to a power


    // Summary operations
    matrix<Number>& sum(DIM axis=DIM::ROW);
    matrix<Number>& min(DIM axis=DIM::ROW);
    matrix<Number>& max(DIM axis=DIM::ROW);
    matrix<int>& argmax(DIM axis=DIM::ROW);
    matrix<int>& argmin(DIM axis=DIM::ROW);
    matrix<double>& mean(DIM axis=DIM::ROW);


    // output matrices
    template <class T>
    friend std::ostream& operator<<(std::ostream&, matrix<T>&);
};

template class matrix<int>;
template class matrix<double>;

template <typename T>
std::ostream& operator<<(std::ostream& os, matrix<T>& m) {
    std::stringstream s;
    s << m._numrows;
    int label_width = s.str().length();
    for (int row = 0; row < m._numrows; row++) {
        for (int col = 0; col < m._numcols; col++) {
            os << std::setprecision(2) << std::fixed << std::setw(8)
               << m._data[row][col];
        }
        os << std::endl;
    }

    os << "[" << m._numrows << " x " << m._numcols << "] matrix" << std::endl;

    return os;
}


#endif //PROJECT_MATRIX_H
