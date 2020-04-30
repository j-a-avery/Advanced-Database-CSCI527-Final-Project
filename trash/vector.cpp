//
// Created by james on 4/11/20.
//

#include <cmath>
#include "vector.h"

//#include <iostream> // for debugging. TODO: DELETE THIS

using namespace jaa;


//////////////////// Vector Construction Methods //////////////////////////////

// Default constructor
// Constructs an empty vector
// By default, the underlying array holds 10 elements
// and doubles in size whenever it runs out of space
template <class T>
vector<T>::vector() {
    _capacity = 10;
    _contents = new T[_capacity];
    _size = 0;
}


// Alternate constructor
// Constructs the underlying array to a specified dimension
template <class T>
vector<T>::vector(int size) {
    _capacity = size;
    _contents = new T[size];
    _size = 0;
}


// Initializer list constructor
// Creates a vector filled with the items from the initializer_list
template <class T>
vector<T>::vector(std::initializer_list<T> that) {
    // Add a little wiggle room to avoid resizing right away
    this->_capacity = that.size() * 1.5;
    this->_contents = new T[this->_capacity];

    int idx = 0;
    for (auto i : that) {
        this->_contents[idx++] = i;
    }

    _size = idx;
}


// Copy constructor
template <class T>
vector<T>::vector(vector<T> &that) {
    this->_capacity = that._capacity;
    this->_size = that._size;

    this->_contents = new T[this->_capacity];

    for (int i = 0; i < this->_size; i++) {
        this->_contents[i] = that._contents[i];
    }
}


// Destructor
template <class T>
vector<T>::~vector() {
    delete[] this->_contents;
}


// Provides a copy of the current vector
// Useful for math operations when in-place operation is undesirable
template <class T>
vector<T> vector<T>::copy() {
    vector<T> newvec{*this};
    return newvec;
}


// Create a vector of a certain size, filled with the specified value
template <class T>
vector<T> vector<T>::filled_with(T value, int size) {
    vector<T> newvec(size);
    for (int i = 0; i < size; i++) {
        newvec.push_back(value);
    }

    return newvec;
}


// Get the size of the vector
template <class T> inline
int vector<T>::size() const {
    return this->_size;
}


// Append an item to the end of the vector
template <class T>
void vector<T>::push_back(T item) {
    // If the underlying array is full
    if (_size == _capacity) {
        // then increase its capacity...
        int expanded_capacity = 1.5 * _size;
        T* expanded_contents = new T[expanded_capacity];

        // and copy its contents to the new array
        for (int i = 0; i < _size; i++) {
            expanded_contents[i] = _contents[i];
        }

        _capacity = expanded_capacity;
        _contents = expanded_contents;
    }

    // Add the new item to the back of the array
    _contents[_size++] = item;
}


// Append several items to the end of the vector
template <class T>
void vector<T>::push_back(std::initializer_list<T> items) {
    // If the underlying array is full,
    // or if it is about to become full...
    if (_size + items.size() >= _capacity) {
        // then increase its capacity...
        int expanded_capacity = _capacity + (int)(1.5 * items.size());
        T* expanded_contents = new T[expanded_capacity];

        // and copy its contents to the new array
        for (int i = 0; i < _size; i++) {
            expanded_contents[i] = _contents[i];
        }

        _capacity = expanded_capacity;
        _contents = expanded_contents;
    }

    // Add the new item to the back of the array
    for (auto i : items) {
        _contents[_size++] = i;
    }
}


//////////////////// Cast Methods /////////////////////////////////////////////
// In a library situation there would be a lot more of these,
// but for this project I only need one

// Cast vector<T> to vector<double>
template <class T>
vector<double> vector<T>::as_double() {
    vector<double> newvec(_size);

    for (int i = 0; i < _size; i++) {
        newvec.push_back((double)_contents[i]);
    }

    return newvec;
}



/////////////////// Vector Summary Operations /////////////////////////////////

// Sum of the items in the vector
template <class T>
T vector<T>::sum() {
    T _sum = 0;

    for (int i = 0; i < _size; i++) {
        _sum += _contents[i];
    }

    return _sum;
}

// Arithmetic mean of the items in the vector
template <class T> inline
double vector<T>::mean() {
    return (double)this->sum() / (double)this->size();
}

// Largest item in the vector
template <class T>
T vector<T>::max() {
    T _max = _contents[0];

    for (int i = 0; i < _size; i++) {
        if (_contents[i] > _max) {
            _max = _contents[i];
        }
    }

    return _max;
}

// Index of the largest item in the vector
template <class T>
int vector<T>::argmax() {
    int _max = 0;

    for (int i = 0; i < _size; i++) {
        if (_contents[i] > _contents[_max]) {
            _max = i;
        }
    }

    return _max;
}

// element-wise absolute value of the vector
template <class T>
vector<T> vector<T>::abs() {
    vector<T> that = this->copy();

    for (int i = 0; i < _size; i++) {
        that[i] = std::abs(that[i]);
    }

    return that;
}

// L2 norm of the vector
template <class T>
double vector<T>::norm() {
    return std::sqrt( (double)(*this).dot(*this) );
}

// "de-mean" the vector, that is, produce a new vector whose contents are
// the previous vector minus the mean
// This is a very useful concept from Boyd * Vandenberghe (2018)
// Introduction to Applied Linear Algebra: Vectors, Matrices, and Least Squares
// for xi ∈ x, x̃i = xi - x̄
template <class T>
vector<double> vector<T>::demean() {
    vector<double> x_bar = vector<double>::filled_with(this->mean(), this->size());
    //std::cout << "debug: x̄=" << x_bar << std::endl;
    return (this->as_double() - x_bar);
}

// dot (inner) product of two vectors
// A.B = sum( [a1b1, a2b2, ... anbn] ) = sum( A * B ), where * is the Hadamard product
template <class T>
T vector<T>::dot(vector<T>& that) {
    // TODO: Bounds checking
    return ( (*this) * that ).sum();
}


/////////////////// Vector Arithmetic Operations //////////////////////////////

template <class T>
T& vector<T>::operator[](int idx) const {
    // TODO: Bounds checking
    if (idx > _size) idx = _size;

    return _contents[idx];
}

// vector subtraction
// A - B = [a1-b1, a2-b2, ..., an-bn]
template <class T>
vector<T> vector<T>::operator-(vector<T> that) {
    // TODO: Bounds checking
    vector<T> newvec = this->copy();

    for (int i = 0; i < _size; i++) {
        newvec[i] -= that[i];
    }

    return newvec;
}

// Hadamard product
// In a different setting, there would probably be a better choice for the *
// operator than the Hadamard product, but since it comes up so much in this
// project, and other forms of vector multiplication do not (aside from the
// inner product, which uses the Hadamard product), * is defined as
// A * B = [a1*b1, a2*b2, ..., an*bn]
template <class T>
vector<T> vector<T>::operator*(vector<T> that) {
    // TODO: Bounds checking
    vector<T> newvec = this->copy();

    for (int i = 0; i < _size; i++) {
        newvec[i] *= that[i];
    }

    return newvec;
}

// Scalar multiplication
// a*V = [a*v1, a*v2, ..., a*vn]
template <class T>
vector<T> vector<T>::operator*(T scalar) {
    // TODO: bounds checking
    vector<T> newvec = this->copy();

    for (int i = 0; i < _size; i++) {
        newvec[i] *= scalar;
    }

    return newvec;
}

template <class T>
vector<T> operator*(T scalar, vector<T> v) {
    vector<T> v_copy = v.copy();
    return v_copy * scalar;
}