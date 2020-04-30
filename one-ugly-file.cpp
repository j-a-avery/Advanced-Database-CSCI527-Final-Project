//
// Created by james on 4/10/20.
//

#include <vector>
#include <exception>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>

using namespace std;


// vector subtraction
template <class T>
vector<T> operator-(vector<T> v, vector<T> w) {
    if (v.size() != w.size()) {
        throw exception();
    }

    vector<T> result(v.size());
    for (int i = 0; i < v.size(); i++) {
        result.push_back(v[i] - w[i]);
    }

    return result;
}

// Hadamard product
template <class T>
vector<T> operator*(vector<T> v, vector<T> w) {
    if (v.size() != w.size()) {
        throw exception();
    }

    vector<T> result;
    for (int i = 0; i < v.size(); i++) {
        result.push_back(v[i] * w[i]);
    }

    return result;
}

// Hadamard product with assignment
template <class T>
void operator*=(vector<T> v, vector<T> w) {
    if (v.size() != w.size()) {
        throw exception();
    }

    for (int i = 0; i < v.size(); i++) {
        v[i] *= w[i];
    }
}

// dot (inner) product of two vectors
template <class T>
T dot(vector<T> v, vector<T> w) {
    return inner_product(v.begin(), v.end(), w.begin(), 0);
}

// sum of a vector
template <class T>
T sum(vector<T> v) {
    return accumulate(v.begin(), v.end(), 0);
}

// returns a vector containing the absolute value of every element of v
template <class T>
vector<T> v_abs(const vector<T> &v) {
    vector<T> v_copy(v.size());
    copy(v.begin(), v.end(), back_inserter(v_copy));
    for_each(v_copy.begin(), v_copy.end(), [](T &x) { return abs(x); });
    return v_copy;
}

// L2 norm of a vector
template <class T>
double norm(const vector<T> &v) {
    // Copy v into a new vector v'
    vector<T> v_copy;
    copy(v.begin(), v.end(), back_inserter(v_copy));

    // Square each element of v'
    for_each(v_copy.begin(), v_copy.end(), [](T &x){return x*x;});
    return sqrt(sum(v_copy));
}

// arithmetic mean of a vector
double mean(const vector<int>& v) {
    return (double)(sum(v)) / (double)(v.size());
}

// max of a vector
// It's just a wrapper for some really ugly-looking iterator-based code
int max(vector<int> v) {
    auto argmax_index = max_element(v.begin(), v.end());
    auto result_index = distance(v.begin(), argmax_index);
    return v.at(result_index);
}

// returns a copy v as v - mean(v)
vector<double> demean(vector<int> v) {
    vector<double> v_copy;
    copy(v.begin(), v.end(), back_inserter(v_copy));

    // Subtract mean(v) from each element of v'
    double vbar = mean(v);
    for_each(v_copy.begin(), v_copy.end(), [vbar](double &n){ n = n - vbar; });

    return v_copy;
}

template <class T>
double Eucildean_distance(const vector<T> &v, const vector<T> &w) {
    return norm(v - w);
}

template <class T>
double Manhattan_distance(const vector<T> &v, const vector<T> &w) {
    return (double)(sum(v_abs(v - w)));
}

template <class T>
double Chebyshev_distance(const vector<T> &v, const vector<T> &w) {
    return max(v_abs(v - w));
}

template <class T>
double cosine_similarity(const vector<T> &v, const vector<T> &w) {
    return dot(v, w) / (norm(v) * norm(w));
}

template <class T>
double Pearson_r(const vector<T> &v, const vector<T> &w) {
    return cosine_similarity<double>(demean(v), demean(w));
}

template <class T>
double real_valued_Jaccard_similarity(const vector<T> &v, const vector<T> &w) {
    return dot(v, w) / (dot(v, v) + dot(w, w) - dot(v, w));
}