//
// Created by james on 3/29/20.
//

#ifndef PROJECT_VMATH_H
#define PROJECT_VMATH_H

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <iterator>

#include "nmmintrin.h"

using std::vector;

namespace jaa::vmath {
    //// vector subtraction
    // subtraction in place v -= w == v = v - w
    template <class T>
    void operator-= (vector<T>& v, const vector<T>& w) {
        for (int i = 0; i < v.size(); i++) {
            v[i] -= w[i];
        }
    }

    // subtraction
    template <class T>
    vector<T>& operator- (vector<T>& v, vector<T>& w) {
        auto *out = new vector<T>();

        for (int i = 0; i < v.size(); i++) {
            out->push_back(v[i] - w[i]);
        }

        return *out;
    }


    //// Hadamard (element-wise) product
    // in-place
    template <class T>
    vector<T>& operator*= (vector<T>& v, const vector<T>& w) {
        for (int i = 0; i < v.size(); i++) {
            v[i] *= w[i];
        }
    }

    template <class T>
    vector<T>& operator* (vector<T>& v, vector<T>& w) {
        auto *out = new vector<T>();

        for (int i = 0; i < v.size(); i++) {
            out->push_back(v[i] * w[i]);
        }

        return *out;
    }

    // dot (inner) product
    template <class T>
    T dot(const vector<T>& v, const vector<T>& w) {
        return std::inner_product(v.begin(), v.end(), w.begin(), 0);
    }

    // sum of a vector
    template <class T>
    T sum(const vector<T>& v) {
        return std::accumulate(v.begin(), v.end(), 0);
    }

    // a vector containing the absolute value of every element of v
    template <class T>
    vector<T>& v_abs(vector<T> &v) {
        auto *v_copy = new vector<T>();
        std::copy(v.begin(), v.end(), back_inserter(*v_copy));

        for_each(v_copy->begin(), v_copy->end(), [](T &n){return (n < 0? -n : n);});
        return *v_copy;
    }

    // L2 norm of a vector
    template <class T>
    double norm(vector<T>& v) {
        vector<T> v_copy(v.size());
        std::copy(v.begin(), v.end(), back_inserter(v_copy));

        return std::sqrt(sum(v*v));
    }

    // arithmetic mean of a vector
    template <class T>
    double mean(const vector<T>& v) {
        return (double)sum(v) / v.size();
    }

    // max of a vector
    // This is just a wrapper for some code that C++ makes unnecessarily ugly
    template <class T>
    T max(const vector<T>& v) {
        auto argmax_iterator = std::max_element(v.begin(), v.end());
        auto result_index = std::distance(v.begin(), argmax_iterator);
        return v.at(result_index);
    }

    // de-meaned vector
    // Returns a new vector where each vi = vi - mean(v)
    template <class T>
    vector<double>& demean(const vector<T>& v) {
        vector<double> *v_copy;
        std::copy(v.begin(), v.end(), std::back_inserter(*v_copy));

        double vbar = mean(v);
        std::for_each(v_copy->begin(), v_copy->end(), [vbar](double &x) {
            x = x - vbar;
        });

        return *v_copy;
    }

    // for printing vector contents
    template <class T>
    std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s, " "));
        return s;
    }
}

#endif //PROJECT_VMATH_H
