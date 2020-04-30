//
// Created by james on 4/11/20.
//

#ifndef PROJECT_VECTOR_H
#define PROJECT_VECTOR_H


#include <ostream>
#include <initializer_list>

namespace jaa {
    template <class T>
    class vector {
    private:
        T* _contents;   // underlying array
        int _size{};      // number of elements in the array
        int _capacity{};  // number of elements that can be stored
                        // without resizing

    public:
        vector<T>();
        explicit vector<T>(int size);
        vector<T>(std::initializer_list<T> that);
        vector<T>(vector<T> &that);
        ~vector<T>();
        vector<T> copy();

        static vector<T> filled_with(T value, int size);

        vector<double> as_double(); // recast to double

        int size() const;
        void push_back(T item);
        void push_back(std::initializer_list<T> items);

        T sum();
        double mean();
        T max();
        int argmax();
        vector<T> abs();
        double norm(); // L2 norm of the vector
        vector<double> demean();

        T dot(vector<T>& that);

        T& operator[](int idx) const;
        vector<T> operator-(vector<T> that);
        vector<T> operator*(vector<T> that);
        vector<T> operator*(T scalar);

        template <class U>
        friend vector<U> operator*(U scalar, vector<U> that);
    };

    template class vector<int>;
    template class vector<double>;

    // Stream insertion operator for printing vectors
    template <class T> inline
    std::ostream& operator<<(std::ostream& str, const jaa::vector<T> &v) {
        str << "( ";

        int idx = 0;
        for (; idx < v.size() - 1; ) {
            str << v[idx++] << ", ";
        }
        str << v[idx];

        str << " )";
        return str;
    }
}



#endif //PROJECT_VECTOR_H
