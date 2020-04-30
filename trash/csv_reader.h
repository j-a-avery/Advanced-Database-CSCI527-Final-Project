//
// Created by james on 3/28/20.
//

#ifndef PROJECT_CSV_READER_H
#define PROJECT_CSV_READER_H

#include <fstream>
#include "vector.h"
#include <string>

using jaa::vector, std::ifstream, std::string;

template <class T>
class csv_reader {
private:
    vector<vector<T>> data;

public:
    explicit csv_reader(const string&);
    vector<vector<T>>& get_data();
};

template class csv_reader<int>;
template class csv_reader<double>;
template class csv_reader<long long>;


#endif //PROJECT_CSV_READER_H
