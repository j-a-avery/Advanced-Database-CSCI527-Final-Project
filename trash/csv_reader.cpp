//
// Created by james on 3/28/20.
//

#include <iostream>
#include <string>
#include <sstream>
#include "csv_reader.h"


using std::getline,
      std::string,
      std::istringstream,
      std::cout,
      std::cerr,
      std::endl,
      std::ws,
      jaa::vector;

template <class T>
csv_reader<T>::csv_reader(const string &sourcefile) {
    ifstream datafile(sourcefile);
    string line;
    int next_item;

    if (datafile.is_open()) {
        while (getline(datafile, line)) {
            istringstream buffer{line};

            vector<T> row;

            while (buffer >> next_item) {
                row.push_back(next_item);

                if ((buffer >> ws).peek() == ',') {
                    buffer.ignore();
                }
            }

            data.push_back(row);
        }

        datafile.close();
    } else {
        cerr << "There was an error opening " << sourcefile << endl;
    }
}

template <class T>
vector<vector<T>>& csv_reader<T>::get_data() {
    return this->data;
}
