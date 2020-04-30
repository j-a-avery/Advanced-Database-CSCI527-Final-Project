//
// Created by james on 4/25/20.
//

#include "k_means.h"
#include "matrix.h"


std::string describe(measurement m)  {
    switch(m) {
        case EUCLIDEAN_DISTANCE:
            return "Euclidean distance";
        case MANHATTAN_DISTANCE:
            return "Manhattan distance";
        case CHEBYSHEV_DISTANCE:
            return "Chebyshev distance";
        case COSINE_SIMILARITY:
            return "cosine similarity";
        case REAL_VALUED_JACCARD_SIMILARITY:
            return "real-valued Jaccard similarity";
        case CORRELATION:
            return "correlation";
        default:
            return "unknown";
    }
}

std::string describe(centroid_policy c) {
    switch(c) {
        case FULLY_RANDOM:
            return "fully randomized";
        case RANDOM_SELECTION:
            return "random selection of existing data";
        case ASSIGNED:
            return "assigned at model creation";
        default:
            return "unknown";
    }
}


// Constructor for pre-selected initial centroids
template <typename Number>
k_means<Number>::k_means(
        int k,
        matrix<Number>& data,
        matrix<Number>& initial_centroids,
        measurement measurement_metric,
        int max_epochs)
{
    // todo: data validation, esp k >= 1, centroids.rows() == k
    this->_k = k;
    this->_data = data;
    this->_initial_centroids = initial_centroids;
    this->_measurement = measurement_metric;
    this->_centroid_policy = centroid_policy::ASSIGNED;
    this->_max_epochs = max_epochs;
}


// Constructor for randomly assigned initial centroids
template <typename Number>
k_means<Number>::k_means(
        int k,
        matrix<Number>& data,
        centroid_policy assignment,
        measurement measurement_metric,
        int max_epochs)
{
    // Todo: data validation, esp k >= 1
    this->_k = k;
    this->_data = data;
    this->_measurement = measurement_metric;
    this->_centroid_policy = assignment;
    this->_max_epochs = max_epochs;
}




template <typename Number>
matrix<Number>& k_means<Number>::centroids() {
    if (!_calculated) {
        this->calculate();
    }

    return this->_centroids;
}

template <typename Number>
float k_means<Number>::SSE() {
    if (!_calculated) {
        this->calculate();
    }

    return this->_sse;
}


// Select centroids at random from the existing data
template <typename Number>
void k_means<Number>::assign_random_centroids() {
    this->_initial_centroids = new matrix<Number>(this->_k, this->_data.cols());

    std::uniform_int_distribution<int> uniform_dist(0, _data.rows());

    for (int i = 0; i < this->_k; i++) {
        this->_initial_centroids.row(i) = this->_data.row(uniform_dist(rand_dev));
    }
}


// Create centroids at random from the domain of the existing data
template <typename Number>
void k_means<Number>::create_random_centroids() {
    this->_initial_centroids = new matrix<Number>(this->_k, this->_data.cols());
    matrix<Number> colmins = this->_data.min(DIM::COL);
    matrix<Number> colmaxs = this->_data.max(DIM::COL);

    for (int row = 0; row < this->_k; row++) {
        for (int col = 0; col < this->_data.cols(); col++) {
            std::uniform_real_distribution<Number> dist((Number)colmins.at(0, col), (Number)colmaxs.at(0, col));
            this->_initial_centroids.at(row, col) = dist(rand_dev);
        }
    }
}


// Calculate Euclidean distances between two matrices
template <typename Number>
matrix<double>& euclidean_distances(matrix<Number> Data, matrix<Number> Cents) {
    // TODO: Ensure that Data and Cents have the same number of columns
    auto *distances = new matrix<double>(Data.rows(), Cents.rows());

    for (int drow = 0; drow < Data.rows(); drow++) {
        for (int crow = 0; crow < Cents.rows(); crow++) {
            distances->at(drow, crow) =
                    (Data.row(drow) - Cents.row(crow))
                    .as_double()
                    .norm()
                    .value();
        }
    }

    return *distances;
}

// Calculate Manhattan distances between two matrices
template <typename Number>
matrix<double>& manhattan_distances(matrix<Number> Data, matrix<Number> Cents) {
    // TODO: Ensure that Data and Cents hhave the same number of columns
    auto *dists = new matrix<double>(Data.rows(), Cents.rows());

    for (int drow = 0; drow < Data.rows(); drow++) {
        for (int crow = 0; crow < Cents.rows(); crow++) {
            distances->at(drow, crow) =
                    (Data.row(drow) - Cents.row(crow))
                    .as_double()
                    .abs()
                    .sum()
                    .value();
        }
    }

    return *dists;
}


template <typename Number>
void k_means<Number>::calculate() {
    return; /*
    // First, determine the centroids, if they weren't assigned at model creation
    if (_centroid_policy != centroid_policy::ASSIGNED) {
        switch (_centroid_policy) {
            case centroid_policy::FULLY_RANDOM:
                // assign random numbers within the domain of the data as centroids
                this->create_random_centroids();
                break;
            case centroid_policy::RANDOM_SELECTION:
                // randomly assign existing points as centroids
                this->assign_random_centroids();
                break;
            default:
                // This should never happen, but just in case...
                // TODO: Create an exception
                throw std::exception();
        }
    }

    this->_centroids = this->_initial_centroids.copy();
    auto previous_centroids = this->_centroids.copy();
    int epoch = 0;

    do {
        // Assign each item to its nearest centroid
        for (int i = 0; i < this->_data.rows(); i++) {
            switch (this->_measurement) {
                case measurement::EUCLIDEAN_DISTANCE:
                    this->_classes = (_centroids - _data.row(i)).square().sum(DIM::ROW).argmin(DIM::COL).value();
                    break;
                case measurement::MANHATTAN_DISTANCE:
                    this->_classes = (_centroids - _data.row(i)).abs().sum(DIM::ROW).argmin(DIM::COL).value();
                    break;
                case measurement::CHEBYSHEV_DISTANCE:
                    this->_classes = (_centroids - _data.row(i)).abs().max(DIM::ROW).argmin8(DIM::COL).value();
                    break;
                case measurement::COSINE_SIMILARITY:
                    // There's probably an even more vectorized way to do this
                    for (int j = 0; j < _k; j++) {
                        _classes.at()
                    }
                    // todo: This is wrong
                    this->_classes = _centroids.dot(_data.row(i)).value() / (_centroids.norm(DIM::ROW) * _data.norm(DIM::ROW)).value();
                    break;
                case measurement::REAL_VALUED_JACCARD_SIMILARITY:
                    break;
                case measurement::CORRELATION:
                    // Spearman correlation is the cosine similarity between A.tilde() and B.tilde()
                    break;
            }
        }

        // Recalculate centroids as the average of each cluster

    } while ((epoch < _max_epochs) && (_centroids == previous_centroids));




*/
}