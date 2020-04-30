//
// Created by james on 4/25/20.
//

#ifndef PROJECT_K_MEANS_H
#define PROJECT_K_MEANS_H

#include "matrix.h"
#include <string>
#include <random>

enum measurement {
    EUCLIDEAN_DISTANCE,
    MANHATTAN_DISTANCE,
    CHEBYSHEV_DISTANCE,
    COSINE_SIMILARITY,
    REAL_VALUED_JACCARD_SIMILARITY,
    CORRELATION
};

enum centroid_policy {
    FULLY_RANDOM,
    RANDOM_SELECTION,
    ASSIGNED
};

std::string describe(measurement m);


template <typename Number>
class k_means {
private:
    bool _calculated = false;
    int _k;
    matrix<Number> _data;
    matrix<Number>* _initial_centroids;
    matrix<Number> _centroids;
    measurement _measurement;
    centroid_policy _centroid_policy;
    matrix<Number> _distances;
    matrix<Number> _classes;
    double _sse;
    int _max_epochs;
    static std::random_device rand_dev;
    static std::default_random_engine engine;

    void create_random_centroids();
    void assign_random_centroids();

    void euclidean_distances();
    void manhattan_distances();
    void chebyshev_distances();
    void cosine_similarities();
    void correlations();
    void jaccard_similarities();

public:
    // Constructor for pre-selected initial centroids
    k_means<Number>(
            int k,
            matrix<Number>& data,
            matrix<Number>& initial_centroids,
            measurement measurement_metric,
            int max_epochs = 20);

    // Constructor for randomly assigned initial centroids
    k_means<Number>(
            int k,
            matrix<Number>& data,
            centroid_policy assignment,
            measurement measurement_metric,
            int max_epochs = 20);

    matrix<Number>& centroids();
    float SSE();
    void calculate();
};

template class k_means<int>;
template class k_means<double>;


#endif //PROJECT_K_MEANS_H
