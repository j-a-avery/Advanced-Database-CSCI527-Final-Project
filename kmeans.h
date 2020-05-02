// Functions for calculating k-means
#pragma once

#include <random>
#include <cmath>
#include <map>
#include <string>

#include "matrix.h"


const double PI = std::acos(-1);	// per cppreference.com


// Typedefs to make functional-style programming less verbose
typedef matrix& measure_function(matrix&, matrix&);
typedef measure_function* pmeasure_func;
typedef std::pair<std::string, pmeasure_func> named_measure_function;
typedef matrix& centroid_generator(int, matrix&);
typedef centroid_generator* pcentroid_generator;
typedef std::pair<std::string, pcentroid_generator> named_centroid_generator;



// Generate random centroids that fall within the range of the data,
// but are not necessarily points found within the data
matrix& generate_random_centroids(int k, matrix& data) {
	// Set up random number generator
	static std::random_device r;
	static std::default_random_engine e(r());
	
	// Find the min and max of each column
	matrix colmins(1, data.cols());
	matrix colmaxs(1, data.cols());

	for (int col = 0; col < data.cols(); col++) {
		colmins.at(0, col) = std::numeric_limits<double>::max();
		colmaxs.at(0, col) = std::numeric_limits<double>::min();

		for (int row = 0; row < data.rows(); row++) {
			if (data.at(row, col) > colmaxs.at(0, col)) {
				colmaxs.at(0, col) = data.at(row, col);
			}

			if (data.at(row, col) < colmins.at(0, col)) {
				colmins.at(0, col) = data.at(row, col);
			}
		}
	}

	// Generate centroids from within the ranges
	matrix* centroids = new matrix(k, data.cols());

	for (int row = 0; row < k; row++) {
		for (int col = 0; col < data.cols(); col++) {
			std::uniform_real_distribution<double> 
				rand(colmins.at(0, col), colmaxs.at(0, col));
			centroids->at(row, col) = rand(e);
		}
	}

	centroids->setprecision(data.precision());
	return *centroids;
}


// Select k random data points to serve as centroids
matrix& select_random_centroids(int k, matrix& data) {
	// Set up random number generator
	static std::random_device r;
	static std::default_random_engine e(r());
	std::uniform_int_distribution<int> randint(0, data.rows()-1);

	matrix* centroids = new matrix(k, data.cols());

	for (int i = 0; i < k; i++) {
		int row = randint(e);

		// Copy the data row into the centroid
		for (int col = 0; col < data.cols(); col++) {
			centroids->at(i, col) = data.at(row, col);
		}
	}

	centroids->setprecision(data.precision());
	return *centroids;
}


// Calculate the euclidean distances between two matrices
// Given matrix(m, n) A and matrix(p, q) B,
// the result will be a matrix(m, p) whose contents
// for row r, column c, are Euclidean_distance(A[row r], B[row c])
matrix& Euclidean_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* dists = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			double dist = 0;
			for (int col = 0; col < a.cols(); col++) {
				dist += std::pow(a.at(arow, col) - b.at(brow, col), 2.0);
			}
			dists->at(arow, brow) = std::sqrt(dist);
		}
	}

	dists->setprecision(std::max({ a.precision(), b.precision(), 2 }));
	return *dists;
}


// Calculate Manhattan distances between two matrices, as above
// Manattan Distance = sum(abs(a-b)) for a, b in A, B
matrix& Manhattan_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* dists = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			dists->at(arow, brow) = 0;
			for (int col = 0; col < a.cols(); col++) {
				dists->at(arow, brow) += std::abs(a.at(arow, col) - b.at(brow, col));
			}
		}
	}

	dists->setprecision(0);
	return *dists;
}


// Calculate Chebyshev distances between two matrices, as above
// Chebyshev distance = max(abs(a-b)) for a, b in A, B
matrix& Chebyshev_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* dists = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			dists->at(arow, brow) = std::abs(a.at(arow, 0) - b.at(brow, 0));

			for (int col = 1; col < a.cols(); col++) {
				double newdist = std::abs(a.at(arow, col) - b.at(brow, col));
				if (newdist > dists->at(arow, brow)) {
					dists->at(arow, brow) = newdist;
				}
			}
		}
	}

	dists->setprecision(0);
	return *dists;
}


// Calculate angular cosine distance between the two matrices, as above.
// Angular-cosine-distance = (1 / cosine-similarity) / pi,
// where cosine-similarity = dot(A-row, B-row) / (norm(A-row) * norm(B-row))
// where dot(V, W) is the dot product between V and W,
//		sum(vw) for v, w in V, W
// and norm(V) is the L2 norm, sqrt(sum(vv)) for v in V
matrix& angular_cosine_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* sims = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			double dotprod = 0;
			double arownorm = 0;
			double brownorm = 0;
			for (int col = 0; col < a.cols(); col++) {
				dotprod += a.at(arow, col) * b.at(brow, col);
				arownorm += std::pow(a.at(arow, col), 2.0);
				brownorm += std::pow(b.at(brow, col), 2.0);
			}

			arownorm = std::sqrt(arownorm);
			brownorm = std::sqrt(brownorm);

			double sim = dotprod / (arownorm * brownorm);

			sims->at(arow, brow) = (1 / sim) / PI;
		}
	}

	sims->setprecision(std::max({ a.precision(), b.precision(), 2 }));
	return *sims;
}


// Calculate real-valued Jaccard distances (Aggarwal, 2018), as above
// Jaccard similarity is ordinarily a set-based similarity metric
// defined for sets S and T as 
//		|S intersect T| / |S union T|
// but this can be abstracted for real-valued vectors V and W as
//		dot(V, W) / (dot(V, V) + dot(W, W) - dot(V, W))
// Jaccard distance is (1 - Jaccard-similarity)
matrix& Jaccard_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* sims = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			double adotb = 0;
			double adota = 0;
			double bdotb = 0;

			for (int col = 0; col < a.cols(); col++) {
				adota += std::pow(a.at(arow, col), 2.0);
				bdotb += std::pow(b.at(brow, col), 2.0);
				adotb += a.at(arow, col) * b.at(brow, col);
			}

			sims->at(arow, brow) = adotb / (adota + bdotb - adotb);
			sims->at(arow, brow) = 1 - sims->at(arow, brow);
		}
	}

	sims->setprecision(std::max({ a.precision(), b.precision(), 2 }));
	return *sims;
}


// Calculate correlation distances
// Per Boyd & Vandenberghe (2018), correlation coefficient can also be used
// as an effective similarity metric between vectors.
// They introduce a non-standard notation, A-tilde, 
// which is simply A - mean(A) (i.e. for a in A, A-tilde = a-mean(A)).
// Under this notation, correlation coefficient rho
// between two vectors A and B is simply
// the cosine similarity between A-tilde and B-tilde
// Because correlation coefficients sit in -1 <= rho <= 1,
// Converting correlation coefficient to distance requires
// d_rho = 1 - ( (rho / 2) + (1/2) )
matrix& correlation_coefficient_distances(matrix& a, matrix& b) {
	if (a.cols() != b.cols()) {
		// TODO: Make this throw an exception
		return *(new matrix());
	}

	matrix* corrs = new matrix(a.rows(), b.rows());

	for (int arow = 0; arow < a.rows(); arow++) {
		for (int brow = 0; brow < b.rows(); brow++) {
			double abar = 0;
			double bbar = 0;

			for (int col = 0; col < a.cols(); col++) {
				abar += a.at(arow, col);
				bbar += b.at(brow, col);
			}

			abar /= a.cols();
			bbar /= b.cols();

			double dotprod = 0;
			double arownorm = 0;
			double brownorm = 0;

			for (int col = 0; col < a.cols(); col++) {
				double atilde = a.at(arow, col) - abar;
				double btilde = b.at(brow, col) - bbar;

				dotprod += atilde * btilde;
				arownorm += atilde * atilde;
				brownorm += btilde * btilde;
			}

			arownorm = std::sqrt(arownorm);
			brownorm = std::sqrt(brownorm);

			double sim = dotprod / (arownorm * brownorm);
			double dist = 1 - (sim / 2 + 0.5);

			corrs->at(arow, brow) = dist;
		}
	}

	corrs->setprecision(std::max({a.precision(), b.precision(), 2 }));
	return *corrs;
}

named_centroid_generator random_generator("generate random centroids", generate_random_centroids);
named_centroid_generator random_selector("select random centroids", select_random_centroids);

named_measure_function euclidean("Euclidean distance", Euclidean_distances);
named_measure_function manhattan("Manhattan distance", Manhattan_distances);
named_measure_function chebyshev("Chebyshev distance", Chebyshev_distances);
named_measure_function cosine("angular cosine distance", angular_cosine_distances);
named_measure_function jaccard("Jaccard distance", Jaccard_distances);
named_measure_function correlation("correlation coefficient distance", correlation_coefficient_distances);

std::string function_name(named_measure_function f) {
	return f.first;
}

std::string function_name(named_centroid_generator f) {
	return f.first;
}

pmeasure_func function_call(named_measure_function f) {
	return f.second;
}

pcentroid_generator function_call(named_centroid_generator f) {
	return f.second;
}


// Kmeans class
// Calculates k-means using provided data and selected distance metric.
// Stores SSE and clusters for later access.
class kmeans {
private:
	int k;
	matrix* data;
	named_centroid_generator centroid_policy;
	named_measure_function measurement_function;
	int max_epochs;
	int epochs;
	double sse = 0;
	matrix* centroids;
	matrix** clusters;

public:
	kmeans(
		int k,
		matrix* data,
		named_centroid_generator centroid_policy,
		named_measure_function measurement_function,
		int max_epochs = 100);
	void calculate();
	inline int epochs_to_convergence() const { return epochs; }
	inline double SSE() const { return sse; }

};

kmeans::kmeans(
	int k,
	matrix* data,
	named_centroid_generator centroid_policy,
	named_measure_function measurement_function,
	int max_epochs)
{
	this->k = k;
	this->data = data;
	this->centroid_policy = centroid_policy;
	this->measurement_function = measurement_function;
	this->max_epochs = max_epochs;
	this->centroids = &function_call(centroid_policy)(k, *data);
}


void kmeans::calculate() {
	matrix previous_centroids(k, centroids->cols());	// Hold previous centroids for comparison
	epochs = 0;											// Count the number of epochs to convergence
	matrix cluster_sizes(1, k);							// The size of each cluster, for averaging centroids
	bool centroids_changed;								// since this is used twice (once for the algorithm and once for output), and matrix==matrix is O(mn)
	matrix cluster_membership(data->rows(), 1);			// Hold which cluster each data point belongs to

	do {
		std::cout << "\rCalculating clusters: Epoch " << std::right << std::setw(4) << ++epochs << std::flush;
		// Back up and zero out the centroids matrix
		for (int row = 0; row < centroids->rows(); row++) {
			for (int col = 0; col < centroids->cols(); col++) {
				previous_centroids.at(row, col) = centroids->at(row, col);
				centroids->at(row, col) = 0;
			}
		}

		// Zero out the cluster_sizes matrix
		for (int col = 0; col < k; col++) {
			cluster_sizes.at(0, col) = 0;
		}


		// Calculate the distances between each data point and each centroid
		matrix distances = function_call(measurement_function)(*data, *centroids);
		// This could also have been written as
		// matrix distances = measurement_function.second(*data, *centroids),
		// but the wrapper function makes the intent clearer, I think.
		
		// Based on these distances, determine cluster membership
		for (int row = 0; row < distances.rows(); row++) {
			// Find the index of the smallest distance for this row
			double rowmin = std::numeric_limits<double>::max();
			for (int col = 0; col < distances.cols(); col++) {
				if (distances.at(row, col) < rowmin) {
					rowmin = col;
				}
			}

			cluster_membership.at(row, 0) = rowmin;

			// Update the centroid with the smallest distance
			for (int col = 0; col < centroids->cols(); col++) {
				centroids->at(rowmin, col) += data->at(row, col);
			}

			// Increment the cluster size
			cluster_sizes.at(0, rowmin)++;
		}

		// Average the centroids by dividing each centroid data point by the size of that centroids
		for (int centroid = 0; centroid < k; centroid++) {
			for (int col = 0; col < centroids->cols(); col++) {
				centroids->at(centroid, col) /= cluster_sizes.at(0, centroid);
			}
		}

		// Are we finished yet?
		centroids_changed = *centroids != previous_centroids;
	} while (centroids_changed && epochs < max_epochs);

	// Report the results
	if (centroids_changed) {
		std::cout << std::endl << "\rMaximum number of epochs reached without convergence." << std::endl;
	}
	else {
		std::cout << std::endl << "\rConvergence reached after " << epochs << " epochs." << std::endl;
	}

	// Create the clusters matrix
	std::cout << "Building clusters. " << std::flush;
	clusters = new matrix*[k];
	for (int i = 0; i < k; i++) {
		clusters[i] = new matrix(cluster_sizes.at(0, i), data->cols());
	}

	// A matrix to track indices while building each cluster
	matrix cluster_row(1, k);
	for (int i = 0; i < k; i++) {
		cluster_row.at(0, i) = 0;
	}

	// Copy data to the clusters
	for (int data_row = 0; data_row < data->rows(); data_row++) {
		int cluster = cluster_membership.at(data_row, 0);
		int this_cluster_row = cluster_row.at(0, cluster);
		for (int col = 0; col < data->cols(); col++) {
			clusters[cluster]->at(this_cluster_row, col) = data->at(data_row, col);
		}
		cluster_row.at(0, cluster) += 1;
	}

	std::cout << "Done." << std::endl << "Calculating Sum of Squared Error. " << std::flush;

	// Calculate each cluster's mean
	matrix cluster_means(k, data->cols());

	for (int cluster = 0; cluster < k; cluster++) {
		// Calculate the column-wise sum of each cluster
		for (int row = 0; row < clusters[cluster]->rows(); row++) {
			for (int col = 0; col < data->cols(); col++) {
				cluster_means.at(cluster, col) += clusters[cluster]->at(row, col);
			}
		}

		// Divide the column-wise sums by the number of items in the cluster
		for (int col = 0; col < data->cols(); col++) {
			cluster_means.at(cluster, col) /= clusters[cluster]->rows();
		}
	}

	// Calculate Sum of Squared Error
	sse = 0;
	for (int row = 0; row < data->rows(); row++) {
		// I did not consider SSE when I created the distance functions,
		// so the distance functions as written are not suited for one-row calculations.
		// The easiest way to take care of this will be to pass 1-row matrices
		// to the distance functions, and use the value of the resulting 1x1 matrix
		// That will be slow, with all the extra copying and destroying,
		// but it will get the job done
		matrix datarow(1, data->cols());
		matrix meanrow(1, data->cols());

		for (int col = 0; col < data->cols(); col++) {
			datarow.at(0, col) = data->at(row, col);
			meanrow.at(0, col) = cluster_means.at(cluster_membership.at(row, 0), col);
		}

		double dist = function_call(measurement_function)(datarow, meanrow).at(0, 0);

		sse += dist * dist;
	}

	std::cout << "Done." << std::endl;
}