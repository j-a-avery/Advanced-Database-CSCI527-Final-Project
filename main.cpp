// main.cpp : This file contains the 'main' function. 
// Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <limits>
#include <vector>
#include <map>
#include <chrono>

#include "matrix.h"
#include "kmeans.h"



inline std::string boolstr(bool b) {
	return b ? "true" : "false";
}


int main(int argc, char** argv) {
	matrix data = load_csv(argv[1]);

	std::map<std::string, pmeasure_func> distance_measures;
		distance_measures.insert(euclidean);
		distance_measures.insert(manhattan);
		distance_measures.insert(chebyshev);
		distance_measures.insert(cosine);
		distance_measures.insert(jaccard);
		distance_measures.insert(correlation);

	std::map<std::string, pcentroid_generator> centroid_policies;
		centroid_policies.insert(random_generator);
		centroid_policies.insert(random_selector);

	int trial = 1;
	for (int k = 3; k < 6; k++) {
		for (auto distance_measure : distance_measures) {
			for (auto centroid_policy : centroid_policies) {
				for (int run = 0; run < 5; run++) {
					std::cout << "********** TRIAL " << trial++ << " **********" << std::endl;
					std::cout << "K: " << k << std::endl
						<< "Distance Measure: " << function_name(distance_measure) << std::endl
						<< "Centroid Policy: " << function_name(centroid_policy) << std::endl;

					kmeans km(k, &data, centroid_policy, distance_measure);

					// Time execution
					auto start = std::chrono::steady_clock::now();
					km.calculate();
					auto end = std::chrono::steady_clock::now();
					std::cout << "Execution time: "
						<< std::chrono::duration<unsigned long long, std::nano>(end - start).count()
						<< "ns" << std::endl
						<< "SSE: " << km.SSE() << std::endl
						<< "Calculated in " << km.epochs_to_convergence() << " epochs."
						<< std::endl << std::endl;
				}
			}
		}
	}

	return 0;
}