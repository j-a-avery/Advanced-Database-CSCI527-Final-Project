#include <iostream>
#include <vector>

#include "matrix.h"
//#include "k_means.h"

using std::cout,
      std::setw,
      std::setprecision,
      std::showpoint,
      std::fixed,
      std::endl,
      std::flush;


matrix<double>& euclidean_distance(matrix<int>& Data, matrix<int>& Cents) {
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

matrix<double>& manhattan_distances(matrix<int>& D, matrix<int>& C) {
    auto *dists = new matrix<double>(D.rows(), C.rows());
    for (int drow = 0; drow < D.rows(); drow++) {
        for (int crow = 0; crow < C.rows(); crow++) {
            dists->at(drow, crow) =
                    (D.row(drow) - C.row(crow))
                    .as_double()
                    .abs()
                    .sum()
                    .value();
        }
    }

    return *dists;
}

matrix<double>& Chebyshev_distances(matrix<int>& D, matrix<int>& C) {
    auto *dists = new matrix<double>(D.rows(), C.rows());

    for (int drow = 0; drow < D.rows(); drow++) {
        for (int crow = 0; crow < C.rows(); crow++) {
            dists->at(drow, crow) =
                    (D.row(drow) - C.row(crow))
                    .as_double()
                    .abs()
                    .max()
                    .value();
        }
    }

    return *dists;
}


matrix<double>& cosine_similarity(matrix<int>* D, matrix<int>* C) {
    std::cout << "Allocating space on the heap for the distances. " << std::flush;
    auto *dists = new matrix<double>(D->rows(), C->rows());
    std::cout << "Done." << std::endl;
    //auto doubleD = D->as_double();
    auto doubleC = C->as_double();

    std::cout << "Calculating D ⋅ C. " << std::flush;
    matrix<double>& D_dot_C = D->as_double().dot(doubleC.T());
    std::cout << "Done." << std::endl;

    std::cout << "Calculating ||D||. " << std::flush;
    auto D_norm = D->as_double().norm(DIM::ROW);
    std::cout << "Done." << std::endl;

    std::cout << "Calculating ||C||. " << std::flush;
    auto C_norm = doubleC.norm(DIM::ROW);
    std::cout << "Done." << std::endl;

    // Surely there's a more abstract way to do this....
    for (int row = 0; row < D->rows(); row++) {
       for (int col = 0; col < C->rows(); col++) {
           dists->at(row, col) = D_dot_C.at(row, col) / (D_norm.at(row, 0) * C_norm.at(col, 0));
       }
    }

    return *dists;
}


int main() {
    matrix<int> A(2, 3);
    matrix<int> B(3, 2);
    matrix<int> Z = matrix<int>::zero(2, 3);

    A.at(0, 0) = 1; A.at(0, 1) = 2; A.at(0, 2) = 3;
    A.at(1, 0) = 4; A.at(1, 1) = 5; A.at(1, 2) = 6;

    B.at(0, 0) = 7; B.at(0, 1) = 8;
    B.at(1, 0) = 9; B.at(1, 1) = 10;
    B.at(2, 0) = 11; B.at(2, 1) = 12;

    std::cout << A << std::endl;
    std::cout << A.mean(DIM::ROW) << std::endl;
    std::cout << A.mean(DIM::COL) << std::endl;
    std::cout << A.mean(DIM::WHOLE) << std::endl;

    std::cout << A.square() << std::endl;
    std::cout << A.pow(2) << std::endl;
    std::cout << A - B.T() << std::endl;
    std::cout << (A - B.T()).abs() << std::endl;

    std::cout << A.row(0).dot(A.row(0).T()).value() << std::endl;

    std::cout << "Hadamard product between A and B.T" << std::endl;
    std::cout << A.Hadamard(B.T()) << std::endl;

    std::cout << "||A|| and Ã" << std::endl;
    std::cout << A.norm() << std::endl << A.tilde() << std::endl;

    std::cout << A.as_double().pow(-1) << std::endl;

    auto p_q = matrix<double>(2, 2);
    auto norms = matrix<double>(2,1);

    p_q.at(0, 0) = 25; p_q.at(0, 1) = 34;
    p_q.at(1, 0) = 44; p_q.at(1, 1) = 56;

    norms.at(0,0) = 0.0319;
    norms.at(1,0) = 0.01697;

    std::cout << norms.T();

    std::cout << p_q.Hadamard(norms.T());

    matrix<int> D(4, 3);
    matrix<int> C(2, 3);
    matrix<double>* E;

    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 3; col++) {
            D.at(row, col) = 1 + row + col;
        }
    }
    std::cout << D << std::endl;

    for (int row = 0; row < 2; row++) {
        for (int col = 0; col < 3; col++) {
            C.at(row, col) = 2 + row + col;
        }
    }
    std::cout << C << std::endl;

    std::cout << "====================== NORMS ======================" << std::endl;
    std::cout << A.norm(DIM::ROW) << std::endl
              << A.norm(DIM::COL) << std::endl
              << A.norm(DIM::WHOLE) << std::endl;
    std::cout << B.norm(DIM::ROW) << std::endl
              << B.norm(DIM::COL) << std::endl
              << B.norm(DIM::WHOLE) << std::endl;
    std::cout << C.norm(DIM::ROW) << std::endl
              << C.norm(DIM::COL) << std::endl
              << C.norm(DIM::WHOLE) << std::endl;
    std::cout << D.norm(DIM::ROW) << std::endl
              << D.norm(DIM::COL) << std::endl
              << D.norm(DIM::WHOLE) << std::endl;


    std::cout << "Euclidean Distance" << std::endl;
    E = &euclidean_distance(D, C);
    std::cout << *E << E->argmin(DIM::ROW) << std::endl;

    std::cout << "Manhattan Distance" << std::endl;
    E = &manhattan_distances(D, C);
    std::cout << *E << E->argmin(DIM::ROW);

    std::cout << "Chebyshev Distance" << std::endl;
    E = &Chebyshev_distances(D, C);
    std::cout << *E << E->argmin(DIM::ROW);

    std::cout << "Cosine Similarity" << std::endl;
    E = &cosine_similarity(&D, &C);
    std::cout << *E << E->argmin(DIM::ROW);


    return 0;
}
