#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <omp.h>
#include "LinearAlgebra.h"
#include "MHelmholtzSolver.h"
#include "PoissonSolver.h"
#include "FreeBoundaryModel.h"

// CMAME Example 1. The solver is copied from the pinned upstream revision;
// only this parameter entry point and output/endpoint handling are adapted.
int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("Usage: t01 GRID DT FINAL_TIME THREADS");
        const int grid = std::stoi(argv[1]);
        const double dt = std::stod(argv[2]);
        const double end = std::stod(argv[3]);
        const int threads = std::stoi(argv[4]);
        const int intervals = static_cast<int>(std::llround(end / dt));
        if (grid < 32 || (grid & (grid - 1)) || dt <= 0 || end <= 0 ||
            threads < 1 || std::abs(intervals * dt - end) > 1e-12)
            throw std::runtime_error("Invalid grid/time/thread configuration");
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
        lambda = 0.5;
        update_sqrt_lambda();
        c_B = 10.0;
        G0 = 0.2;
        const int n = 32;
        double** points = createMatrix(n, 2);
        for (int i = 0; i < n; ++i) {
            const double theta = i * 2.0 * M_PI / n;
            points[i][0] = 0.5 * std::cos(theta);
            points[i][1] = 0.5 * std::sin(theta);
        }
        std::cout << "T01 Example 1: grid=" << grid << ", dt=" << dt
                  << ", final_time=" << end << ", threads=" << threads << std::endl;
        testModel(-1.5, 1.5, -1.5, 1.5, grid, grid, points, n,
                  4 * grid, dt, intervals + 1);
        freeMatrix(points, n);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
