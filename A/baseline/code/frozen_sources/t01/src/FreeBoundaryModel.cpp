#include <fstream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>

#include "LinearAlgebra.h"
#include "PoissonSolver.h"
#include "MHelmholtzSolver.h"
#include "CartesianGridAndControlPoints.h"
#include "FreeBoundaryModel.h"
#include "CleanFiles.h"
#include "CubicSpline.h"

// Move each control point in the normal direction using the supplied boundary velocity.
void evolveCtrlPointsStrategy1(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step)
{
    for (int i = 0; i < n_ctrl; i++) {
        if (v[i] * time_step < 10) {
            ctrl_points[i][0] += v[i] * time_step * nxny[i][0];
            ctrl_points[i][1] += v[i] * time_step * nxny[i][1];
        } else {
            std::cout << "Invalid velocity!" << std::endl;
            exit(EXIT_FAILURE); 
        }
    }
}

// Advance the boundary and resample it so control points remain evenly distributed.
void evolveCtrlPointsStrategy2(double** ctrl_points, int& n_ctrl, double** nxny, double* v, double time_step)
{
    evolveCtrlPointsStrategy1(ctrl_points, n_ctrl, nxny, v, time_step);

    double* M_x = createVector(n_ctrl);
    double* alpha_x = createVector(n_ctrl);
    double* beta_x = createVector(n_ctrl);
    double* M_y = createVector(n_ctrl);
    double* alpha_y = createVector(n_ctrl);
    double* beta_y = createVector(n_ctrl);

    periodicCubicSplineCurveInterpolation(ctrl_points, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, n_ctrl);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "radius at next time step: " << computeTotalLength(n_ctrl, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y) / M_2PI << std::endl;
    SelectPointsUniformly(n_ctrl, M_x, alpha_x, beta_x, M_y, alpha_y, beta_y, ctrl_points);
    
    freeVector(M_x);
    freeVector(alpha_x);
    freeVector(beta_x);
    freeVector(M_y);
    freeVector(alpha_y);
    freeVector(beta_y);
}


// Main time-stepping driver for the one-interface tumor model.
void testModel(double x_min_, double x_max_, double y_min_, double y_max_, int I_, int J_, double** ctrl_points_, int n_ctrl_, int M_, double time_step, int steps)
{  
    double** c = createMatrix(I_ + 1, J_ + 1);
    double** p = createMatrix(I_ + 1, J_ + 1);
    double* v = createVector(n_ctrl_);

    std::ofstream metrics("radius_diagnostics.csv");
    metrics.precision(17);
    metrics << "step,time,radius_spline,radius_mean,radius_min,radius_max,velocity_mean,nutrient_min_inside,nutrient_max_inside,pressure_min_inside,pressure_max_inside\n";
    if (!metrics) throw std::runtime_error("Cannot create diagnostics");
    
    for (int i = 0; i < steps; i++) {
        
        std::cout << "Computing the cell population density and limit pressure at step " << i << ", time = " << i * time_step << std::endl;
        std::string filename_ctrlPts_ = "./result/points/ctrl_points_" + std::to_string(i) + ".txt";
        std::string filename_ctrl_nxny_ = "./result/points/ctrl_nxny_" + std::to_string(i) + ".txt";
        std::string filename_bdryPts_ = "./result/points/bdry_points_" + std::to_string(i) + ".txt";
        std::string filename_c_ = "./result/functions/c_" + std::to_string(i) + ".txt";
        std::string filename_p_ = "./result/functions/p_" + std::to_string(i) + ".txt";
        std::string filename_v_ = "./result/functions/v_" + std::to_string(i) + ".txt";
        
        CartesianGridAndControlPoints* G = new CartesianGridAndControlPoints(x_min_, x_max_, y_min_, y_max_, I_, J_, ctrl_points_, n_ctrl_, M_, filename_ctrlPts_, filename_ctrl_nxny_, filename_bdryPts_);

        solveModifiedHelmholtz(G, c, filename_c_);

        solvePoisson(G, c, p, v, filename_p_, filename_v_);

        const double radius = computeTotalLength(n_ctrl_, G->M_x, G->alpha_x, G->beta_x,
                                                 G->M_y, G->alpha_y, G->beta_y) / M_2PI;
        double rmean = 0, rmin = 1e100, rmax = 0, vmean = 0;
        for (int j = 0; j < n_ctrl_; ++j) {
            double r = std::hypot(ctrl_points_[j][0], ctrl_points_[j][1]);
            if (!std::isfinite(r) || !std::isfinite(v[j])) throw std::runtime_error("Non-finite boundary/velocity");
            rmean += r / n_ctrl_; rmin = std::min(rmin, r); rmax = std::max(rmax, r);
            vmean += v[j] / n_ctrl_;
        }
        double cmin = 1e100, cmax = -1e100, pmin = 1e100, pmax = -1e100;
        for (int j = 0; j <= I_; ++j) for (int k = 0; k <= J_; ++k) {
            if (!std::isfinite(c[j][k]) || !std::isfinite(p[j][k])) throw std::runtime_error("Non-finite field");
            if (G->interior[j][k]) {
                cmin = std::min(cmin, c[j][k]); cmax = std::max(cmax, c[j][k]);
                pmin = std::min(pmin, p[j][k]); pmax = std::max(pmax, p[j][k]);
            }
        }
        if (!std::isfinite(radius)) throw std::runtime_error("Non-finite spline radius");
        metrics << i << ',' << i * time_step << ',' << radius << ',' << rmean << ','
                << rmin << ',' << rmax << ',' << vmean << ',' << cmin << ',' << cmax << ',' << pmin << ',' << pmax << '\n';
        // State i is saved at i*dt; do not advance beyond the requested endpoint.
        if (i + 1 < steps) evolveCtrlPointsStrategy2(ctrl_points_, n_ctrl_, G->ctrlPts_nxny, v, time_step);

        std::cout << std::endl;

        delete G;
        
    }

    freeVector(v);
    freeMatrix(p, I_ + 1);
    freeMatrix(c, I_ + 1);
}
