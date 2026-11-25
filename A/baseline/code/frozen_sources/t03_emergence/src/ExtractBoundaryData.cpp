#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "ExtractBoundaryData.h"
#include "CartesianGridAndControlPoints.h"
#include "LU.h"

void extractBoundaryDataForContinousFunction(CartesianGridAndControlPoints* G, double** u, double* bdry_u, double* bdry_un)
{
    double x0 = G->Grid_x[0];
    double y0 = G->Grid_y[0]; 

    int dir[2];
    double coord[6][2];
    double grid_u[6]; 
    double q[2]; 

    double mat[6][6]; 
    double sol[6]; 
    double c[6];
    double b[6];
    c[0] = 1.0;

    for (int flag = 0; flag < G->M; flag++) {
        double xi = G->xy[flag][0];
        double eta = G->xy[flag][1];

        int i = static_cast<int>((xi - x0) / G->dx + 0.5); 
        int j = static_cast<int>((eta - y0) / G->dy + 0.5); 

        double xi0 = x0 + i * G->dx; 
        double eta0 = y0 + j * G->dy;

        if (xi > xi0) {
            dir[0] = 1; 
        } else {
            dir[0] = - 1; 
        }
        if (eta > eta0) {
            dir[1] = 1;
        } else {
            dir[1] = - 1; 
        }

        coord[0][0] = 0.0;      coord[0][1] = 0.0; 

        coord[1][0] = 1.0;      coord[1][1] = 0.0; 
        coord[2][0] = - 1.0;    coord[2][1] = 0.0;

        coord[3][0] = 0.0;      coord[3][1] = 1.0;
        coord[4][0] = 0.0;      coord[4][1] = - 1.0;

        coord[5][0] = dir[0];   coord[5][1] = dir[1]; 


        grid_u[0] = u[i][j];
        grid_u[1] = u[i + 1][j];
        grid_u[2] = u[i - 1][j];
        grid_u[3] = u[i][j + 1];
        grid_u[4] = u[i][j - 1];
        grid_u[5] = u[i + dir[0]][j + dir[1]]; 


        q[0] = (xi - xi0) / G->dx;
        q[1] = (eta - eta0) / G->dy;

        for (i = 0; i < 6; i++) {

            b[i] = grid_u[i];

            double dx_ = coord[i][0] - q[0];
            double dy_ = coord[i][1] - q[1]; 

            c[1] = dx_; 
            c[2] = dy_; 

            c[3] = 0.5 * dx_ * dx_; 
            c[4] = 0.5 * dy_ * dy_; 

            c[5] = dx_ * dy_;

            for (j = 0; j < 6; j++) {
                mat[i][j] = c[j];
            }
        }

        bool status = LUsolve(mat, b, sol); 

        bdry_u[flag] = sol[0];
        bdry_un[flag] = sol[1] /G->dx * G->nxny[flag][0] + sol[2] /G->dy * G->nxny[flag][1];

    }
}











