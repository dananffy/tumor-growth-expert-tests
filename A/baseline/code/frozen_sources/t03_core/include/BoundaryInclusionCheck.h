#ifndef BOUNDARYINCLUSIONCHECK_H
#define BOUNDARYINCLUSIONCHECK_H

#include <algorithm>
#include <cmath>

bool isPointInPolygon(double x, double y, double** polygon, int n);
bool isBoundaryInside(double** ctrl_points_0, int n_ctrl_0,
                      double** ctrl_points_1, int n_ctrl_1);

#endif