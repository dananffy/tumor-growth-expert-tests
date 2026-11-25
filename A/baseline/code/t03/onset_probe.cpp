#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include "LinearAlgebra.h"
#include "Config.h"
#include "CartesianGridAndControlPoints0.h"
#include "MHelmholtzSolver.h"
#include "PoissonSolver.h"
int main(int argc,char** argv) {
    if(argc!=3)return 1;
    int N=std::stoi(argv[1]);double R=std::stod(argv[2]);
    lambda=1;c_B=10;n_c=1e-3;model_G0=1;c_bar=5;sqrt_lambda=1;
    omp_set_dynamic(0);omp_set_num_threads(2);
    int n=64;double** points=createMatrix(n,2);
    for(int i=0;i<n;++i){double th=2*M_PI*i/n;points[i][0]=R*cos(th);points[i][1]=R*sin(th);}
    auto* G=new CartesianGridAndControlPoints0(-5,5,-5,5,N,N,points,n,N,"points.txt","normals.txt","samples.txt",1);
    double** c=createMatrix(N+1,N+1);double** p=createMatrix(N+1,N+1);double* v=createVector(n);
    solveModifiedHelmholtz(G,c,"nutrient.txt");
    for(int i=0;i<=N;++i)for(int j=0;j<=N;++j)c[i][j]=-(c[i][j]-5);
    solvePoisson(G,c,p,v,"pressure.txt","velocity.txt");
    double minp=1e100;
    for(int i=0;i<=N;++i)for(int j=0;j<=N;++j)if(G->interior[i][j]&&!G->near_bdry[i][j]){
        if(!std::isfinite(p[i][j]))return 2;
        minp=std::min(minp,p[i][j]);
    }
    std::ofstream f("probe.json");f<<std::setprecision(17)<<"{\"grid\":"<<N<<",\"radius\":"<<R<<",\"center_pressure\":"<<p[N/2][N/2]<<",\"min_pressure_away_from_boundary\":"<<minp<<"}";
    delete G;freeMatrix(points,n);freeMatrix(c,N+1);freeMatrix(p,N+1);freeVector(v);return 0;
}
