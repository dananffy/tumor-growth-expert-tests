#include <iostream>
#include <cmath>
#include <stdexcept>
#include <omp.h>
#include "LinearAlgebra.h"
#include "FreeBoundaryModel.h"
#ifdef T03_EMERGENCE
#include "Config.h"
#else
#include "InterfaceBvpSolver.h"
#endif
int main(int argc,char** argv) {
    try {
        if(argc!=4) throw std::runtime_error("Usage: t03 GRID DT TEND");
        int grid=std::stoi(argv[1]); double dt=std::stod(argv[2]),end=std::stod(argv[3]);
        if(grid<64 || (grid&(grid-1)) || dt<=0 || end<=0) throw std::runtime_error("Invalid parameters");
        omp_set_dynamic(0); omp_set_num_threads(2);
        lambda=1; c_B=10; n_c=1e-3; model_G0=1;
        int n=64;
        double** outer=createMatrix(n,2);
#ifdef T03_EMERGENCE
        c_bar=5; sqrt_lambda=1;
        double radius=2.2;
#else
        setParameters(); double radius=2.5;
#endif
        for(int i=0;i<n;++i){ double a=2*M_PI*i/n; outer[i][0]=radius*cos(a);outer[i][1]=radius*sin(a); }
#ifdef T03_EMERGENCE
        testModel(-5,5,-5,5,grid,grid,outer,n,grid,dt,end,5);
#else
        double** inner=createMatrix(n,2);
        for(int i=0;i<n;++i){double a=2*M_PI*i/n;inner[i][0]=0.5751983378*cos(a);inner[i][1]=0.5751983378*sin(a);}
        int intervals=std::llround(end/dt);
        if(std::abs(intervals*dt-end)>1e-12) throw std::runtime_error("TEND must be a multiple of DT");
        testModel3(-5,5,-5,5,grid,grid,inner,n,outer,n,grid,dt,intervals+1);
        freeMatrix(inner,n);
#endif
        freeMatrix(outer,n); return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
