#ifndef MESH_HPP
#define MESH_HPP

#include <math.h>
#include <string>

#include "Array.hpp"
using namespace Array;
//using namespace std;

class Mesh{
  public:
    int N;
    int N_s;
    int M;
    int M_s;
    double Ly;
    double Lx;
    double beta;
    double dx_min;
    double dx_max;
    double dy_min;
    double dz;
    double dn;

    array1<double>::opt x_vect_sX;
    array1<double>::opt x_vect_cX;
    array1<double>::opt y_vect_sY;
    array1<double>::opt y_vect_cY;
    
    array1<double>::opt dxdz_sX;
    array1<double>::opt dxdz_cX;
    array1<double>::opt dzdx_sX;
    array1<double>::opt dzdx_cX;
    array1<double>::opt dydn_sY;
    array1<double>::opt dydn_cY;
    array1<double>::opt dndy_sY;
    array1<double>::opt dndy_cY;
    
    void genXmesh(std::string type);
    void genYmesh(std::string type);
    void printXmesh(std::string type);
    void printYmesh(std::string type);
    private:
      class GammaFunctorTwoSided {
        public:
          double dF_min_des;
          double n;
          GammaFunctorTwoSided(double init_val,int nx) : dF_min_des(init_val), n(nx) {} 
          double operator()(double gam) {
            return gam/(n*(tanh(gam/2)-tanh(-gam/2)))*pow(1/cosh(-gam/2),2)-dF_min_des;
          }
      };
      class GammaFunctorOneSided {
        public:
          double dF_min_des;
          double n;
          GammaFunctorOneSided(double init_val,int nx) : dF_min_des(init_val), n(nx) {}
          double operator()(double gam) {
            return tanh(gam/2.*(1./n - 1.))/tanh(gam/2.) + 1. - dF_min_des;
          }
      }; 

};

#endif /*MESH_HPP*/
