#ifndef NSEQNS_HPP
#define NSEQNS_HPP

#include "mesh.hpp"
#include "Array.hpp"
using namespace Array;

class NSEqns2D {

  public:
    //double Sc;
    double kappa;
    
    array2<double> u_star;
    array2<double> u_n;
    array2<double> u_nMinus1;
    array2<double> v_star;
    array2<double> v_n;
    array2<double> v_nMinus1;

    void setUp(void);
    NSEqns2D(Mesh& ref) : mesh(ref){this->mesh = mesh;}
  private:
    Mesh& mesh;
};
#endif /*NSEQNS_HPP*/
