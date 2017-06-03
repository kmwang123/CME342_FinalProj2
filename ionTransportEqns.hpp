#ifndef IONTRANSPORTEQNS_HPP
#define IONTRANSPORTEQNS_HPP

#include "mesh.hpp"
#include "Array.hpp"
using namespace Array;

class IonTransportEqns2D {

  public:

    double D1;
    double D2;
    double epsilon;
 
    array2<double> C1_star;
    array2<double> C1_n; 
    array2<double> C1_nMinus1;
    array2<double> C2_star;
    array2<double> C2_n;
    array2<double> C2_nMinus1;
    array2<double> phi;
    
    void setUp(void);
    IonTransportEqns2D(Mesh& ref) : mesh(ref){this-> mesh = mesh;} 
  private:
    Mesh& mesh;
};

#endif /*IONTRANSPORTEQNS_HPP*/
