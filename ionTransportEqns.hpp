#ifndef IONTRANSPORTEQNS_HPP
#define IONTRANSPORTEQNS_HPP

#include <string>

#include "mesh.hpp"
#include "mpiWrapper.hpp"
#include "Array.hpp"
using namespace Array;
using namespace std;
class IonTransportEqns2D {

  public:
    int ndim = 2;
    double D1;
    double D2;
    double epsilon;
    double Ey_NBC_sX;
    double Ey_SBC_sX;
    double Phi_LHS_BC_sX;
    double Phi_RHS_BC_sX;
    int solver_id_gauss = 0;
 
    array2<double> C1_star;
    array2<double> C1_n; 
    array2<double> C1_nMinus1;
    array2<double> C2_star;
    array2<double> C2_n;
    array2<double> C2_nMinus1;
    array2<double> phi;
    
    void setUp(bool restart, bool perturb);
    double frand(double fMin, double fMax);
    void perturbOneConcentration(array2<double> &data);
    void printOneConcentration(string type);
    void GaussLawSolveStruct(void);
    void GaussLawSolveSStruct(void);
    IonTransportEqns2D(Mesh& ref, MPI_Wrapper ref2) : mesh(ref),mpi(ref2) {this-> mesh = mesh; this->mpi = mpi;} 
  private:
    Mesh& mesh;
    MPI_Wrapper& mpi;
};

#endif /*IONTRANSPORTEQNS_HPP*/
