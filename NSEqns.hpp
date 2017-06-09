#ifndef NSEQNS_HPP
#define NSEQNS_HPP

#include "mesh.hpp"
#include "Array.hpp"
#include "mpiWrapper.hpp"
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

    void setUp(bool restart);
    void updateConvergedValues(void);
    void updateBCs(void);
    NSEqns2D(Mesh& ref, MPI_Wrapper& ref2) : mesh(ref),mpi(ref2) {this->mesh = mesh; this->mpi = mpi;}
  private:
    Mesh& mesh;
    MPI_Wrapper& mpi;
};
#endif /*NSEQNS_HPP*/
