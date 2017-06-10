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
    double C1_LHS_BC_sX;
    double C1_RHS_BC_sX;
    double C2_RHS_BC_sX;
    array2<int> C1_bc_type; 
    array2<int> C2_bc_type;
    array2<int> Ey_bc_type;
    int solver_id_gauss = 0;
 
    array2<double> C1_star;
    array2<double> C1_n; 
    array2<double> C1_nMinus1;
    array2<double> C2_star;
    array2<double> C2_n;
    array2<double> C2_nMinus1;
    array2<double> phi;
   
    array2<double> f1star_flux_sX;
    array2<double> f2star_flux_sX;
    array2<double> g1star_flux_sY;
    array2<double> g2star_flux_sY;
    array2<double> Ex_star_sX;
    array2<double> Ey_star_sY;
 
    int istartc, iendc, jstartc, jendc;
 
    array1<double>::opt halo_m;
    array1<double>::opt halo_n;

    void setUp(bool restart, bool perturb);
    void createBCarrays(int C1_bcs[2], int C2_bcs[2], int Ey_bcs[2]);
    void updateBCs(void);
    void updateFluxes(array2<double> u_star, array2<double> v_star);
    void updateInteriorFluxes(array2<double> u_star, array2<double> v_star);
    void updateBoundaryFluxes(array2<double> u_star, array2<double> v_star);
    double frand(double fMin, double fMax);
    void perturbOneConcentration(array2<double> &data);
    void printOneConcentration(string type);
    void setCstarValuesfrmCn(void);
    void updateConvergedValues(void);
    IonTransportEqns2D(Mesh& ref, MPI_Wrapper& ref2) : mesh(ref),mpi(ref2) {this-> mesh = mesh; this->mpi = mpi;} 
  private:
    Mesh& mesh;
    MPI_Wrapper& mpi;
};

#endif /*IONTRANSPORTEQNS_HPP*/
