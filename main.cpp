#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>

#include "mesh.hpp"
#include "Array.hpp"
#include "ionTransportEqns.hpp"
#include "NSEqns.hpp"
#include "mpiWrapper.hpp"
#include "hypreSolver.hpp"

using namespace std;

Mesh meshSetup(double Lx,double dz,int N,double beta,double dx_max,double dx_min, double Ly, double dn, int M, double dy_min, int p1, int p2);
IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, double Ey_NBC_sX, double Ey_SBC_sX, double Phi_LHS_BC_sX, double Phi_RHS_BC_sX, double C1_LHS_BC_sX, double C1_RHS_BC_sX, double C2_RHS_BC_sX, Mesh mesh, bool restart, bool perturb, MPI_Wrapper mpi);
NSEqns2D nsSetup(double kappa, Mesh mesh, bool restart, MPI_Wrapper mpi);

int main(int argc,char** argv)  {

  /*Setup MPI*/
  const int ndim = 2;
  //const int p1 = 1;
  //const int p2 = 1;
  /*Initialize MPI*/
  //int myid, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Wrapper mpi;
  mpi.p1 = 2;//p1;
  mpi.p2 = 3;//p2;
  mpi.pi = mpi.myid % mpi.p1;
  mpi.pj = mpi.myid / mpi.p1;
  mpi.getNeighbors();
  //MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  //MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  /*Boolean options*/
  bool restart = 0;  
  bool perturb = 0;

  /*Define Constants*/
  //XMesh
  const double Lx = 1;
  const double dz = 1;
  const int N = 16;//0;
  const double beta = 0.03;
  const double dx_max = 1e-2;
  const double dx_min = 1e-4;
  
  //YMesh
  const double Ly = 1;
  const double dn = 1;
  const int M = 18;//0;
  const double dy_min = 1e-2;
 
  //Ion Transport Physical Constants
  const double D1 = 1;
  const double D2 = 1;
  const double epsilon = 0.1;
 
  //NS Physical Constants
  const double kappa = 0; //uncoupled for now
  //const double E = 2.0*pow(epsilon,2); //double check this
  //const double Sc = 1000;

  //Species BC
  const double C1_LHS_BC_sX = 2;
  const double C1_RHS_BC_sX = 1;
  const double C2_RHS_BC_sX = 1;
  const double Ey_NBC_sX = 1./epsilon; //top wall
  const double Ey_SBC_sX = 1./epsilon; //bottom wall
  const string PHI_BC_TYPE = "constant";
  const double Phi_LHS_BC_sX = -60;
  const double Phi_RHS_BC_sX = 0;
  #define NO_FLUX_BC 0
  #define DIRICHLET_BC 1
  #define WALL_CHARGE_BC 2
  int C1_bcs[2] = {DIRICHLET_BC, DIRICHLET_BC};
  //////// C2 BCS //////////////////
  int C2_bcs[2] = {NO_FLUX_BC, DIRICHLET_BC}; 
  /*NOTE: Side wall bc's automatically enforced (flux = 0 through walls)*/
  /*NOTE: No slip and penetration enforced for velocity bcs*/   
  int Ey_bcs[2] = {WALL_CHARGE_BC,WALL_CHARGE_BC}; //{north,south} 

  //Time Step
  const double dt = pow(epsilon,2);
  const double T_final = pow(epsilon,2);
  const double Tvis = 1;
  const int numOfIterations = 1;//3; NEED TO CHANGE BACK!!!
  const int numOfTimeSteps = int(T_final/dt);

  /*Setup Mesh*/
  Mesh mesh = meshSetup(Lx,dz,N,beta,dx_max,dx_min,Ly,dn,M,dy_min,mpi.p1,mpi.p2);
  mesh.genXmesh("exponential");
  mesh.genYmesh("twoSided"); 

  //mesh.printYmesh("y_vect_sY");
  //mesh.printYmesh("dndy_sY");
  
 
  /*Setup system*/
  /*Perturb C1 and C2 ions*/
  IonTransportEqns2D ionSys = ionSetup(D1,D2,epsilon,Ey_NBC_sX,Ey_SBC_sX,Phi_LHS_BC_sX,Phi_RHS_BC_sX,C1_LHS_BC_sX,C1_RHS_BC_sX,C2_RHS_BC_sX,mesh,restart,perturb,mpi);
  NSEqns2D nsSys = nsSetup(kappa,mesh,restart,mpi);
  //ionSys.printOneConcentration("C1");

  //setup bc arrays
  ionSys.createBCarrays(C1_bcs,C2_bcs,Ey_bcs);

  GaussLawSolveStruct(mesh,mpi,ndim,ionSys.C1_n, ionSys.C2_n,ionSys.phi,epsilon, Ey_SBC_sX, Ey_NBC_sX, Phi_LHS_BC_sX,Phi_RHS_BC_sX,1);

  // set first star (guess) value with the initial concentration
  ionSys.setCstarValuesfrmCn();

  /////////////////////// Begin Time Stepping ///////////////////
  for (int time_i=1; time_i<numOfTimeSteps+1; time_i++) {
    ////////////// Iterate within a timestep ////////////////
    for (int k=0; k<numOfIterations; k++) {
      ionSys.updateBCs();
      //nsSys.updateBCs();
      
      //send halo cell-centered values, update interior flux, update boundary flux
      ionSys.sendCenters_updateFluxes(nsSys.u_star,nsSys.v_star);
      //send halo fluxes, update interior rhs, update boundary rhs
      //ionSys.sendFluxes_updateRHS(time_i);
      //ionSys.updateBoundaryFluxes(nsSys.u_star,nsSys.v_star);
    }
    //////////////////// End Iteration ///////////////////////
    ionSys.updateConvergedValues();
    //nsSys.updateConvergedValues();
  }


  MPI_Finalize();
  return 0;
}

Mesh meshSetup(double Lx,double dz,int N,double beta,double dx_max,double dx_min, double Ly, double dn, int M, double dy_min, int p1, int p2) {
  Mesh mesh;
  mesh.Lx = Lx;
  mesh.dz = dz;
  mesh.N = N;
  mesh.N_s = N+1;
  mesh.beta = beta;
  mesh.dx_max = dx_max;
  mesh.dx_min = dx_min;
  mesh.Ly = Ly;
  mesh.dn = dn;
  mesh.M = M;
  mesh.M_s = M+1;
  mesh.dy_min = dy_min;
  mesh.n = N/p1;
  mesh.n_s = mesh.n+1;
  mesh.m = M/p2;
  mesh.m_s = mesh.m+1;
  Allocate(mesh.x_vect_cX,mesh.N);
  Allocate(mesh.x_vect_sX,mesh.N_s);
  Allocate(mesh.y_vect_cY,mesh.M);
  Allocate(mesh.y_vect_sY,mesh.M_s);

  Allocate(mesh.dxdz_sX,mesh.N_s);
  Allocate(mesh.dxdz_cX,mesh.N);
  Allocate(mesh.dzdx_sX,mesh.N_s);
  Allocate(mesh.dzdx_cX,mesh.N);

  Allocate(mesh.dydn_sY,mesh.M_s);
  Allocate(mesh.dydn_cY,mesh.M);
  Allocate(mesh.dndy_sY,mesh.M_s);
  Allocate(mesh.dndy_cY,mesh.M);
  return mesh;
}

IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, double Ey_NBC_sX, double Ey_SBC_sX, double Phi_LHS_BC_sX, double Phi_RHS_BC_sX, double C1_LHS_BC_sX, double C1_RHS_BC_sX, double C2_RHS_BC_sX, Mesh mesh, bool restart, bool perturb, MPI_Wrapper mpi) {
  IonTransportEqns2D ionSys(mesh,mpi);
  ionSys.D1 = D1;
  ionSys.D2 = D2;
  ionSys.epsilon = epsilon;
  ionSys.Ey_NBC_sX = Ey_NBC_sX;
  ionSys.Ey_SBC_sX = Ey_SBC_sX;
  ionSys.Phi_LHS_BC_sX = Phi_LHS_BC_sX;
  ionSys.Phi_RHS_BC_sX = Phi_RHS_BC_sX;
  ionSys.C1_LHS_BC_sX = C1_LHS_BC_sX;
  ionSys.C1_RHS_BC_sX = C1_RHS_BC_sX;
  ionSys.C2_RHS_BC_sX = C2_RHS_BC_sX;
  ionSys.setUp(restart,perturb);
  return ionSys;
}

NSEqns2D nsSetup(double kappa, Mesh mesh, bool restart, MPI_Wrapper mpi) {
  NSEqns2D nsSys(mesh,mpi);
  nsSys.kappa = kappa;
  nsSys.setUp(restart);
  return nsSys;
}
