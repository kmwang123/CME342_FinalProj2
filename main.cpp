#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>

#include "mesh.hpp"
#include "Array.hpp"
#include "ionTransportEqns.hpp"
#include "NSEqns.hpp"
#include "mpiWrapper.hpp"
using namespace std;

Mesh meshSetup(double Lx,double dz,int N,double beta,double dx_max,double dx_min, double Ly, double dn, int M, double dy_min);
IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, double Ey_NBC_sX, double Ey_SBC_sX, double Phi_LHS_BC_sX, double Phi_RHS_BC_sX, Mesh mesh, bool restart, bool perturb, MPI_Wrapper mpi);
NSEqns2D nsSetup(double kappa, Mesh mesh, bool restart);

int main(int argc,char** argv)  {

  /*Setup MPI*/
  const int ndim = 2;
  //const int p1 = 1;
  //const int p2 = 1;
  /*Initialize MPI*/
  //int myid, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Wrapper mpi;
  mpi.p1 = 1;//p1;
  mpi.p2 = 1;//p2;
  //cout << "mpi.myid: " << mpi.myid << endl;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  //MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  /*Boolean options*/
  bool restart = 0;  
  bool perturb = 0;

  /*Define Constants*/
  //XMesh
  const double Lx = 1;
  const double dz = 1;
  const int N = 8;//0;
  const double beta = 0.03;
  const double dx_max = 1e-2;
  const double dx_min = 1e-4;
  
  //YMesh
  const double Ly = 1;
  const double dn = 1;
  const int M = 5;//0;
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


  //Time Step
  const double dt = pow(epsilon,2);
  const double T_final = 1;
  const double Tvis = 1;
  const int numOfIterations = 3;

  /*Setup Mesh*/
  Mesh mesh = meshSetup(Lx,dz,N,beta,dx_max,dx_min,Ly,dn,M,dy_min);
  mesh.genXmesh("exponential");
  mesh.genYmesh("twoSided"); 
  mesh.printYmesh("dydn_sY"); 
  mesh.printYmesh("dydn_cY");
  mesh.printXmesh("dxdz_sX");
  mesh.printYmesh("dxdz_cX");
  /*Setup system*/
  /*Perturb C1 and C2 ions*/
  IonTransportEqns2D ionSys = ionSetup(D1,D2,epsilon,Ey_NBC_sX,Ey_SBC_sX,Phi_LHS_BC_sX,Phi_RHS_BC_sX,mesh,restart,perturb,mpi);
  NSEqns2D nsSys = nsSetup(kappa,mesh,restart);
  ionSys.printOneConcentration("C1");

  ionSys.GaussLawSolveStruct();

  MPI_Finalize();
  return 0;
}

Mesh meshSetup(double Lx,double dz,int N,double beta,double dx_max,double dx_min, double Ly, double dn, int M, double dy_min) {
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

IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, double Ey_NBC_sX, double Ey_SBC_sX, double Phi_LHS_BC_sX, double Phi_RHS_BC_sX, Mesh mesh, bool restart, bool perturb, MPI_Wrapper mpi) {
  IonTransportEqns2D ionSys(mesh,mpi);
  ionSys.D1 = D1;
  ionSys.D2 = D2;
  ionSys.epsilon = epsilon;
  ionSys.Ey_NBC_sX = Ey_NBC_sX;
  ionSys.Ey_SBC_sX = Ey_SBC_sX;
  ionSys.Phi_LHS_BC_sX = Phi_LHS_BC_sX;
  ionSys.Phi_RHS_BC_sX = Phi_RHS_BC_sX;
  ionSys.setUp(restart,perturb);
  return ionSys;
}

NSEqns2D nsSetup(double kappa, Mesh mesh, bool restart) {
  NSEqns2D nsSys(mesh);
  nsSys.kappa = kappa;
  nsSys.setUp(restart);
  return nsSys;
}
