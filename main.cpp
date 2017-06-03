#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>

#include "mesh.hpp"
#include "Array.hpp"
#include "ionTransportEqns.hpp"
#include "NSEqns.hpp"
using namespace std;

Mesh meshSetup(double Lx,double dz,int N,double beta,double dx_max,double dx_min, double Ly, double dn, int M, double dy_min);
IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, Mesh mesh);
NSEqns2D nsSetup(double kappa, Mesh mesh);

int main(int argc,char** argv)  {

  /*Setup MPI*/
  const int ndim = 2;
  const int p1 = 1;
  const int p2 = 1;
  /*Initialize MPI*/
  int myid, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  

  /*Define Constants*/
  //XMesh
  const double Lx = 1;
  const double dz = 1;
  const int N = 80;
  const double beta = 0.03;
  const double dx_max = 1e-2;
  const double dx_min = 1e-4;
  
  //YMesh
  const double Ly = 1;
  const double dn = 1;
  const int M = 50;
  const double dy_min = 1e-2;
 
  //Ion Transport Physical Constants
  const double D1 = 1;
  const double D2 = 1;
  const double epsilon = 0.01;
 
  //NS Physical Constants
  const double kappa = 0; //uncoupled for now
  //const double E = 2.0*pow(epsilon,2); //double check this
  //const double Sc = 1000;

  //Species BC
  const double C1_LHS_BC_sX = 2;
  const double C1_RHS_BC_sX = 1;
  const double C2_RHS_BC_sX = 1;
  const double Ey_wall1 = 1./epsilon; //bottom wall
  const double Ey_wall2 = 1./epsilon; //top wall
  const string PHI_BC_TYPE = "constant";
  const double Phi_LHS_BC_sX = 0;
  const double Phi_RHS_BC_sX = 60;


  //Time Step
  const double dt = pow(epsilon,2);
  const double T_final = 1;
  const double Tvis = 1;
  const int numOfIterations = 3;

  /*Setup Mesh*/
  Mesh mesh = meshSetup(Lx,dz,N,beta,dx_max,dx_min,Ly,dn,M,dy_min);
  mesh.genXmesh("exponential");
  mesh.genYmesh("uniform"); 
  //mesh.printXmesh("x_vect_sX");
  //mesh.printYmesh("y_vect_sY"); 
  /*Setup system*/
  IonTransportEqns2D ionSys = ionSetup(D1,D2,epsilon,mesh);
  NSEqns2D nsSys = nsSetup(kappa,mesh);
  /*int status = sys.readInputFile(inputfile);
  if (status) {
    cerr << "ERROR: System setup was unsuccessful!" << endl;
    return 1;
  }*/

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

IonTransportEqns2D ionSetup(double D1, double D2, double epsilon, Mesh mesh) {
  IonTransportEqns2D ionSys(mesh);
  ionSys.D1 = D1;
  ionSys.D2 = D2;
  ionSys.epsilon = epsilon;
  ionSys.setUp();
  return ionSys;
}

NSEqns2D nsSetup(double kappa, Mesh mesh) {
  NSEqns2D nsSys(mesh);
  nsSys.kappa = kappa;
  nsSys.setUp();
  return nsSys;
}
