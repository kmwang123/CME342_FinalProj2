#include <stdlib.h>
#include <string>
#include <iomanip>
#include <iostream> 
#include <math.h>
#include <time.h>

#include "Array.hpp"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "ionTransportEqns.hpp"

#define NO_FLUX_BC 0
#define DIRICHLET_BC 1
#define WALL_CHARGE_BC 2
#define LEFT_BC 0
#define RIGHT_BC 1
#define SOUTH_BC 0
#define NORTH_BC 1

using namespace std;

void IonTransportEqns2D::setUp(bool restart, bool perturb) {
  C1_star.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize); //include halo 
  C1_n.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C1_nMinus1.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_star.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_n.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_nMinus1.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  phi.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);

  istartc = mpi.hsize; 
  iendc = mesh.n+mpi.hsize-1;
  jstartc = mpi.hsize;
  jendc = mesh.m+mpi.hsize-1;

  f1star_flux_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  f2star_flux_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  g1star_flux_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));
  g2star_flux_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));

  Ex_star_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  Ey_star_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));

  if (!restart) {
    for(int i=istartc; i<=iendc; i++) {
      for (int j=jstartc; j<=jendc; j++) {
        C1_n[i][j] = 1.0;
        C2_n[i][j] = 1.0;
      }
    }
  }
  //add this later
  /*else {

  }*/

  if (perturb) {
    perturbOneConcentration(C1_n);
    perturbOneConcentration(C2_n);
  }  
  
}

void IonTransportEqns2D::createBCarrays(int C1_bcs[2], int C2_bcs[2], int Ey_bcs[2]) {
  C1_bc_type.Allocate(2,mesh.m);
  C2_bc_type.Allocate(2,mesh.m);
  Ey_bc_type.Allocate(mesh.n,2);

  // arrays for determining the bc type of the current y-location
  for(int bc_idx = 0; bc_idx < 2; bc_idx++){
    int bc_type = C1_bcs[bc_idx];
    if (bc_type == DIRICHLET_BC || bc_type == NO_FLUX_BC) {
      for(int j = 0; j < mesh.m; j++){
       C1_bc_type[bc_idx][j] = bc_type;
      }
    }
    else {
      cout << "unknown boundary condition type: " << bc_type << endl;
      cout << "defa ulting to electrode boundary condition";
      for(int j = 0; j < mesh.m; j++){
        C1_bc_type[bc_idx][j] = NO_FLUX_BC;
      }
    }
    // currently C2 doesn't have any patterning so every cell has either
    // a no flux or dirichlet b.c.
    bc_type = C2_bcs[bc_idx];
    for(int j = 0; j < mesh.m; j++){
      C2_bc_type[bc_idx][j] = bc_type;
    }
  }

  // arrays for determining the bc type of current x-location
  for (int bc_idx = 0; bc_idx < 2; bc_idx++) {
    int bc_type = Ey_bcs[bc_idx];
    if (bc_type == WALL_CHARGE_BC) {
      for (int i=0; i<mesh.n; i++) {
        Ey_bc_type[i][bc_idx] = bc_type;
      }
    }
  }

}

void IonTransportEqns2D::updateBCs(void) {
  int pi, pj;
  pi = mpi.myid % mpi.p1;
  pj = mpi.myid / mpi.p1;
  //LHS and RHS BCs
  for (int j=0; j<mesh.m; j++) {
    if (pi == mpi.p1-1) {
      // dirichlet right BC for C+ (C1)
      if (C1_bc_type[RIGHT_BC][j] == DIRICHLET_BC) {
        f1star_flux_sX[mesh.n_s-1][j] = f1star_flux_sX[mesh.n_s-1][j]
                                        -D1*mesh.dzdx_sX[mesh.N_s-1]*(2*(C1_RHS_BC_sX-C1_star[iendc][jstartc+j])/mesh.dz); //diffusion
        f1star_flux_sX[mesh.n_s-1][j] = f1star_flux_sX[mesh.n_s-1][j]
                                        -D1*C1_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(2*(Phi_RHS_BC_sX-phi[iendc][jstartc+j])/mesh.dz); //em
                                        //-D1*C1_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(-2*phi[iendc][jstartc+j])/mesh.dz; //em
      }
      // dirichlet right BC for C- (C2)
      if (C2_bc_type[RIGHT_BC][j] == DIRICHLET_BC) {
        f2star_flux_sX[mesh.n_s-1][j] = f2star_flux_sX[mesh.n_s-1][j]
                                       -D2*mesh.dzdx_sX[mesh.N_s-1]*(2*(C2_RHS_BC_sX-C2_star[iendc][jstartc+j])/mesh.dz); //diffusion
        f2star_flux_sX[mesh.n_s-1][j] = f2star_flux_sX[mesh.n_s-1][j]
                                       +D2*C2_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(2*(Phi_RHS_BC_sX-phi[iendc][jstartc+j])/mesh.dz); //em
                                       //+D2*C2_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(-2*phi[iendc][jstartc+j])/mesh.dz; //em
      }
    }
    if (pi == 0) {
      // dirichlet left bc for C+ (C1)
      if (C1_bc_type[LEFT_BC][j] == DIRICHLET_BC) {
        f1star_flux_sX[0][j] = f1star_flux_sX[0][j]
                              -D1*mesh.dzdx_sX[0]*(2*(C1_star[istartc][jstartc+j]-C1_LHS_BC_sX)/mesh.dz); //diffusion
        f1star_flux_sX[0][j] = f1star_flux_sX[0][j]
                              -D1*C1_LHS_BC_sX*mesh.dzdx_sX[0]*(2*(phi[istartc][jstartc+j]-Phi_LHS_BC_sX)/mesh.dz); //em 
      }
      // no flux left bc for C- (C2)
      if (C2_bc_type[LEFT_BC][j] == NO_FLUX_BC) {
        f2star_flux_sX[0][j] = 0.0; 
        //if (mpi.myid == 0) {
        //  cout << "fluxes: " << setprecision(15) << setw(19) << f2star_flux_sX[0][j] << endl;
        //}
     } 
    } 
  }
  //North and South BCs
  for (int i=0; i<mesh.n; i++) {
    // no flux bc for north bc C+ (C1)
    if (pj == mpi.p2-1) {
      g1star_flux_sY[i][mesh.m_s-1] = 0.0;
      g2star_flux_sY[i][mesh.m_s-1] = 0.0;  
      if (Ey_bc_type[i][NORTH_BC] == WALL_CHARGE_BC) {
        Ey_star_sY[i][mesh.m_s-1] = Ey_NBC_sX;
      }
    }
    if (pj == 0) {
      g1star_flux_sY[i][0] = 0.0;
      g2star_flux_sY[i][0] = 0.0;
      if (Ey_bc_type[i][SOUTH_BC] == WALL_CHARGE_BC) {
        Ey_star_sY[i][0] = -Ey_SBC_sX;
      }
    }
  }
}

void IonTransportEqns2D::updateInteriorFluxes(array2<double> u_star, array2<double> v_star) { 
  if (mpi.myid == 0)
  cout << "interior fluxes: " << endl;

  int pi, pj;
  pi = mpi.myid % mpi.p1;
  pj = mpi.myid / mpi.p1;

  for (int i=1; i<mesh.n_s-1; i++) {
    for (int j=0; j<mesh.m; j++) {
 
      //electric field in x-dir
      Ex_star_sX[i][j] = -(1/mesh.dz)*mesh.dzdx_sX[i]*(phi[istartc+i][jstartc+j]-phi[istartc+i-1][jstartc+j]); 

      //flux of positive ions in x-dir
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            -D1*(1/mesh.dz)*(mesh.dzdx_sX[pi*mesh.n+i])*(C1_star[istartc+i][jstartc+j]-C1_star[istartc+i-1][jstartc+j]); //diffusion
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i-1][jstartc+j])*(D1*Ex_star_sX[i][j]);//em
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i-1][jstartc+j])*(D1*u_star[i][j]);//advection
      //flux of negative ions in y-dir
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            -D1*(1/mesh.dz)*(mesh.dz)*(mesh.dzdx_sX[pi*mesh.n+i])*(C2_star[istartc+i][jstartc+j]-C2_star[istartc+i-1][jstartc+j]); //diffusion
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i-1][jstartc+j])*(D2*Ex_star_sX[i][j]);//em
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i-1][jstartc+j])*(D2*u_star[i][j]);//advection
    }
  }

  for (int i=0; i<mesh.n; i++) {
    for (int j=1; j<mesh.m_s-1; j++) {
      //electric field in y-dir
      Ey_star_sY[i][j] = -(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(phi[istartc+i][jstartc+j]-phi[istartc+i][jstartc+j-1]); 

      //flux of positive ions in y-dir
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                             -D1*(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(C1_star[istartc+i][jstartc+j]-C1_star[istartc+i][jstartc+j-1]); //diffusion
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i][jstartc+j-1])*(D1*Ey_star_sY[i][j]);//em
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i][jstartc+j-1])*(D1*v_star[i][j]);//advection 
 
      //flux of negative ions in y-dir
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           -D2*(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(C2_star[istartc+i][jstartc+j]-C2_star[istartc+i][jstartc+j-1]); //diffusion
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i][jstartc+j-1])*(D2*Ey_star_sY[i][j]);//em
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i][jstartc+j-1])*(D2*v_star[i][j]);//advection 
    }
  }
  /*for (int j=0; j<mesh.m_s; j++) {
    for (int i=0; i<mesh.n; i++) {
      if (mpi.myid == 2){
        //cout << setprecision(15) << setw(19) << Ey_star_sY[i][j] << " ";
        cout << setprecision(15) << setw(19) << g1star_flux_sY[i][j] << " ";
    }}
    cout << endl;
  }*/
  /*for (int j=0; j<mesh.m; j++) { 
    for (int i=0; i<mesh.n_s; i++) {
     if (mpi.myid == 2) {
        //cout << setprecision(15) << setw(19) << Ex_star_sX[i][j] << " ";
        cout << setprecision(15) << setw(19) << f1star_flux_sX[i][j] << " ";
      }
    }
    cout << endl;
  }*/

}

double IonTransportEqns2D::frand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f *(fMax-fMin);
}

void IonTransportEqns2D::perturbOneConcentration(array2<double> &data) {
  srand(1);
  double onePercentOfLocalValue;
  array1<double>::opt netPerturbation(mesh.n);
  array2<double> localPerturbation(mesh.n,mesh.m);
  for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      onePercentOfLocalValue = 0.00001*data[i][j];
      localPerturbation[i][j] = frand(-onePercentOfLocalValue, onePercentOfLocalValue);
      netPerturbation[i] = netPerturbation[i] + localPerturbation[i][j];
    }
  }
  for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      data[i][j] = data[i][j] + localPerturbation[i][j] - netPerturbation[i]/mesh.M; 
    }
  }

}

void IonTransportEqns2D::printOneConcentration(string type) {
  if (type == "C1") {
    cout << "C1: " << endl;
    for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
      for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
        cout << setprecision(15) << setw(19) << C1_n[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "C2") {
    cout << "C2: " << endl;
    for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
      for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
        cout << setprecision(15) << setw(19) << C2_n[i][j] << " ";
      }
      cout << endl;
    }
  }
  else { 
    cout << "Invalid type" << endl;
  }
}

void IonTransportEqns2D::setCstarValuesfrmCn(void) {
  for(int i = istartc; i <= iendc; i++){
    for(int j = jstartc; j <= jendc; j++){
      C1_star[i][j] = C1_n[i][j];
      C2_star[i][j] = C2_n[i][j];
    }
  }
}

void IonTransportEqns2D::updateConvergedValues(void) {
  for(int i = mpi.hsize; i < mesh.n+mpi.hsize; i++){
    for(int j = mpi.hsize; j < mesh.m+mpi.hsize; j++){
      C1_nMinus1[i][j] = C1_n[i][j];
      C2_nMinus1[i][j] = C2_n[i][j];
      C1_n[i][j] = C1_star[i][j];
      C2_n[i][j] = C2_star[i][j];
    }
  }
}

