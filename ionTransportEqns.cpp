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
  C1_star.Allocate(mesh.N,mesh.M);
  C1_n.Allocate(mesh.N,mesh.M);
  C1_nMinus1.Allocate(mesh.N,mesh.M);
  C2_star.Allocate(mesh.N,mesh.M);
  C2_n.Allocate(mesh.N,mesh.M);
  C2_nMinus1.Allocate(mesh.N,mesh.M);
  phi.Allocate(mesh.N,mesh.M);

  f1star_flux_sX.Allocate(mesh.N_s,mesh.M);
  f2star_flux_sX.Allocate(mesh.N_s,mesh.M);
  g1star_flux_sY.Allocate(mesh.N,mesh.M_s);
  g2star_flux_sY.Allocate(mesh.N,mesh.M_s);

  if (!restart) {
    for(int i=0; i<mesh.N; i++) {
      for (int j=0; j<mesh.M; j++) {
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
  C1_bc_type.Allocate(2,mesh.M);
  C2_bc_type.Allocate(2,mesh.M);
  Ey_bc_type.Allocate(mesh.N,2);

  // arrays for determining the bc type of the current y-location
  for(int bc_idx = 0; bc_idx < 2; bc_idx++){
    int bc_type = C1_bcs[bc_idx];
    if (bc_type == DIRICHLET_BC || bc_type == NO_FLUX_BC) {
      for(int j = 0; j < mesh.M; j++){
       C1_bc_type[bc_idx][j] = bc_type;
      }
    }
    else {
      cout << "unknown boundary condition type: " << bc_type << endl;
      cout << "defa ulting to electrode boundary condition";
      for(int j = 0; j < mesh.M; j++){
        C1_bc_type[bc_idx][j] = NO_FLUX_BC;
      }
    }
    // currently C2 doesn't have any patterning so every cell has either
    // a no flux or dirichlet b.c.
    bc_type = C2_bcs[bc_idx];
    for(int j = 0; j < mesh.M; j++){
      C2_bc_type[bc_idx][j] = bc_type;
    }
  }

  // arrays for determining the bc type of current x-location
  for (int bc_idx = 0; bc_idx < 2; bc_idx++) {
    int bc_type = Ey_bcs[bc_idx];
    if (bc_type == WALL_CHARGE_BC) {
      for (int i=0; i<mesh.N; i++) {
        Ey_bc_type[i][bc_idx] = bc_type;
      }
    }
  }

}

void IonTransportEqns2D::updateBCs(void) {
  //LHS and RHS BCs
  /*for (int j=0; j<mesh.M; j++) {
    // dirichlet right BC for C+ (C1)
    if (C1_bc_type[RIGHT_BC][j] == DIRICHLET_BC) {
      f1star_flux_sX[mesh.N_s-1][j] = -D1*mesh.dzdx_sX[mesh.N_s-1]*(2*(C1_RHS_BC_sX-C1_star[mesh.N_s-1][j])/mesh.dz)
                                  -D1*C1_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(-2*phi[mesh.N-1][j])/mesh.dz;
    }  
  }*/
}

double IonTransportEqns2D::frand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f *(fMax-fMin);
}

void IonTransportEqns2D::perturbOneConcentration(array2<double> &data) {
  srand(1);
  double onePercentOfLocalValue;
  array1<double>::opt netPerturbation(mesh.N);
  array2<double> localPerturbation(mesh.N,mesh.M);
  for (int i=0; i<mesh.N; i++) {
    for (int j=0; j<mesh.M; j++) {
      onePercentOfLocalValue = 0.00001*data[i][j];
      localPerturbation[i][j] = frand(-onePercentOfLocalValue, onePercentOfLocalValue);
      netPerturbation[i] = netPerturbation[i] + localPerturbation[i][j];
    }
  }
  for (int i=0; i<mesh.N; i++) {
    for (int j=0; j<mesh.M; j++) {
      data[i][j] = data[i][j] + localPerturbation[i][j] - netPerturbation[i]/mesh.M; 
    }
  }

}

void IonTransportEqns2D::printOneConcentration(string type) {
  if (type == "C1") {
    cout << "C1: " << endl;
    for (int i=0; i<mesh.N; i++) {
      for (int j=0; j<mesh.M; j++) {
        cout << setprecision(15) << setw(19) << C1_n[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "C2") {
    cout << "C2: " << endl;
    for (int i=0; i<mesh.N; i++) {
      for (int j=0; j<mesh.M; j++) {
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
  for(int n = 0; n < mesh.N; n++){
    for(int m = 0; m < mesh.M; m++){
      C1_star[n][m] = C1_n[n][m];
      C2_star[n][m] = C2_n[n][m];
    }
  }
}

void IonTransportEqns2D::updateConvergedValues(void) {
  for(int n = 0; n < mesh.N; n++){
    for(int m = 0; m < mesh.M; m++){
      C1_nMinus1[n][m] = C1_n[n][m];
      C2_nMinus1[n][m] = C2_n[n][m];
      C1_n[n][m] = C1_star[n][m];
      C2_n[n][m] = C2_star[n][m];
    }
  }
}

