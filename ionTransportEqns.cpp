#include <stdlib.h>
#include <string>
#include <iomanip>
#include <iostream> 
#include <time.h>

#include "HYPRE_sstruct_ls.h"
#include "ionTransportEqns.hpp"

using namespace std;

void IonTransportEqns2D::setUp(bool restart, bool perturb) {
  C1_star.Allocate(mesh.N,mesh.M);
  C1_n.Allocate(mesh.N,mesh.M);
  C1_nMinus1.Allocate(mesh.N,mesh.M);
  C2_star.Allocate(mesh.N,mesh.M);
  C2_n.Allocate(mesh.N,mesh.M);
  C2_nMinus1.Allocate(mesh.N,mesh.M);
  phi.Allocate(mesh.N,mesh.M);

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

void IonTransportEqns2D::GaussLawSolve(void) {
  /*Setup HYPRE*/
  HYPRE_SStructGrid     grid;
  HYPRE_SStructGraph    graph;
  HYPRE_SStructStencil  stencil;
  HYPRE_SStructMatrix   A;
  HYPRE_SStructVector   b;
  HYPRE_SStructVector   x;

  /*Using Struct Solvers*/
  HYPRE_StructSolver solver;
  HYPRE_StructSolver precond;
  
  int object_type;
} 
