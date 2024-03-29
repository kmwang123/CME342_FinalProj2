
#include "NSEqns.hpp"

using namespace std;

void NSEqns2D::setUp(bool restart) {
  u_star.Allocate(mesh.N_s,mesh.M);
  u_n.Allocate(mesh.N_s,mesh.M);
  u_nMinus1.Allocate(mesh.N_s,mesh.M);
  v_star.Allocate(mesh.N,mesh.M_s);
  v_n.Allocate(mesh.N,mesh.M_s);
  v_nMinus1.Allocate(mesh.N,mesh.M_s);

  //add this later
  /*if (restart) {
    
  }*/
}

void NSEqns2D::updateConvergedValues(void) {
  for(int n = 0; n < mesh.N; n++){
    for(int m = 0; m < mesh.M; m++){
      u_nMinus1[n][m] = u_n[n][m];
      v_nMinus1[n][m] = v_n[n][m];
      u_n[n][m] = u_star[n][m];
      v_n[n][m] = v_star[n][m];
    }
  }
  for(int m = 0 ; m < mesh.M; m++){
    u_nMinus1[mesh.N][m] = u_n[mesh.N][m];
    u_n[mesh.N][m] = u_star[mesh.N][m];
  }
  for(int n = 0 ; n <mesh.N; n++){
    v_nMinus1[n][mesh.M] = v_n[n][mesh.M];
    v_n[n][mesh.M] = v_star[n][mesh.M];
  }

}

void NSEqns2D::updateBCs(void) {
  /*LHS and RHS BCs*/
  for (int j=0; j<mesh.M; j++) {
    u_star[0][j] = 0.0;
    u_star[mesh.N_s-1][j] = 0.0;
  }
  /*North and South BCs*/
  for (int i=0; i<mesh.N; i++) {
    v_star[i][0] = 0.0;
    v_star[i][mesh.M_s-1] = 0.0;
  }
}
