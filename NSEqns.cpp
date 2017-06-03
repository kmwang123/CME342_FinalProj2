
#include "NSEqns.hpp"

using namespace std;

void NSEqns2D::setUp(void) {
  u_star.Allocate(mesh.N_s,mesh.M);
  u_n.Allocate(mesh.N_s,mesh.M);
  u_nMinus1.Allocate(mesh.N_s,mesh.M);
  v_star.Allocate(mesh.N,mesh.M_s);
  v_n.Allocate(mesh.N,mesh.M_s);
  v_nMinus1.Allocate(mesh.N,mesh.M_s);
}
