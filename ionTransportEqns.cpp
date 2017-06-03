
#include "ionTransportEqns.hpp"

using namespace std;

void IonTransportEqns2D::setUp(void) {
    C1_star.Allocate(mesh.N,mesh.M);
    C1_n.Allocate(mesh.N,mesh.M);
    C1_nMinus1.Allocate(mesh.N,mesh.M);
    C2_star.Allocate(mesh.N,mesh.M);
    C2_n.Allocate(mesh.N,mesh.M);
    C2_nMinus1.Allocate(mesh.N,mesh.M);
    phi.Allocate(mesh.N,mesh.M);
}
