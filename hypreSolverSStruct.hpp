#ifndef HYPRESOLVERSSTRUCT_HPP
#define HYPRESOLVERSSTRUCT_HPP

#include "Array.hpp"
#include "mesh.hpp"
#include "mpiWrapper.hpp"
#include "HYPRE_sstruct_ls.h"

using namespace Array;
using namespace std;

class HypreSolverSStruct {

  public:
   HYPRE_SStructGrid     grid;
   HYPRE_SStructGraph    graph;
   HYPRE_SStructStencil  stencil_c1;
   HYPRE_SStructStencil  stencil_c2;
   HYPRE_SStructStencil  stencil_phi;
   HYPRE_SStructMatrix   A;
   HYPRE_SStructVector   b;
   HYPRE_SStructVector   x;
   HYPRE_SStructSolver   solver;
   HYPRE_SStructSolver   precond;
   HYPRE_Solver          par_solver;
   HYPRE_Solver          par_precond;
   int ilower[2], iupper[2];
   int bc_ilower[2];
   int bc_iupper[2];
   int ndim;
   int object_type;
   int solver_id = 0;

   void IonSystemSStructInit_Matrix(int ndim);
   void IonSystemSStruct_C1(array2<double> phi);
   void IonSystemSStruct_Matrix(double epsilon,
                                array2<double> C1,
                                array2<double> C2,
                                array2<double> phi);
   void IonSystemSStruct_Gauss(double epsilon);
   void IonSystemSStruct_RHS(array2<double> RHS_C1_star,
                             array2<double> RHS_C2_star,
                             array2<double> RHS_phi_star);
  void IonSystemSStruct_Solve(array2<double> C1,
                              array2<double> C2,
                              array2<double> phi);
   
  HypreSolverSStruct(Mesh& ref, MPI_Wrapper& ref2) : mesh(ref),mpi(ref2) {
    this-> mesh = mesh; 
    this->mpi = mpi;
    ilower[0] = mpi.pi*mesh.n;
    ilower[1] = mpi.pj*mesh.m;
    iupper[0] = ilower[0] + mesh.n-1;
    iupper[1] = ilower[1] + mesh.m-1;
  }
  ~HypreSolverSStruct(); 
  private:
    Mesh& mesh;
    MPI_Wrapper& mpi;
};

#endif /*HYPRESOLVERSSTRUCT_HPP*/
