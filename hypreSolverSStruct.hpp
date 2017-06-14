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
   int n_pre = 1;
   int n_post = 1;
   int istartc, iendc, jstartc, jendc;
   int solver_id = 0;
   bool restart;

   void IonSystemSStructInit_Matrix(int ndim, bool restart);
   void IonSystemSStruct_C1(array2<double> phi,
                            array2<double> C1,
                            array2<double> phiM_sX_cY,
                            array2<double> C1star_sX_cY,
                            array2<double> phiM_cX_sY,
                            array2<double> C1star_cX_sY,
                            double C1_LHS_BC_sX, 
                            double C1_RHS_BC_sX,
                            int time_i, double dt);
   void IonSystemSStruct_C2(array2<double> phi,
                            array2<double> C2,
                            array2<double> phiM_sX_cY,
                            array2<double> C2star_sX_cY,
                            array2<double> phiM_cX_sY,
                            array2<double> C2star_cX_sY,
                            double C2_RHS_BC_sX,
                            int time_i, double dt);
   void IonSystemSStruct_Matrix(double epsilon,
                                array2<double> C1,
                                array2<double> C2,
                                array2<double> phi,
                                double C1_LHS_BC_sX,
                                double C1_RHS_BC_sX,
                                double C2_RHS_BC_sX,
                                int time_i, double dt);
   void IonSystemSStruct_Gauss(double epsilon);
   void IonSystemSStruct_RHS(array2<double> RHS_C1_star,
                             array2<double> RHS_C2_star,
                             array2<double> RHS_phi_star);
  void IonSystemSStruct_Solve(void);
  
  HypreSolverSStruct(Mesh& ref, MPI_Wrapper& ref2) : mesh(ref),mpi(ref2) {
    this-> mesh = mesh; 
    this->mpi = mpi;
    ilower[0] = mpi.pi*mesh.n;
    ilower[1] = mpi.pj*mesh.m;
    iupper[0] = ilower[0] + mesh.n-1;
    iupper[1] = ilower[1] + mesh.m-1;
  }
  void CleanUp(void); 
  private:
    Mesh& mesh;
    MPI_Wrapper& mpi;
};

#endif /*HYPRESOLVERSSTRUCT_HPP*/
