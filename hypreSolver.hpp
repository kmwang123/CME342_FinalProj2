#ifndef HYPRESOLVER_HPP
#define HYPRESOLVER_HPP

#include "Array.hpp"
#include "mesh.hpp"
#include "mpiWrapper.hpp"
#include "HYPRE_struct_ls.h"

using namespace Array;
using namespace std;

class HypreSolver {

  public:
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   A;
    HYPRE_StructVector   b;
    HYPRE_StructVector   x;
    HYPRE_StructVector bigphi;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    int ilower[2];
    int iupper[2];
    int bc_ilower[2];
    int bc_iupper[2];
    double *values;
    double *values_vec;
    double *xvalues;
    double *yvalues;
    double *xclearval;
    double *yclearval;
    int num_iterations;
    double final_res_norm;
    HypreSolver();
    void GaussLawSolveStruct_Matrix(Mesh mesh, 
                         MPI_Wrapper mpi,
                         int ndim);
    
    void GaussLawSolveStruct_RHS(Mesh mesh,
                         MPI_Wrapper mpi,
                         int ndim,
                         array2<double> C1,
                         array2<double> C2,
                         double epsilon,
                         double Ey_SBC_sX,
                         double Ey_NBC_sX,
                         double Phi_LHS_BC_sX,
                         double Phi_RHS_BC_sX);
    void updateRHSPhi(Mesh mesh,
                         MPI_Wrapper mpi,
                         int ndim,
                         array2<double> C1,
                         array2<double> C2,
                         array2<double> phi,
                         array2<double> RHS_phi_star,
                         double epsilon,
                         double Ey_SBC_sX,
                         double Ey_NBC_sX,
                         double Phi_LHS_BC_sX,
                         double Phi_RHS_BC_sX);
    void GaussLawSolveStruct_Solve(Mesh mesh, 
                         MPI_Wrapper mpi, 
                         int ndim, 
                         array2<double> phi, 
                         int solver_id);
    void CleanUp(void);
};


#endif /*HYPRESOLVER_HPP*/
