#ifndef HYPRESOLVER_HPP
#define HYPRESOLVER_HPP

#include "Array.hpp"
#include "mesh.hpp"
#include "mpiWrapper.hpp"

void GaussLawSolveStruct(Mesh mesh, 
                         MPI_Wrapper mpi,
                         int ndim,
                         array2<double> C1, 
                         array2<double> C2, 
                         array2<double> phi,
                         double epsilon, 
                         double Ey_SBC_sX, 
                         double Ey_NBC_sX, 
                         double Phi_LHS_BC_sX, 
                         double Phi_RHS_BC_sX, 
                         int solver_id);




#endif /*HYPRESOLVER_HPP*/
