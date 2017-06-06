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

void IonTransportEqns2D::GaussLawSolveStruct(void) {
  /*Setup HYPRE*/
  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   A;
  HYPRE_StructVector   b;
  HYPRE_StructVector   x;
  HYPRE_StructSolver   solver;
  HYPRE_StructSolver   precond;
  int num_iterations;
  double final_res_norm;
  
  /*Grid Setup*/
  HYPRE_StructGridCreate(mpi.comm, ndim, &grid);
  /*Set grid extents*/
  //hardcode everything in
  int pi = mpi.myid % mesh.N;
  int pj = mpi.myid / mesh.N;
  int ilower[2];
  int iupper[2];
  ilower[0] = pi*mesh.N;
  ilower[1] = pj*mesh.M;
  iupper[0] = ilower[0] + mesh.N-1;
  iupper[1] = ilower[1] + mesh.M-1;
  HYPRE_StructGridSetExtents(grid, ilower, iupper);
  /*Finalize grid assembly*/
  HYPRE_StructGridAssemble(grid); 

  /*Define Stencil*/
  int nentries = 5; //5 stencil entries per gridpoint
  HYPRE_StructStencilCreate(ndim, nentries, &stencil);
  int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
  /* Assign each of the 5 stencil entries */
  for (int entry = 0; entry < nentries; entry++)
    HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);

  /*Set up a Struct Matrix*/
  HYPRE_StructMatrixCreate(mpi.comm, grid, stencil, &A);
  /* Indicate that the matrix coefficients are ready to be set */
  HYPRE_StructMatrixInitialize(A);
  /*Set up matrix coefficients for stencil entries over all gridpoints*/
  int nvalues = mesh.N*mesh.M*nentries;
  double *values = new double[nvalues];
  int stencil_indices[nentries];
  //label stencil indices- these correspond to offsets defined above
  for (int j=0; j<nentries; j++) 
    stencil_indices[j] = j;
  //set interior values
  int elnum;
  for (int i = 0; i < nvalues; i += nentries) {
    elnum = i/nentries;
    values[i] = -((1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnum%mesh.N]*mesh.dzdx_sX[elnum%mesh.N] +
                  (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnum%mesh.N]*mesh.dzdx_sX[elnum%mesh.N+1]) 
                - ((1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnum/mesh.N]*mesh.dndy_sY[elnum/mesh.N] +
                  (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnum/mesh.N]*mesh.dndy_sY[elnum/mesh.N+1]);
    values[i+1] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnum%mesh.N]*mesh.dzdx_sX[elnum%mesh.N]; //west
    values[i+2] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnum%mesh.N]*mesh.dzdx_sX[elnum%mesh.N+1]; //east
    //cout << "values west: " << values[i+1] << endl;
    //cout << "values east: " << values[i+2] << endl;
    values[i+3] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnum/mesh.N]*mesh.dndy_sY[elnum/mesh.N]; //south
    //cout << "values south: " << values[i+3] << endl;
    values[i+4] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnum/mesh.N]*mesh.dndy_sY[elnum/mesh.N+1]; //north
    //cout << "values north: " << values[i+4] << endl;
  } 
  HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                 stencil_indices, values); 
  /*Set Boundary Conditions*/
  int bc_ilower[2];
  int bc_iupper[2];
  nentries = 1;
  int nxvalues = nentries*mesh.N;
  int nyvalues = nentries*mesh.M;
  double *xvalues = new double[nxvalues](); 
  double *yvalues = new double[nyvalues]();
  double *xclearval = new double[nxvalues]();
  double *yclearval = new double[nyvalues]();
  //set bc 
  for (int j = 0; j < nyvalues; j++)
    yclearval[j] = 0.0;
  for (int i=0; i<nxvalues; i++) 
    xclearval[i] = 0.0;
  int stencil_indices_bc[1];
  if (pj == 0) {
    //Bottom row of grid points 
    bc_ilower[0] = pi*mesh.N;
    bc_ilower[1] = 0;

    bc_iupper[0] = bc_ilower[0] + mesh.N-1;
    bc_iupper[1] = bc_ilower[1];

    stencil_indices_bc[0] = 3; //3 means south values
    HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices_bc, xclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[0];
    HYPRE_StructMatrixAddToBoxValues(A, bc_ilower, bc_iupper, nentries,
                                     stencil_indices_bc, xvalues);
  }
  if (pj == mpi.p2-1) {
    /* upper row of grid points */
    bc_ilower[0] = pi*mesh.N;
    bc_ilower[1] = pj*mesh.M + mesh.M-1;

    bc_iupper[0] = bc_ilower[0] + mesh.N-1;
    bc_iupper[1] = bc_ilower[1];

    stencil_indices_bc[0] = 4; //4 means north values

    HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                   stencil_indices_bc, xclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = +(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M_s-1];
    HYPRE_StructMatrixAddToBoxValues(A, bc_ilower, bc_iupper, nentries,
                                     stencil_indices_bc, xvalues);
  }
  if (pi == 0) {
    // Left row of grid points 
    bc_ilower[0] = 0;
    bc_ilower[1] = pj*mesh.M;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.M-1;

    stencil_indices_bc[0] = 1; //1 means west values

    HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                   stencil_indices_bc, yclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int j=0; j<nyvalues; j++)
      yvalues[j] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[0]; 
    HYPRE_StructMatrixAddToBoxValues(A, bc_ilower, bc_iupper, nentries,
                                     stencil_indices_bc, yvalues);
    
  }
  if (pi == mpi.p1-1) {
    // Right row of grid points 
    bc_ilower[0] = pi*mesh.N + mesh.N-1;
    bc_ilower[1] = pj*mesh.M;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.M-1;

    stencil_indices_bc[0] = 2;

    HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                   stencil_indices_bc, yclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int j=0; j<nyvalues; j++)
      yvalues[j] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N_s-1];
    HYPRE_StructMatrixAddToBoxValues(A, bc_ilower, bc_iupper, nentries,
                                     stencil_indices_bc, yvalues);

  }
  /*Finalized Matrix Assembly*/
  HYPRE_StructMatrixAssemble(A);
  //print matrix
  int ret = HYPRE_StructMatrixPrint("A.dat",A,0);
  
  /*Set up Sruct vectors b and x*/
  nvalues = mesh.M*mesh.N;
  double *values_vec = new double[nvalues];
  /* Create an empty vector object */
  HYPRE_StructVectorCreate(mpi.comm, grid, &b);
  HYPRE_StructVectorCreate(mpi.comm, grid, &x);
  /* Indicate that the vector coefficients are ready to be set */
  HYPRE_StructVectorInitialize(b);
  HYPRE_StructVectorInitialize(x);
  /*Set vector values*/
  for (int j=0; j<mesh.M; j++) {
    for (int i=0; i<mesh.N; i++) {
      values_vec[j*mesh.N+i] = (-1./(2*epsilon*epsilon))*(C1_n[i][j]-C2_n[i][j]);
    }
  }
  HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values_vec);
  /*Add BCs*/
  if (pj == 0) {
    // Bottom row of grid points 
    bc_ilower[0] = pi*mesh.N;
    bc_ilower[1] = 0;

    bc_iupper[0] = bc_ilower[0] + mesh.N-1;
    bc_iupper[1] = bc_ilower[1];
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = +mesh.dndy_sY[0]*Ey_SBC_sX/mesh.dn;
    
    HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, xvalues);

  }
  if (pj == mpi.p2-1) {
    //upper row of grid points
    bc_ilower[0] = pi*mesh.N;
    bc_ilower[1] = pj*mesh.M + mesh.M-1;

    bc_iupper[0] = bc_ilower[0] + mesh.N-1;
    bc_iupper[1] = bc_ilower[1];
    
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = +mesh.dndy_sY[mesh.M_s-1]*Ey_NBC_sX/mesh.dn;

    HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, xvalues);
  }
  if (pi == 0) {
    //Left row of grid points 
    bc_ilower[0] = 0;
    bc_ilower[1] = pj*mesh.M;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.M-1;

    for (int j=0; j < nyvalues; j++)
      yvalues[j] = -(2*Phi_LHS_BC_sX/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[0];
    HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, yvalues);
  }
  if (pi == mpi.p1-1) {
    // Right row of grid points 
    bc_ilower[0] = pi*mesh.N + mesh.N-1;
    bc_ilower[1] = pj*mesh.M;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.M-1;
    for (int j=0; j < nyvalues; j++)
      yvalues[j] = -(2*Phi_RHS_BC_sX/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N_s-1];
    HYPRE_StructVectorAddToBoxValues(b, bc_ilower, bc_iupper, yvalues);
  }
  //Debug
  int ret2 = HYPRE_StructVectorPrint("bvec", b, 0);
  /*Construct x-vector*/
  for (int i = 0; i < nvalues; i ++)
    values_vec[i] = 0.0;
  HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values_vec);

  delete [] xvalues;
  delete [] yvalues;
  delete [] xclearval;
  delete [] yclearval;
  delete [] values_vec;
  /*Finalize vector assembly*/
  HYPRE_StructVectorAssemble(b);
  HYPRE_StructVectorAssemble(x);
 
  if (solver_id_gauss == 0) {
    /*Set up and use a struct solver*/
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructPCGSetMaxIter(solver, 50 );
    HYPRE_StructPCGSetTol(solver, 1.0e-06 );
    HYPRE_StructPCGSetTwoNorm(solver, 1 );
    HYPRE_StructPCGSetRelChange(solver, 0 );
    HYPRE_StructPCGSetPrintLevel(solver, 2 ); /* print each CG iteration */
    HYPRE_StructPCGSetLogging(solver, 1);

    /* Use symmetric SMG as preconditioner */
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
    HYPRE_StructSMGSetMemoryUse(precond, 0);
    HYPRE_StructSMGSetMaxIter(precond, 1);
    HYPRE_StructSMGSetTol(precond, 0.0);
    HYPRE_StructSMGSetZeroGuess(precond);
    HYPRE_StructSMGSetNumPreRelax(precond, 1);
    HYPRE_StructSMGSetNumPostRelax(precond, 1);

    /* Set the preconditioner and solve */
    HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
                              HYPRE_StructSMGSetup, precond);
    HYPRE_StructPCGSetup(solver, A, b, x);
    HYPRE_StructPCGSolve(solver, A, b, x);

    /* Get some info on the run */
    HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
    HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    /* Clean up */
    HYPRE_StructPCGDestroy(solver);
  }
  if (solver_id_gauss == 1) {
    //Set default parameters
    int n_pre = 1;
    int n_post = 1;
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructSMGSetMemoryUse(solver, 0);
    HYPRE_StructSMGSetMaxIter(solver, 50);
    HYPRE_StructSMGSetTol(solver, 1.0e-06);
    HYPRE_StructSMGSetRelChange(solver, 0);
    HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
    HYPRE_StructSMGSetNumPostRelax(solver, n_post);
    /* Logging must be on to get iterations and residual norm info below */
    HYPRE_StructSMGSetLogging(solver, 1);

    /* Setup and solve */
    HYPRE_StructSMGSetup(solver, A, b, x);
    HYPRE_StructSMGSolve(solver, A, b, x);

    /* Get some info on the run */
    HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
    HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

    /* Clean up */
    HYPRE_StructSMGDestroy(solver);
  }
  if (solver_id_gauss == 2) {
     //Set default parameters
     int n_pre = 1;
     int n_post = 1;
     int skip = 0;
     int rap = 0;
     int relax = 1;    
       
     /* Options and setup */
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
      HYPRE_StructPFMGSetMaxIter(solver, 50);
      HYPRE_StructPFMGSetTol(solver, 1.0e-06);
      HYPRE_StructPFMGSetRelChange(solver, 0);
      HYPRE_StructPFMGSetRAPType(solver, rap);
      HYPRE_StructPFMGSetRelaxType(solver, relax);
      HYPRE_StructPFMGSetNumPreRelax(solver, n_pre);
      HYPRE_StructPFMGSetNumPostRelax(solver, n_post);
      HYPRE_StructPFMGSetSkipRelax(solver, skip);
      HYPRE_StructPFMGSetPrintLevel(solver, 1);
      HYPRE_StructPFMGSetLogging(solver, 1);
      HYPRE_StructPFMGSetup(solver, A, b, x);
      /* Solve */
      HYPRE_StructPFMGSolve(solver, A, b, x);
      /* Get info and release memory */
      HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      HYPRE_StructPFMGDestroy(solver);

  }
  //Debug
  int ret3 = HYPRE_StructVectorPrint("xvec", x, 0);
  //reshape and store in phi
  HYPRE_StructVectorGetBoxValues(x,ilower,iupper,values);
  for (int j=0; j<mesh.M; j++) {
    for (int i=0; i<mesh.N; i++) {
      phi[i][j] = values[j*mesh.N+i];
    }
  }
  for (int i=0; i<mesh.N; i++) {
    for (int j=0; j<mesh.M; j++) {
      cout << phi[i][j] << " ";
    }
    cout << endl;
  }

  /* Free memory */
  delete [] values;
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);
}

void IonTransportEqns2D::GaussLawSolveSStruct(void) {
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

  //hardcode everything in
  int ilower[2] = {0,0};
  int iupper[2] = {mesh.N,mesh.M};
} 
