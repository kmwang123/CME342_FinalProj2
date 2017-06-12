#include "hypreSolverSStruct.hpp"
#include "mesh.hpp"
#include "mpiWrapper.hpp"
#include "HYPRE_sstruct_ls.h"

using namespace Array;
using namespace std;

void HypreSolverSStruct::IonSystemSStructInit_Matrix(int ndim) {
  this->ndim = ndim;
  int nparts = 1;
  int part = 0;
  int stencil_size;
  int var;
  int entry;
  //Create an empty 2D grid object
  HYPRE_SStructGridCreate(mpi.comm, ndim, nparts, &grid);
  //add a new box to the grid
  HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
  //Set variable type and number of variables on each part
  int nvars = 3;
  HYPRE_SStructVariable vartypes[3] = {HYPRE_SSTRUCT_VARIABLE_CELL,
                                       HYPRE_SSTRUCT_VARIABLE_CELL,
                                       HYPRE_SSTRUCT_VARIABLE_CELL};
  for (int i=0; i<nparts; i++) 
    HYPRE_SStructGridSetVariables(grid, i, nvars, vartypes);
  //finalize grid assembly
  HYPRE_SStructGridAssemble(grid);

  //Stencil object for variable C1 (labeled as variable 0)
  int offsets[10][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}, {0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
  stencil_size = 10;
  HYPRE_SStructStencilCreate(ndim, stencil_size, &stencil_c1);
  //the first 5 entries are for the c1-c1 connections
  var = 0; //connect to variable 0
  for (entry=0; entry < stencil_size-5; entry++) 
    HYPRE_SStructStencilSetEntry(stencil_c1, entry, offsets[entry], var);
  //the last 5 entries are for the c1-phi connections
  var = 1; //connect to variable 1
  for (entry=5; entry < stencil_size; entry++) 
    HYPRE_SStructStencilSetEntry(stencil_c1, entry, offsets[entry], var);

  //Stencil object for variable phi (labeled as variable 1)
  int offsets1[7][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}, {0,0}, {0,0}};
  stencil_size = 7;
  HYPRE_SStructStencilCreate(ndim, stencil_size, &stencil_phi);
  //first 5 entries are for the phi-phi connections
  var = 1; //connect to variable 1
  for (entry=0; entry < stencil_size-5; entry++) 
    HYPRE_SStructStencilSetEntry(stencil_phi, entry, offsets1[entry], var);
  var = 0; //connect to variable 0
  entry = 5;
  HYPRE_SStructStencilSetEntry(stencil_phi, entry, offsets1[entry], var);
  var = 2; //connect to variable 2
  entry = 6;
  HYPRE_SStructStencilSetEntry(stencil_phi, entry, offsets1[entry], var);

  //Stencil object for variable c2 (labeled as variable 2)
  int offsets2[10][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}, {0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
  stencil_size = 10;
  HYPRE_SStructStencilCreate(ndim, stencil_size, &stencil_c2);
  var = 1; //connect to variable 1 (phi)
  for (entry=0; entry < stencil_size-5; entry++) 
    HYPRE_SStructStencilSetEntry(stencil_c2, entry, offsets2[entry], var);
  var = 2; //connect to variable 2 (c2)
  for (entry=5; entry < stencil_size; entry++)  
    HYPRE_SStructStencilSetEntry(stencil_c2, entry, offsets2[entry], var);

  //Set up the Graph-this determines non-zero structure of the matrix
  HYPRE_SStructGraphCreate(mpi.comm, grid, &graph);
  if (solver_id == 0) {
    object_type = HYPRE_SSTRUCT;
  }
  HYPRE_SStructGraphSetObjectType(graph, object_type);

  //Assign the c1-stencil we created for variable c1 (var 0)
  var = 0;
  HYPRE_SStructGraphSetStencil(graph, part, var, stencil_c1);
  //Assign the phi-stencil we created for variable phi (var 1)
  var = 1;
  HYPRE_SStructGraphSetStencil(graph, part, var, stencil_phi); 
  //Assign the c2-stencil we created for variable c2 (var 2)
  var = 2;
  HYPRE_SStructGraphSetStencil(graph, part, var, stencil_c2); 

  //Assemble the graph
  HYPRE_SStructGraphAssemble(graph); 

}

void HypreSolverSStruct::IonSystemSStruct_Matrix(double epsilon,
                                                 array2<double> C1,
                                                 array2<double> C2,
                                                 array2<double> phi) {

  //Create an empty matrix object
  HYPRE_SStructMatrixCreate(mpi.comm, graph, &A);
  //Set the object type (by default HYPRE_SSTRUCT). This determines the
  //data structure used to store the matrix.  If you want to use
  //unstructured solvers, e.g. BoomerAMG, the object type should be
  //HYPRE_PARCSR. If the problem is purely structured (with one part), you
  //may want to use HYPRE_STRUCT to access the structured solvers.
  HYPRE_SStructMatrixSetObjectType(A, object_type);
  HYPRE_SStructMatrixInitialize(A);

  IonSystemSStruct_Gauss(epsilon);
  IonSystemSStruct_C1(phi);
} 
void HypreSolverSStruct::IonSystemSStruct_C1(array2<double> phi) {
  int part = 0;
  int var;
  //first set 
}


void HypreSolverSStruct::IonSystemSStruct_Gauss(double epsilon) {


  int part = 0;
  int var;
  //first set phi-stencil entries for Gauss's law
  int c1_phi_indices[1] = {5};
  int phi_indices[5] = {0, 1, 2, 3 ,4};
  int c2_phi_indices[1] = {6};

  //set the phi-c1 connections
  var = 0;
  int nentries = 1;
  int nvalues = nentries*mesh.m*mesh.n;
  double *phi_c_values = new double[nvalues];
  for (int i=0; i<nvalues; i++) {
    phi_c_values[i] = 1/(2*epsilon*epsilon);
  }
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c1_phi_indices, phi_c_values);
  //set the phi-c2 connections
  var = 2;
  for (int i=0; i<nvalues; i++) {
    phi_c_values[i] = -1/(2*epsilon*epsilon);
  }
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c2_phi_indices, phi_c_values);
  //set the phi-phi connections
  var = 1;
  nentries = 5;
  nvalues = nentries*mesh.m*mesh.n;
  double *phi_values = new double[nvalues];
  //set interior values
  int elnumX, elnumY;
  for (int i = 0; i < nvalues; i += nentries) {
    elnumX = ((i/nentries)%mesh.n) + mpi.pi*mesh.n;
    elnumY = ((i/nentries)/mesh.n) + mpi.pj*mesh.m;
    phi_values[i] = -((1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX] +
                  (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX+1])
                - ((1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY] +
                  (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY+1]);
    phi_values[i+1] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX]; //west
    phi_values[i+2] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX+1]; //east
    phi_values[i+3] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY]; //south
    phi_values[i+4] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY+1]; //north
  }
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  phi_indices, phi_values);
  delete [] phi_values;
  delete [] phi_c_values;


  //incorporate boundary conditions
  //for Gauss's law
  nentries = 1;
  int nxvalues = nentries*mesh.n;
  int nyvalues = nentries*mesh.m;
  double *xvalues = new double[nxvalues]();
  double *yvalues = new double[nyvalues]();
  double *xclearval = new double[nxvalues]();
  double *yclearval = new double[nyvalues]();
  int stencil_indices_bc[1];

  for (int j = 0; j < nyvalues; j++)
    yclearval[j] = 0.0;
  for (int i=0; i<nxvalues; i++)
    xclearval[i] = 0.0;

  if (mpi.pj == 0) {
    //Bottom row of grid points 
    bc_ilower[0] = mpi.pi*mesh.n;
    bc_ilower[1] = mpi.pj*mesh.m;

    bc_iupper[0] = bc_ilower[0] + mesh.n-1;
    bc_iupper[1] = bc_ilower[1];

    stencil_indices_bc[0] = 3; //3 means south values
    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, 
                                    var, nentries,
                                    stencil_indices_bc,xclearval);      
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[0];
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, 
                                    var, nentries,
                                    stencil_indices_bc,xvalues); 
  }
  if (mpi.pj == mpi.p2-1) {
    // upper row of grid points 
    bc_ilower[0] = mpi.pi*mesh.n;
    bc_ilower[1] = mpi.pj*mesh.m + mesh.m-1;

    bc_iupper[0] = bc_ilower[0] + mesh.n-1;
    bc_iupper[1] = bc_ilower[1];

    stencil_indices_bc[0] = 4; //4 means north values
    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int i=0; i<nxvalues; i++)
      xvalues[i] = +(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M_s-1];
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xvalues);
  }
  if (mpi.pi == 0) {
    // Left row of grid points
    bc_ilower[0] = mpi.pi*mesh.n;
    bc_ilower[1] = mpi.pj*mesh.m;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.m-1;

    stencil_indices_bc[0] = 1; //1 means west values
    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int j=0; j<nyvalues; j++)
      yvalues[j] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[0];
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yvalues);
  }
  if (mpi.pi == mpi.p1-1) {
    // Right row of grid points 
    bc_ilower[0] = mpi.pi*mesh.n + mesh.n-1;
    bc_ilower[1] = mpi.pj*mesh.m;

    bc_iupper[0] = bc_ilower[0];
    bc_iupper[1] = bc_ilower[1] + mesh.m-1;

    stencil_indices_bc[0] = 2;

    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yclearval);
    //need to change center value
    stencil_indices_bc[0] = 0; //0 means center values
    for (int j=0; j<nyvalues; j++)
      yvalues[j] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N_s-1];
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yvalues);
  }

  delete [] xvalues;
  delete [] yvalues;
  delete [] xclearval;
  delete [] yclearval;
}

void HypreSolverSStruct::IonSystemSStruct_RHS(array2<double> RHS_C1_star,
                                              array2<double> RHS_C2_star,
                                              array2<double> RHS_phi_star) {
  int nvalues = mesh.n*mesh.m;
  double *values = new double[nvalues];
  double *values1 = new double[nvalues];
  double *values2 = new double[nvalues];
  int part = 0;
  int var;
  
  //create an empty vector object
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
  //set object type
  HYPRE_SStructVectorSetObjectType(b, object_type);
  HYPRE_SStructVectorInitialize(b);
  //set values for b
  for (int j=0; j<mesh.m; j++) {
    for (int i=0; i<mesh.n; i++) {
      values[j*mesh.n+i] = RHS_C1_star[i][j];
      values1[j*mesh.n+i] = RHS_phi_star[i][j];
      values2[j*mesh.n+i] = RHS_C2_star[i][j];
    }
  }
 
  var = 0;
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, values);
  var = 1;
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, values1);
  var = 2;
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, values2);

  HYPRE_SStructVectorAssemble(b);
 
  delete [] values;
  delete [] values1;
  delete [] values2;
}

void HypreSolverSStruct::IonSystemSStruct_Solve(array2<double> C1,
                                                array2<double> C2,
                                                array2<double> phi) {

}

HypreSolverSStruct::~HypreSolverSStruct() {
  //free memory to prevent leakage
  HYPRE_SStructGridDestroy(grid);
  HYPRE_SStructStencilDestroy(stencil_phi);
  HYPRE_SStructStencilDestroy(stencil_c1);
  HYPRE_SStructStencilDestroy(stencil_c2);
  HYPRE_SStructGraphDestroy(graph);
   HYPRE_SStructMatrixDestroy(A);
   HYPRE_SStructVectorDestroy(b);
}