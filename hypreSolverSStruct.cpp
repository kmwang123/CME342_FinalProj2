#include "hypreSolverSStruct.hpp"
#include "mesh.hpp"
#include "mpiWrapper.hpp"
#include "HYPRE_sstruct_ls.h"

using namespace Array;
using namespace std;

void HypreSolverSStruct::IonSystemSStructInit_Matrix(int ndim, bool restart) {
  this->ndim = ndim;
  this->restart = restart;
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
                                                 array2<double> phi,
                                                 double C1_LHS_BC_sX,
                                                 double C1_RHS_BC_sX,
                                                 double C2_RHS_BC_sX,
                                                 int time_i, double dt) {

  //Create an empty matrix object
  HYPRE_SStructMatrixCreate(mpi.comm, graph, &A);
  //Set the object type (by default HYPRE_SSTRUCT). This determines the
  //data structure used to store the matrix.  If you want to use
  //unstructured solvers, e.g. BoomerAMG, the object type should be
  //HYPRE_PARCSR. If the problem is purely structured (with one part), you
  //may want to use HYPRE_STRUCT to access the structured solvers.
  HYPRE_SStructMatrixSetObjectType(A, object_type);
  HYPRE_SStructMatrixInitialize(A);

  //initialize and compute shared data structures
  array2<double> phiM_sX_cY(mesh.n-1,mesh.m);
  array2<double> C1star_sX_cY(mesh.n-1,mesh.m); 
  array2<double> C2star_sX_cY(mesh.n-1,mesh.m);
  array2<double> phiM_cX_sY(mesh.n,mesh.m-1);
  array2<double> C1star_cX_sY(mesh.n,mesh.m-1);
  array2<double> C2star_cX_sY(mesh.n,mesh.m-1);

  istartc = mpi.hsize;
  iendc = mesh.n+mpi.hsize-1;
  jstartc = mpi.hsize;
  jendc = mesh.m+mpi.hsize-1;

  for (int i=0; i<mesh.n-1; i++) {
    for (int j=0; j<mesh.m; j++) {
      phiM_sX_cY[i][j] = 0.5*(phi[istartc+i+1][jstartc+j] - phi[istartc+i][jstartc+j]); 
      C1star_sX_cY[i][j] = 0.5*(C1[istartc+1+i][jstartc+j] + C1[istartc+i][jstartc+j]);
      C2star_sX_cY[i][j] = 0.5*(C2[istartc+1+i][jstartc+j] + C2[istartc+i][jstartc+j]);
      //if (mpi.myid == 0)
      //  cout << phiM_sX_cY[i][j] << " ";
    }
    //cout << endl;
  }

  for (int i=0; i<mesh.n; i++) {
    for (int j=0; j<mesh.m-1; j++) {
      phiM_cX_sY[i][j] = 0.5*(phi[istartc+i][jstartc+j+1] - phi[istartc+i][jstartc+j]);
      C1star_cX_sY[i][j] = 0.5*(C1[istartc+i][jstartc+j+1] + C1[istartc+i][jstartc+j]);
      C2star_cX_sY[i][j] = 0.5*(C2[istartc+i][jstartc+j+1] + C2[istartc+i][jstartc+j]);
      //if (mpi.myid == 0)
        //cout << phiM_cX_sY[i][j] << " ";
    }
    //cout << endl;
  }

  IonSystemSStruct_Gauss(epsilon);
  IonSystemSStruct_C1(phi,C1,phiM_sX_cY,C1star_sX_cY,phiM_cX_sY,C1star_cX_sY,C1_LHS_BC_sX,C1_RHS_BC_sX,time_i,dt);
  IonSystemSStruct_C2(phi,C2,phiM_sX_cY,C2star_sX_cY,phiM_cX_sY,C2star_cX_sY,C2_RHS_BC_sX,time_i,dt);
} 
void HypreSolverSStruct::IonSystemSStruct_C1(array2<double> phi,
                                             array2<double> C1, 
                                             array2<double> phiM_sX_cY,
                                             array2<double> C1star_sX_cY,
                                             array2<double> phiM_cX_sY,
                                             array2<double> C1star_cX_sY,
                                             double C1_LHS_BC_sX, 
                                             double C1_RHS_BC_sX,
                                             int time_i, double dt) {

  int part = 0;
  int var;
  //first set c1 stencil entries 
  int c1_indices[5] = {0, 1, 2, 3 ,4};
  int c1_phi_indices[5] = {5, 6, 7 , 8, 9};
  
  int nentries = 5;
  int nvalues = nentries*mesh.m*mesh.n;
  double *c1_values = new double[nvalues];
  double *phi_values = new double[nvalues];
  double *c1_phi_values = new double[nvalues];
  //set interior values 
  int elnumX, elnumY,ix,jy;
  for (int i = 0; i < nvalues; i += nentries) {
    ix = ((i/nentries)%mesh.n);
    jy = ((i/nentries)/mesh.n);
    elnumX = ((i/nentries)%mesh.n) + mpi.pi*mesh.n;
    elnumY = ((i/nentries)/mesh.n) + mpi.pj*mesh.m;
    //diffusion
    c1_values[i] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*(mesh.dzdx_sX[elnumX]+mesh.dzdx_sX[elnumX+1]) 
                  +(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(mesh.dndy_sY[elnumY]+mesh.dndy_sY[elnumY+1]);
    c1_values[i+1] =-(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX]; //west
    c1_values[i+2] =-(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX+1]; //east
    c1_values[i+3] =-(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY];//south
    c1_values[i+4] =-(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY+1];//north
    //em
    if (ix != 0 && ix != mesh.n-1) {
      phi_values[i] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*( phiM_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX]
                                                                  -phiM_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1] );
      c1_phi_values[i] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*(C1star_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX]
                                                                    +C1star_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1] );
    }
    if (jy != 0 && jy != mesh.m-1) {
      phi_values[i] = phi_values[i] + (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(phiM_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY] 
                                                                                     -phiM_cX_sY[ix][jy]  *mesh.dndy_sY[elnumY+1] );
      c1_phi_values[i] = c1_phi_values[i] + (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(C1star_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY]
                                                                                       +C1star_cX_sY[ix][jy]  *mesh.dndy_sY[elnumY+1] );
    }
   if (time_i == 1 && !restart) {
     phi_values[i] = phi_values[i] + 1./dt;
   }
   else {
     phi_values[i] = phi_values[i] + 3./(2*dt);
   }
    if (ix != 0) {
      phi_values[i+1] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*phiM_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX];//west
      c1_phi_values[i+1] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*C1star_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX];//west
    }
    if (ix != mesh.n-1) {
      phi_values[i+2] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*phiM_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1]; //east 
      c1_phi_values[i+2] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*C1star_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1]; //east 
    }
    if (jy != 0 && jy != mesh.m-1) {
      phi_values[i+3] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*phiM_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY];//south
      phi_values[i+4] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*phiM_cX_sY[ix][jy]*mesh.dndy_sY[elnumY+1]; //north
      c1_phi_values[i+3] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*C1star_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY];//south
      c1_phi_values[i+4] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*C1star_cX_sY[ix][jy]*mesh.dndy_sY[elnumY+1]; //north
    }
  }

  //set the c1-c1 connections
  var = 0;
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c1_indices, c1_values);
  //set the c1-c1 connections
  HYPRE_SStructMatrixAddToBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c1_indices, phi_values);  
  //set the c1-phi connections
  HYPRE_SStructMatrixAddToBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c1_phi_indices, c1_phi_values); 
  delete [] c1_values;
  delete [] phi_values;
  delete [] c1_phi_values;

  //incorporate boundary conditions
  nentries = 1;
  int nxvalues = nentries*mesh.n;
  int nyvalues = nentries*mesh.m;
  double *xvalues = new double[nxvalues]();
  double *yvalues = new double[nyvalues]();
  double *xphi_values = new double[nxvalues]();
  double *yphi_values = new double[nyvalues]();
  double *xphic1_values = new double[nxvalues]();
  double *yphic1_values = new double[nyvalues]();
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
    for (int i=0; i<nxvalues; i++) {
      //diff
      xvalues[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[1];
      xphi_values[i] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*phiM_cX_sY[i][0]*mesh.dndy_sY[1];
      xphic1_values[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[1]*C1star_cX_sY[i][0];
      //if (mpi.myid==0) 
        //cout << xphic1_values[i] << endl;
    }
    var = 0; 
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphic1_values);
    //need to change north values
    stencil_indices_bc[0] = 4; //4 means north values
    for (int i=0; i<nxvalues; i++) {
      //diff
      xvalues[i] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[1];
      //em
      xphic1_values[i] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[0]*mesh.dndy_sY[1]*C1star_cX_sY[i][0];
      //if (mpi.myid ==0)
      //  cout << xphic1_values[i] << endl;
    }
    var = 0;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphic1_values);
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
    for (int i=0; i<nxvalues; i++) {
      //diff
      xvalues[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1];
      //em
      xphi_values[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1]*phiM_cX_sY[i][mesh.m-2];
      xphic1_values[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1]*C1star_cX_sY[i][mesh.m-2];
      //if (mpi.myid ==4)
        //cout << xphic1_values[i] << endl;
    }
    var = 0;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphic1_values);
    //need to change south values
    stencil_indices_bc[0] = 3; //3 means south values
    for (int i=0; i<nxvalues; i++) {
      //diff
      xvalues[i] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1];
      //em
      xphi_values[i] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1]*phiM_cX_sY[i][mesh.m-2];
      xphic1_values[i] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[mesh.M-1]*mesh.dndy_sY[mesh.M-1]*C1star_cX_sY[i][mesh.m-2];
      if (mpi.myid == 4)
        cout << xphic1_values[i] << endl;
    }
    var = 0;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,xphic1_values);
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
    stencil_indices_bc[0] = 0; //center value
    for (int j=0; j<nyvalues; j++) {
      yvalues[j] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[0]*2
                  +(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[1];
      yphi_values[j] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*mesh.dzdx_sX[1]*phiM_sX_cY[0][j];
      yphic1_values[j] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[0]*(mesh.dzdx_sX[0]*2*C1_LHS_BC_sX+
                                                                mesh.dzdx_sX[1]*C1star_sX_cY[0][j]);
    }
    //use set to replace values
    var = 0;
    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yphic1_values);
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
    stencil_indices_bc[0] = 0; //center value
    for (int j=0; j<nyvalues; j++) {
      yvalues[j] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N-1]
                  +(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N_s-1]*2;
      yphi_values[j] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*mesh.dzdx_sX[mesh.N-1]*phiM_sX_cY[mesh.n-2][j];
      yphic1_values[j] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[mesh.N-1]*(mesh.dzdx_sX[mesh.N-1]*C1star_sX_cY[mesh.n-2][j]+
                                                                       mesh.dzdx_sX[mesh.N_s-1]*(2*C1_RHS_BC_sX));
      //if (mpi.myid==1) 
      //  cout << yphic1_values[j] << endl;
    }
    //use set to replace values
    var = 0;
    HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yvalues);
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                      var, nentries,
                                      stencil_indices_bc,yphi_values);
    var = 1;
    HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper,
                                    var, nentries,
                                    stencil_indices_bc,yphic1_values);
  }
  delete [] xvalues;
  delete [] yvalues;
  delete [] xphi_values;
  delete [] yphi_values;
  delete [] xphic1_values;
  delete [] yphic1_values;
  delete [] xclearval;
  delete [] yclearval;
}

void HypreSolverSStruct::IonSystemSStruct_C2(array2<double> phi,
                                             array2<double> C2,
                                             array2<double> phiM_sX_cY,
                                             array2<double> C2star_sX_cY,
                                             array2<double> phiM_cX_sY,
                                             array2<double> C2star_cX_sY,
                                             double C2_RHS_BC_sX,
                                             int time_i, double dt) {
  int part = 0;
  int var;
  //first set c2 stencil entries 
  int c2_indices[5] = {0, 1, 2, 3 ,4};
  int c2_phi_indices[5] = {5, 6, 7 , 8, 9};

  //set the c2-c2 connections
  var = 2;
  int nentries = 5;
  int nvalues = nentries*mesh.m*mesh.n;
  double *c2_values = new double[nvalues];
  double *phi_values = new double[nvalues];
  double *c2_phi_values = new double[nvalues];
  //set interior values 
  int elnumX, elnumY, ix, jy;
  for (int i = 0; i < nvalues; i += nentries) {
    ix = ((i/nentries)%mesh.n);
    jy = ((i/nentries)/mesh.n);
    elnumX = ((i/nentries)%mesh.n) + mpi.pi*mesh.n;
    elnumY = ((i/nentries)/mesh.n) + mpi.pj*mesh.m;
    //diffusion
    c2_values[i] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*(mesh.dzdx_sX[elnumX]+mesh.dzdx_sX[elnumX+1])
                  +(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(mesh.dndy_sY[elnumY]+mesh.dndy_sY[elnumY+1]);
    c2_values[i+1] =-(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX]; //west
    c2_values[i+2] =-(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*mesh.dzdx_sX[elnumX+1]; //east
    c2_values[i+3] =-(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY];//south
    c2_values[i+4] =-(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*mesh.dndy_sY[elnumY+1];//north
    //em
    if (ix != 0 && ix != mesh.n-1) {
      phi_values[i] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*( phiM_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX]
                                                                  -phiM_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1] );
      c2_phi_values[i] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*(C2star_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX]
                                                                    +C2star_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1] );
  
    }
    if (jy != 0 && jy != mesh.m-1) {
      phi_values[i] = phi_values[i] - (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(phiM_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY] 
                                                                                 -phiM_cX_sY[ix][jy]  *mesh.dndy_sY[elnumY+1] );
      c2_phi_values[i] = c2_phi_values[i] - (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*(C2star_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY]
                                                                                           +C2star_cX_sY[ix][jy]  *mesh.dndy_sY[elnumY+1] );
    }
    if (time_i == 1 && !restart) {
     phi_values[i] = phi_values[i] + 1./dt;
   }
   else {
     phi_values[i] = phi_values[i] + 3./(2*dt);
   }
   if (ix != 0) {
     phi_values[i+1] = -(1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*phiM_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX];//west
     c2_phi_values[i+1] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*C2star_sX_cY[ix-1][jy]*mesh.dzdx_sX[elnumX];//west
   }
   if (ix != mesh.n-1) {
      phi_values[i+2] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*phiM_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1]; //east 
      c2_phi_values[i+2] = (1/(mesh.dz*mesh.dz))*mesh.dzdx_cX[elnumX]*C2star_sX_cY[ix][jy]*mesh.dzdx_sX[elnumX+1]; //east 
   }
   if (jy != 0 && jy != mesh.m-1) {
     phi_values[i+3] = -(1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*phiM_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY];//south
     phi_values[i+4] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*phiM_cX_sY[ix][jy]*mesh.dndy_sY[elnumY+1]; //north
     c2_phi_values[i+3] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*C2star_cX_sY[ix][jy-1]*mesh.dndy_sY[elnumY];//south
     c2_phi_values[i+4] = (1/(mesh.dn*mesh.dn))*mesh.dndy_cY[elnumY]*C2star_cX_sY[ix][jy]*mesh.dndy_sY[elnumY+1]; //north
   }
  }
  //set the c2-c2 connections
  var = 2;
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c2_indices, c2_values);
  //set the c2-c2 connections
  HYPRE_SStructMatrixAddToBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c2_indices, phi_values);
  //set the c2-phi connections
  HYPRE_SStructMatrixAddToBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c2_phi_indices, c2_phi_values); 
  
  delete [] c2_values;
  delete [] phi_values;
  delete [] c2_phi_values;
  
  //incorporate boundary conditions
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
  }   


  delete [] xvalues;
  delete [] yvalues;
  delete [] xclearval;
  delete [] yclearval;
}


void HypreSolverSStruct::IonSystemSStruct_Gauss(double epsilon) {


  int part = 0;
  int var;
  //first set phi-stencil entries for Gauss's law
  int c1_phi_indices[1] = {5};
  int phi_indices[5] = {0, 1, 2, 3 ,4};
  int c2_phi_indices[1] = {6};

  //set the phi-c1 connections
  var = 2; //set values for phi connections
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
  var = 2; //set values for phi connections
  for (int i=0; i<nvalues; i++) {
    phi_c_values[i] = -1/(2*epsilon*epsilon);
  }
  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper,
                                  var, nentries,
                                  c2_phi_indices, phi_c_values);
  //set the phi-phi connections
  var = 2; //set values for phi connections
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
