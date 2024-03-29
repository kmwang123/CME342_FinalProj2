#include <stdlib.h>
#include <string>
#include <iomanip>
#include <iostream> 
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "Array.hpp"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "ionTransportEqns.hpp"

#define LEFT 0
#define RIGHT 1
#define SOUTH 2
#define NORTH 3

#define NO_FLUX_BC 0
#define DIRICHLET_BC 1
#define WALL_CHARGE_BC 2
#define LEFT_BC 0
#define RIGHT_BC 1
#define SOUTH_BC 0
#define NORTH_BC 1

using namespace std;

void IonTransportEqns2D::setUp(bool restart, bool perturb) {
  C1_star.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize); //include halo 
  C1_n.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C1_nMinus1.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_star.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_n.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  C2_nMinus1.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);
  phi.Allocate(mesh.n+2*mpi.hsize,mesh.m+2*mpi.hsize);

  istartc = mpi.hsize; 
  iendc = istartc + mesh.n-1;
  jstartc = mpi.hsize;
  jendc = jstartc + mesh.m-1;

  f1star_flux_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  f2star_flux_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  g1star_flux_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));
  g2star_flux_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));

  Ex_star_sX.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  Ey_star_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));

  RHS_C1_star.Allocate(mesh.n,mesh.m);
  RHS_C2_star.Allocate(mesh.n,mesh.m);
  RHS_phi_star.Allocate(mesh.n,mesh.m);

  phiM_sX_cY.Allocate(mesh.n_s+2*(mpi.hsize-1),mesh.m);
  C1star_sX_cY.Allocate(mesh.n_s+2*(mpi.hsize-1)+1,mesh.m);
  C2star_sX_cY.Allocate(mesh.n_s+2*(mpi.hsize-1)+1,mesh.m);

  phiM_cX_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));
  C1star_cX_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));
  C2star_cX_sY.Allocate(mesh.n,mesh.m_s+2*(mpi.hsize-1));

  this->restart = restart;

  if (!restart) {
    for(int i=istartc; i<=iendc; i++) {
      for (int j=jstartc; j<=jendc; j++) {
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

void IonTransportEqns2D::createBCarrays(int C1_bcs[2], int C2_bcs[2], int Ey_bcs[2]) {
  C1_bc_type.Allocate(2,mesh.m);
  C2_bc_type.Allocate(2,mesh.m);
  Ey_bc_type.Allocate(mesh.n,2);

  // arrays for determining the bc type of the current y-location
  for(int bc_idx = 0; bc_idx < 2; bc_idx++){
    int bc_type = C1_bcs[bc_idx];
    if (bc_type == DIRICHLET_BC || bc_type == NO_FLUX_BC) {
      for(int j = 0; j < mesh.m; j++){
       C1_bc_type[bc_idx][j] = bc_type;
      }
    }
    else {
      cout << "unknown boundary condition type: " << bc_type << endl;
      cout << "defa ulting to electrode boundary condition";
      for(int j = 0; j < mesh.m; j++){
        C1_bc_type[bc_idx][j] = NO_FLUX_BC;
      }
    }
    // currently C2 doesn't have any patterning so every cell has either
    // a no flux or dirichlet b.c.
    bc_type = C2_bcs[bc_idx];
    for(int j = 0; j < mesh.m; j++){
      C2_bc_type[bc_idx][j] = bc_type;
    }
  }

  // arrays for determining the bc type of current x-location
  for (int bc_idx = 0; bc_idx < 2; bc_idx++) {
    int bc_type = Ey_bcs[bc_idx];
    if (bc_type == WALL_CHARGE_BC) {
      for (int i=0; i<mesh.n; i++) {
        Ey_bc_type[i][bc_idx] = bc_type;
      }
    }
  }

}

void IonTransportEqns2D::updateBCs(void) {
  int pi, pj;
  pi = mpi.myid % mpi.p1;
  pj = mpi.myid / mpi.p1;
  //LHS and RHS BCs
  for (int j=0; j<mesh.m; j++) {
    if (pi == mpi.p1-1) {
      // dirichlet right BC for C+ (C1)
      if (C1_bc_type[RIGHT_BC][j] == DIRICHLET_BC) {
        f1star_flux_sX[mesh.n_s-1][j] = f1star_flux_sX[mesh.n_s-1][j]
                                        -D1*mesh.dzdx_sX[mesh.N_s-1]*(2*(C1_RHS_BC_sX-C1_star[iendc][jstartc+j])/mesh.dz); //diffusion
        f1star_flux_sX[mesh.n_s-1][j] = f1star_flux_sX[mesh.n_s-1][j]
                                        -D1*C1_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(2*(Phi_RHS_BC_sX-phi[iendc][jstartc+j])/mesh.dz); //em
                                        //-D1*C1_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(-2*phi[iendc][jstartc+j])/mesh.dz; //em
      }
      // dirichlet right BC for C- (C2)
      if (C2_bc_type[RIGHT_BC][j] == DIRICHLET_BC) {
        f2star_flux_sX[mesh.n_s-1][j] = f2star_flux_sX[mesh.n_s-1][j]
                                       -D2*mesh.dzdx_sX[mesh.N_s-1]*(2*(C2_RHS_BC_sX-C2_star[iendc][jstartc+j])/mesh.dz); //diffusion
        f2star_flux_sX[mesh.n_s-1][j] = f2star_flux_sX[mesh.n_s-1][j]
                                       +D2*C2_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(2*(Phi_RHS_BC_sX-phi[iendc][jstartc+j])/mesh.dz); //em
                                       //+D2*C2_RHS_BC_sX*mesh.dzdx_sX[mesh.N_s-1]*(-2*phi[iendc][jstartc+j])/mesh.dz; //em
      }
    }
    if (pi == 0) {
      // dirichlet left bc for C+ (C1)
      if (C1_bc_type[LEFT_BC][j] == DIRICHLET_BC) {
        f1star_flux_sX[0][j] = f1star_flux_sX[0][j]
                              -D1*mesh.dzdx_sX[0]*(2*(C1_star[istartc][jstartc+j]-C1_LHS_BC_sX)/mesh.dz); //diffusion
        f1star_flux_sX[0][j] = f1star_flux_sX[0][j]
                              -D1*C1_LHS_BC_sX*mesh.dzdx_sX[0]*(2*(phi[istartc][jstartc+j]-Phi_LHS_BC_sX)/mesh.dz); //em 
      }
      // no flux left bc for C- (C2)
      if (C2_bc_type[LEFT_BC][j] == NO_FLUX_BC) {
        f2star_flux_sX[0][j] = 0.0; 
        //if (mpi.myid == 0) {
        //  cout << "fluxes: " << setprecision(15) << setw(19) << f2star_flux_sX[0][j] << endl;
        //}
     } 
    } 
  }
  //North and South BCs
  for (int i=0; i<mesh.n; i++) {
    // no flux bc for south bc C+ (C1)
    if (pj == mpi.p2-1) {
      g1star_flux_sY[i][mesh.m_s-1] = 0.0;
      g2star_flux_sY[i][mesh.m_s-1] = 0.0;  
      if (Ey_bc_type[i][SOUTH_BC] == WALL_CHARGE_BC) {
        Ey_star_sY[i][mesh.m_s-1] = Ey_SBC_sX;
      }
    }
    if (pj == 0) {
      g1star_flux_sY[i][0] = 0.0;
      g2star_flux_sY[i][0] = 0.0;
      if (Ey_bc_type[i][NORTH_BC] == WALL_CHARGE_BC) {
        Ey_star_sY[i][0] = -Ey_NBC_sX;
      }
    }
  }
}

void IonTransportEqns2D::sendFluxes_updateRHS(int time_i, double dt) {
  MPI_Request request[2];
  MPI_Request send_req[2];
  MPI_Status status[2];
  double *halo_m = new double[5*mesh.m]();
  double *halo_n = new double[5*mesh.n]();

  //recv info if there is a neighbor to the left
  if (mpi.neighbor[LEFT] != -1) {
    MPI_Irecv(halo_m,5*mesh.m,MPI_DOUBLE,mpi.neighbor[LEFT],0,mpi.comm,&request[0]);

  }
  //otherwise, send info to neighbor on the right
  if (mpi.neighbor[RIGHT] != -1) {
    for (int j=0; j<5*mesh.m; j++) {
      if (j>=0 && j < mesh.m)
        halo_m[j] = f1star_flux_sX[mesh.n_s-1][j];
      if (j>= mesh.m && j <2*mesh.m)
        halo_m[j] = f2star_flux_sX[mesh.n_s-1][j-mesh.m];
      if (j>= 2*mesh.m && j<3*mesh.m)
        halo_m[j] = phiM_sX_cY[mesh.n_s-1][j-2*mesh.m];
      if (j>= 3*mesh.m && j<4*mesh.m)
        halo_m[j] = C1star_sX_cY[mesh.n_s-1][j-3*mesh.m];
      if (j>= 4*mesh.m && j<5*mesh.m)
        halo_m[j] = C2star_sX_cY[mesh.n_s-1][j-4*mesh.m];
    }
      MPI_Isend(halo_m,5*mesh.m,MPI_DOUBLE,mpi.neighbor[RIGHT],0,mpi.comm,&send_req[0]);
  }

  //recv info if there is a neighbor to the north
  if (mpi.neighbor[NORTH] != -1) {
    MPI_Irecv(halo_n,5*mesh.n,MPI_DOUBLE,mpi.neighbor[NORTH],0,mpi.comm,&request[1]);
  }
  //otherwise, send info to neighbor to the south
  if (mpi.neighbor[SOUTH] != -1) {
    for (int i=0; i<5*mesh.n; i++) {
      if (i>=0 && i< mesh.n)
        halo_n[i] = g1star_flux_sY[i][mesh.m_s-1];
      if (i>= mesh.n && i <2*mesh.n)
        halo_n[i] = g2star_flux_sY[i-mesh.n][mesh.m_s-1];
      if (i>=2*mesh.n && i<3*mesh.n)
        halo_n[i] = phiM_cX_sY[i-2*mesh.n][mesh.m_s-1];
      if (i>=3*mesh.n && i<4*mesh.n)
        halo_n[i] = C1star_cX_sY[i-3*mesh.n][mesh.m_s-1];
      if (i>=4*mesh.n && i<5*mesh.n)
        halo_n[i] = C2star_cX_sY[i-4*mesh.n][mesh.m_s-1];
    }
      MPI_Isend(halo_n,5*mesh.n,MPI_DOUBLE,mpi.neighbor[SOUTH],0,mpi.comm,&send_req[1]);
  }

  //update interior rhs
  updateInteriorRHS(time_i,dt); 

  //update L and R
  if (mpi.neighbor[RIGHT] != -1) 
    MPI_Wait(&send_req[0],&status[0]);
  if (mpi.neighbor[LEFT] != -1) {
    MPI_Wait(&request[0],&status[0]);
    //after recv'ed, unpack and store in flux matrices
    for (int j=0; j<5*mesh.m; j++) {
      if (j>=0 && j < mesh.m)
        f1star_flux_sX[0][j] = halo_m[j];
      if (j>= mesh.m && j <2*mesh.m)
        f2star_flux_sX[0][j-mesh.m] = halo_m[j];
      if (j>=2*mesh.m && j<3*mesh.m)
        phiM_sX_cY[0][j-2*mesh.m] = halo_m[j];
      if (j>= 3*mesh.m && j<4*mesh.m)
        C1star_sX_cY[0][j-3*mesh.m] = halo_m[j];
      if (j>= 4*mesh.m && j<5*mesh.m)
        C2star_sX_cY[0][j-4*mesh.m] = halo_m[j];
    }
    /*for (int j=0; j<mesh.m; j++) {
      for (int i=0; i<mesh.n_s; i++) {
      if (mpi.myid == 1) {
        cout << setprecision(15) << setw(19) << phiM_sX_cY[i][j] << " ";
      }
    }
    cout << endl;
  }*/
  }
  
  //update S and N
  if (mpi.neighbor[SOUTH] != -1) 
    MPI_Wait(&send_req[1],&status[1]);
  if (mpi.neighbor[NORTH] != -1) {
    MPI_Wait(&request[1],&status[1]);
    //after recv'ed, unpack and store in matrices
    for (int i=0; i<5*mesh.n; i++) {
      if (i>=0 && i< mesh.n)
        g1star_flux_sY[i][0] = halo_n[i];
      if (i>= mesh.n && i <2*mesh.n)
        g2star_flux_sY[i-mesh.n][0] = halo_n[i];
      if (i>=2*mesh.n && i<3*mesh.n)
        phiM_cX_sY[i-2*mesh.n][0] = halo_n[i];
      if (i>=3*mesh.n && i<4*mesh.n)
        C1star_cX_sY[i-3*mesh.n][0] = halo_n[i];
      if (i>=4*mesh.n && i<5*mesh.n)
        C2star_cX_sY[i-4*mesh.n][0] = halo_n[i];
    }

   /*for (int j=0; j<mesh.m_s; j++) {
    for (int i=0; i<mesh.n; i++) {
      if (mpi.myid == 2){
        cout << setprecision(15) << setw(19) << phiM_cX_sY[i][j] << " ";
    }}
    cout << endl;
    }*/ 
  }

  //update boundary rhs 
  updateBoundaryRHS(time_i,dt);
  delete [] halo_m;
  delete [] halo_n;
}

void IonTransportEqns2D::updateBoundaryRHS(int time_i, double dt) {
  // Add in the appropriate time term depending on which 
  // scheme is being used (1st order step 1, 2nd otherwise)
  if( time_i == 1 && !restart){
    for (int j=0; j<mesh.m; j++) {
      //north
      RHS_C1_star[0][j] = -(C1_star[istartc][jstartc+j] - C1_n[istartc][jstartc+j])/dt
                                - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n]*(f1star_flux_sX[1][j]-f1star_flux_sX[0][j])
                                - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[0][j+1]-g1star_flux_sY[0][j]);
      RHS_C2_star[0][j] = -(C2_star[istartc][jstartc+j] - C2_n[istartc][jstartc+j])/dt
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n]*(f2star_flux_sX[1][j]-f2star_flux_sX[0][j])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[0][j+1]-g2star_flux_sY[0][j]);
      
      //south
      RHS_C1_star[mesh.n-1][j] = -(C1_star[iendc][jstartc+j] - C1_n[iendc][jstartc+j])/dt
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+mesh.n-1]*(f1star_flux_sX[mesh.n_s-1][j]-f1star_flux_sX[mesh.n_s-2][j])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[mesh.n-1][j+1]-g1star_flux_sY[mesh.n-1][j]);
      RHS_C2_star[mesh.n-1][j] = -(C2_star[iendc][jstartc+j] - C2_n[iendc][jstartc+j])/dt
                                - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+mesh.n-1]*(f2star_flux_sX[mesh.n_s-1][j]-f2star_flux_sX[mesh.n_s-2][j])
                                - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[mesh.n-1][j+1]-g2star_flux_sY[mesh.n-1][j]);
    }
    for (int i=1; i<mesh.n-1; i++) {
      //left
      RHS_C1_star[i][0] = -(C1_star[istartc+i][jstartc] - C1_n[istartc+i][jstartc])/dt
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][0]-f1star_flux_sX[i][0])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m]*(g1star_flux_sY[i][1]-g1star_flux_sY[i][0]);
      RHS_C2_star[i][0] = -(C2_star[istartc+i][jstartc] - C2_n[istartc+i][jstartc])/dt
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][0]-f2star_flux_sX[i][0])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m]*(g2star_flux_sY[i][1]-g2star_flux_sY[i][0]);
      //right
      RHS_C1_star[i][mesh.m-1] = -(C1_star[istartc+i][jendc] - C1_n[istartc+i][jendc])/dt
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][mesh.m-1]-f1star_flux_sX[i][mesh.m-1])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+mesh.m-1]*(g1star_flux_sY[i][mesh.m_s-1]-g1star_flux_sY[i][mesh.m_s-2]);
      RHS_C2_star[i][mesh.m-1] = -(C2_star[istartc+i][jendc] - C2_n[istartc+i][jendc])/dt
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][mesh.m-1]-f2star_flux_sX[i][mesh.m-1])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+mesh.m-1]*(g2star_flux_sY[i][mesh.m_s-1]-g2star_flux_sY[i][mesh.m_s-2]);
                                  
    }
  }
  else {
    for (int j=0; j<mesh.m; j++) {
      //north
      RHS_C1_star[0][j] =  -(3*C1_star[istartc][jstartc+j]-4*C1_n[istartc][jstartc+j]+C1_nMinus1[istartc][jstartc+j])/(2*dt)
                           - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n]*(f1star_flux_sX[1][j]-f1star_flux_sX[0][j])
                           - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[0][j+1]-g1star_flux_sY[0][j]);
      RHS_C2_star[0][j] = -(3*C2_star[istartc][jstartc+j]-4*C2_n[istartc][jstartc+j]+C2_nMinus1[istartc][jstartc+j])/(2*dt)
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n]*(f2star_flux_sX[1][j]-f2star_flux_sX[0][j])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[0][j+1]-g2star_flux_sY[0][j]);

      //south
      RHS_C1_star[mesh.n-1][j] = -(3*C1_star[iendc][jstartc+j]-4*C1_n[iendc][jstartc+j]+C1_nMinus1[iendc][jstartc+j])/(2*dt)
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+mesh.n-1]*(f1star_flux_sX[mesh.n_s-1][j]-f1star_flux_sX[mesh.n_s-2][j])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[mesh.n-1][j+1]-g1star_flux_sY[mesh.n-1][j]);
      RHS_C2_star[mesh.n-1][j] = -(3*C2_star[iendc][jstartc+j]-4*C2_n[iendc][jstartc+j]+C2_nMinus1[iendc][jstartc+j])/(2*dt) 
                                - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+mesh.n-1]*(f2star_flux_sX[mesh.n_s-1][j]-f2star_flux_sX[mesh.n_s-2][j])
                                - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[mesh.n-1][j+1]-g2star_flux_sY[mesh.n-1][j]);
    }
    for (int i=1; i<mesh.n-1; i++) {
      //left
      RHS_C1_star[i][0] =  -(3*C1_star[istartc+i][jstartc]-4*C1_n[istartc+i][jstartc]+C1_nMinus1[istartc+i][jstartc])/(2*dt)
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][0]-f1star_flux_sX[i][0])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m]*(g1star_flux_sY[i][1]-g1star_flux_sY[i][0]);
      RHS_C2_star[i][0] = -(3*C2_star[istartc+i][jstartc]-4*C2_n[istartc+i][jstartc]+C2_nMinus1[istartc+i][jstartc])/(2*dt) 
                          - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][0]-f2star_flux_sX[i][0])
                          - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m]*(g2star_flux_sY[i][1]-g2star_flux_sY[i][0]);
      //right
      RHS_C1_star[i][mesh.m-1] = -(3*C1_star[istartc+i][jendc]-4*C1_n[istartc+i][jendc]+C1_nMinus1[istartc+i][jendc])/(2*dt) 
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][mesh.m-1]-f1star_flux_sX[i][mesh.m-1])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+mesh.m-1]*(g1star_flux_sY[i][mesh.m_s-1]-g1star_flux_sY[i][mesh.m_s-2]);
      RHS_C2_star[i][mesh.m-1] = -(3*C2_star[istartc+i][jendc]-4*C2_n[istartc+i][jendc]+C2_nMinus1[istartc+i][jendc])/(2*dt)    
                                 - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][mesh.m-1]-f2star_flux_sX[i][mesh.m-1])
                                 - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+mesh.m-1]*(g2star_flux_sY[i][mesh.m_s-1]-g2star_flux_sY[i][mesh.m_s-2]);

    }
  }

  /*if (mpi.myid ==0) {
    for (int j=0; j<mesh.m; j++) {
      for (int i=0; i<mesh.n; i++) {
        cout << setprecision(10) << setw(19)  << RHS_C2_star[i][j] << " ";
      }
      cout << endl;
    }
  }*/
}

void IonTransportEqns2D::updateInteriorRHS(int time_i, double dt) {
  // Add in the appropriate time term depending on which 
  // scheme is being used (1st order step 1, 2nd otherwise)
  if( time_i == 1 && !restart){  
    for (int i=1; i<mesh.n-1; i++) {
      for (int j=1; j<mesh.m-1; j++) {
        RHS_C1_star[i][j] = -(C1_star[istartc+i][jstartc+j] - C1_n[istartc+i][jstartc+j])/dt 
                            - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][j]-f1star_flux_sX[i][j]) 
                            - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[i][j+1]-g1star_flux_sY[i][j]);
        RHS_C2_star[i][j] = -(C2_star[istartc+i][jstartc+j] - C2_n[istartc+i][jstartc+j])/dt
                            - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][j]-f2star_flux_sX[i][j]) 
                            - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[i][j+1]-g2star_flux_sY[i][j]);
      }
    }
  }
  else {
    for (int i=1; i<mesh.n-1; i++) {
      for (int j=1; j<mesh.m-1; j++) {
        RHS_C1_star[i][j] = -(3*C1_star[istartc+i][jstartc+j]-4*C1_n[istartc+i][jstartc+j]+C1_nMinus1[istartc+i][jstartc+j])/(2*dt)
                            - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f1star_flux_sX[i+1][j]-f1star_flux_sX[i][j])
                            - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g1star_flux_sY[i][j+1]-g1star_flux_sY[i][j]);;
        RHS_C2_star[i][j] = -(3*C2_star[istartc+i][jstartc+j]-4*C2_n[istartc+i][jstartc+j]+C2_nMinus1[istartc+i][jstartc+j])/(2*dt)
                            - (1/mesh.dz)*mesh.dzdx_cX[mpi.pi*mesh.n+i]*(f2star_flux_sX[i+1][j]-f2star_flux_sX[i][j])
                            - (1/mesh.dn)*mesh.dndy_cY[mpi.pj*mesh.m+j]*(g2star_flux_sY[i][j+1]-g2star_flux_sY[i][j]);  
      }
    }
  }

  /*if (mpi.myid ==2) {
  for (int j=0; j<mesh.m; j++) {
    for (int i=0; i<mesh.n; i++) {
      cout << setprecision(10) << setw(19)  << RHS_C2_star[i][j] << " ";
    }
    cout << endl;
  }
  }*/


}

void IonTransportEqns2D::sendCenters_updateFluxes(array2<double> u_star, array2<double> v_star) {
  //int req_num = 0;
  MPI_Request request[2];
  MPI_Request send_req[2];
  MPI_Status status[2];
  double *haloc_m = new double[3*mesh.m]();
  double *haloc_n = new double[3*mesh.n]();
  //recv info if there is a neighbor to the right
  if (mpi.neighbor[RIGHT] != -1) {//(mpi.pi != mpi.p1-1) {
    //cout << "neighbor: " << mpi.neighbor[RIGHT] << endl;
    //recv information
    MPI_Irecv(haloc_m,3*mesh.m,MPI_DOUBLE,mpi.neighbor[RIGHT],0,mpi.comm,&request[0]);
    //if (mpi.myid == 2) {
    //  for (int j=0; j<3*mesh.m; j++) {
    //    cout << "recv: " << haloc_m[j] << endl;
    //  }
    //}
  }
  //otherwise, send info to neighbor on left
  if (mpi.neighbor[LEFT] != -1) {//(mpi.pi != 0)  {
    for (int j=0; j<3*mesh.m; j++) {
      if (j >= 0 && j < mesh.m)
        haloc_m[j] = phi[istartc][jstartc+j]; 
      if (j >= mesh.m && j <2*mesh.m)  
        haloc_m[j] = C1_star[istartc][jstartc+j-mesh.m];
      if(j >= 2*mesh.m && j < 3*mesh.m) 
        haloc_m[j] = C2_star[istartc][jstartc+j-2*mesh.m];
      //if (mpi.myid==3)
     // cout << "sending: " << haloc_m[j] << endl;

    }
    MPI_Isend(haloc_m,3*mesh.m,MPI_DOUBLE,mpi.neighbor[LEFT],0,mpi.comm,&send_req[0]);
  }

  //recv info if there is a neighbor to the south
  if (mpi.neighbor[SOUTH] != -1) {
    MPI_Irecv(haloc_n,3*mesh.n,MPI_DOUBLE,mpi.neighbor[SOUTH],0,mpi.comm,&request[1]);
  }
  //otherwise, send info to neighbor to the north
  if (mpi.neighbor[NORTH] != -1) {
    for (int i=0; i<3*mesh.n; i++) {
      if (i>=0 && i<mesh.n) 
        haloc_n[i] = phi[istartc+i][jstartc];
      if (i>=mesh.n && i<2*mesh.n) 
        haloc_n[i] = C1_star[istartc+i-mesh.n][jstartc];
      if (i>=2*mesh.n && i<3*mesh.n)
        haloc_n[i] = C2_star[istartc+i-2*mesh.n][jstartc];
    }
    MPI_Isend(haloc_n,3*mesh.n,MPI_DOUBLE,mpi.neighbor[NORTH],0,mpi.comm,&send_req[1]);
  }



  updateInteriorFluxes(u_star,v_star);

  //update L and R
  if (mpi.neighbor[LEFT] !=-1)
    MPI_Wait(&send_req[0],&status[0]);
  if (mpi.neighbor[RIGHT] != -1) { 
    MPI_Wait(&request[0],&status[0]);
    //after recv'ed, upack and store in matrices
    for (int j=0; j<3*mesh.m;j++) {
      if (j >= 0 && j < mesh.m)
        phi[iendc+1][jstartc+j] = haloc_m[j];
      if (j >= mesh.m && j <2*mesh.m)  
        C1_star[iendc+1][jstartc+j-mesh.m] = haloc_m[j];
      if(j >= 2*mesh.m && j < 3*mesh.m) 
        C2_star[iendc+1][jstartc+j-2*mesh.m] = haloc_m[j];
    }
    //calculate f_fluxes across iendc+1 and iendc
    
    
    //send these fluxes to right 
    
    /*for (int j=0; j<=jendc+1; j++) {
      for (int i=0; i<=iendc+1; i++) {
        if (mpi.myid == 0){
          cout << setprecision(15) << setw(19) << phi[i][j] << " ";
        }
      }
      cout << endl;
    }*/
  }

  //update S and N
  if (mpi.neighbor[NORTH] != -1)
    MPI_Wait(&send_req[1],&status[1]);
  if (mpi.neighbor[SOUTH] != -1) {
    MPI_Wait(&request[1],&status[1]);
    //after recv'ed, upack and store in matrices
    for (int i=0; i<3*mesh.n; i++) {
      if (i>=0 && i<mesh.n) 
        phi[istartc+i][jendc+1] = haloc_n[i];
      if (i>=mesh.n && i<2*mesh.n) 
        C1_star[istartc+i-mesh.n][jendc+1] = haloc_n[i];
      if (i>=2*mesh.n && i<3*mesh.n)
        C2_star[istartc+i-2*mesh.n][jendc+1] = haloc_n[i];
    }
  }

  //update boundary fluxes
  updateBoundaryFluxes(u_star,v_star);

   delete [] haloc_m;
   delete [] haloc_n;

}

void IonTransportEqns2D::updateBoundaryFluxes(array2<double> u_star, array2<double> v_star) {
   
  int iendX = mesh.n_s-1;
  int jendY = mesh.m_s-1;

  if (mpi.neighbor[RIGHT] != -1) {
    //x-dir
    for (int j=0; j<mesh.m; j++) {
      //e-field in x-dir
      Ex_star_sX[iendX][j] = -(1/mesh.dz)*mesh.dzdx_sX[mpi.pi*mesh.n+iendX]*(phi[istartc+iendX][jstartc+j]-phi[istartc+iendX-1][jstartc+j]);
      //flux of positive ions in x-dir
      f1star_flux_sX[iendX][j] = f1star_flux_sX[iendX][j]
                           -D1*(1/mesh.dz)*(mesh.dzdx_sX[mpi.pi*mesh.n+iendX])*(C1_star[istartc+iendX][jstartc+j]-C1_star[istartc+iendX-1][jstartc+j]); //diffusion 
      f1star_flux_sX[iendX][j] = f1star_flux_sX[iendX][j]
                            + 0.5*(C1_star[istartc+iendX][jstartc+j]+C1_star[istartc+iendX-1][jstartc+j])*(D1*Ex_star_sX[iendX][j]);//em
      f1star_flux_sX[iendX][j] = f1star_flux_sX[iendX][j]
                            + 0.5*(C1_star[istartc+iendX][jstartc+j]+C1_star[istartc+iendX-1][jstartc+j])*(D1*u_star[iendX][j]);//advection
      //flux of negative ions in y-dir
      f2star_flux_sX[iendX][j] = f2star_flux_sX[iendX][j]
                            -D2*(1/mesh.dz)*(mesh.dz)*(mesh.dzdx_sX[mpi.pi*mesh.n+iendX])*(C2_star[istartc+iendX][jstartc+j]-C2_star[istartc+iendX-1][jstartc+j]); //diffusion
      f2star_flux_sX[iendX][j] = f2star_flux_sX[iendX][j]
                            - 0.5*(C2_star[istartc+iendX][jstartc+j]+C2_star[istartc+iendX-1][jstartc+j])*(D2*Ex_star_sX[iendX][j]);//em
      f2star_flux_sX[iendX][j] = f2star_flux_sX[iendX][j]
                            - 0.5*(C2_star[istartc+iendX][jstartc+j]+C2_star[istartc+iendX-1][jstartc+j])*(D2*u_star[iendX][j]);//advection
      phiM_sX_cY[iendX][j] = 0.5*(phi[istartc+iendX][jstartc+j] - phi[istartc+iendX-1][jstartc+j]);
      C1star_sX_cY[iendX][j] = 0.5*(C1_star[istartc+iendX][jstartc+j] + C1_star[istartc+iendX-1][jstartc+j]);
      C2star_sX_cY[iendX][j] = 0.5*(C2_star[istartc+iendX][jstartc+j] + C2_star[istartc+iendX-1][jstartc+j]);
    }
  } 


  if (mpi.neighbor[SOUTH] != -1) {
    //y-dir
    for (int i=0; i<mesh.n; i++) {
      //electric field in y-dir
      Ey_star_sY[i][jendY] = -(1/mesh.dn)*mesh.dndy_sY[mpi.pj*mesh.m+jendY]*(phi[istartc+i][jstartc+jendY]-phi[istartc+i][jstartc+jendY-1]);
      //flux of positive ions in y-dir
      g1star_flux_sY[i][jendY] = g1star_flux_sY[i][jendY]
                             -D1*(1/mesh.dn)*mesh.dndy_sY[mpi.pj*mesh.m+jendY]*(C1_star[istartc+i][jstartc+jendY]-C1_star[istartc+i][jstartc+jendY-1]); //diffusion
      g1star_flux_sY[i][jendY] = g1star_flux_sY[i][jendY]
                            + 0.5*(C1_star[istartc+i][jstartc+jendY]+C1_star[istartc+i][jstartc+jendY-1])*(D1*Ey_star_sY[i][jendY]);//em
      g1star_flux_sY[i][jendY] = g1star_flux_sY[i][jendY]
                            + 0.5*(C1_star[istartc+i][jstartc+jendY]+C1_star[istartc+i][jstartc+jendY-1])*(D1*v_star[i][jendY]);//advection 
      //flux of negative ions in y-dir
      g2star_flux_sY[i][jendY] = g2star_flux_sY[i][jendY]
                           -D2*(1/mesh.dn)*mesh.dndy_sY[mpi.pj*mesh.m+jendY]*(C2_star[istartc+i][jstartc+jendY]-C2_star[istartc+i][jstartc+jendY-1]); //diffusion
      g2star_flux_sY[i][jendY] = g2star_flux_sY[i][jendY]
                           - 0.5*(C2_star[istartc+i][jstartc+jendY]+C2_star[istartc+i][jstartc+jendY-1])*(D2*Ey_star_sY[i][jendY]);//em
      g2star_flux_sY[i][jendY] = g2star_flux_sY[i][jendY]
                           - 0.5*(C2_star[istartc+i][jstartc+jendY]+C2_star[istartc+i][jstartc+jendY-1])*(D2*v_star[i][jendY]);//advection 

      phiM_cX_sY[i][jendY] = 0.5*(phi[istartc+i][jstartc+jendY] - phi[istartc+i][jstartc+jendY-1]);
      C1star_cX_sY[i][jendY] = 0.5*(C1_star[istartc+i][jstartc+jendY] + C1_star[istartc+i][jstartc+jendY-1]);
      C2star_cX_sY[i][jendY] = 0.5*(C2_star[istartc+i][jstartc+jendY] + C2_star[istartc+i][jstartc+jendY-1]); 
    }
  }
  /*for (int j=0; j<mesh.m_s; j++) {
    for (int i=0; i<mesh.n; i++) {
      if (mpi.myid == 2){
        cout << setprecision(15) << setw(19) << g2star_flux_sY[i][j] << " ";
    }}
    cout << endl;
  }*/
  /*for (int j=0; j<mesh.m; j++) { 
    for (int i=0; i<mesh.n_s; i++) {
     if (mpi.myid == 1) {
        cout << setprecision(15) << setw(19) << phiM_sX_cY[i][j] << " ";
      }
    }
    cout << endl;
  }*/

}

void IonTransportEqns2D::updateInteriorFluxes(array2<double> u_star, array2<double> v_star) { 
  //if (mpi.myid == 0)
  //cout << "interior fluxes: " << endl;

  int pi, pj;
  pi = mpi.myid % mpi.p1;
  pj = mpi.myid / mpi.p1;

  for (int i=1; i<mesh.n_s-1; i++) {
    for (int j=0; j<mesh.m; j++) {
 
      //electric field in x-dir
      Ex_star_sX[i][j] = -(1/mesh.dz)*mesh.dzdx_sX[pi*mesh.n+i]*(phi[istartc+i][jstartc+j]-phi[istartc+i-1][jstartc+j]); 

      //flux of positive ions in x-dir
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            -D1*(1/mesh.dz)*(mesh.dzdx_sX[pi*mesh.n+i])*(C1_star[istartc+i][jstartc+j]-C1_star[istartc+i-1][jstartc+j]); //diffusion
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i-1][jstartc+j])*(D1*Ex_star_sX[i][j]);//em
      f1star_flux_sX[i][j] = f1star_flux_sX[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i-1][jstartc+j])*(D1*u_star[i][j]);//advection
      //flux of negative ions in y-dir
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            -D2*(1/mesh.dz)*(mesh.dz)*(mesh.dzdx_sX[pi*mesh.n+i])*(C2_star[istartc+i][jstartc+j]-C2_star[istartc+i-1][jstartc+j]); //diffusion
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i-1][jstartc+j])*(D2*Ex_star_sX[i][j]);//em
      f2star_flux_sX[i][j] = f2star_flux_sX[i][j]
                            - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i-1][jstartc+j])*(D2*u_star[i][j]);//advection
    
      phiM_sX_cY[i][j] = 0.5*(phi[istartc+i][jstartc+j] - phi[istartc+i-1][jstartc+j]);
      C1star_sX_cY[i][j] = 0.5*(C1_star[istartc+i][jstartc+j] + C1_star[istartc+i-1][jstartc+j]);
      C2star_sX_cY[i][j] = 0.5*(C2_star[istartc+i][jstartc+j] + C2_star[istartc+i-1][jstartc+j]);
    }
  }

  for (int i=0; i<mesh.n; i++) {
    for (int j=1; j<mesh.m_s-1; j++) {
      //electric field in y-dir
      Ey_star_sY[i][j] = -(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(phi[istartc+i][jstartc+j]-phi[istartc+i][jstartc+j-1]); 

      //flux of positive ions in y-dir
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                             -D1*(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(C1_star[istartc+i][jstartc+j]-C1_star[istartc+i][jstartc+j-1]); //diffusion
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i][jstartc+j-1])*(D1*Ey_star_sY[i][j]);//em
      g1star_flux_sY[i][j] = g1star_flux_sY[i][j]
                            + 0.5*(C1_star[istartc+i][jstartc+j]+C1_star[istartc+i][jstartc+j-1])*(D1*v_star[i][j]);//advection 
 
      //flux of negative ions in y-dir
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           -D2*(1/mesh.dn)*mesh.dndy_sY[pj*mesh.m+j]*(C2_star[istartc+i][jstartc+j]-C2_star[istartc+i][jstartc+j-1]); //diffusion
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i][jstartc+j-1])*(D2*Ey_star_sY[i][j]);//em
      g2star_flux_sY[i][j] = g2star_flux_sY[i][j]
                           - 0.5*(C2_star[istartc+i][jstartc+j]+C2_star[istartc+i][jstartc+j-1])*(D2*v_star[i][j]);//advection 
 
      phiM_cX_sY[i][j] = 0.5*(phi[istartc+i][jstartc+j] - phi[istartc+i][jstartc+j-1]);
      C1star_cX_sY[i][j] = 0.5*(C1_star[istartc+i][jstartc+j] + C1_star[istartc+i][jstartc+j-1]);
      C2star_cX_sY[i][j] = 0.5*(C2_star[istartc+i][jstartc+j] + C2_star[istartc+i][jstartc+j-1]);
    }
  }

  /*for (int j=0; j<mesh.m_s; j++) {
    for (int i=0; i<mesh.n; i++) {
      if (mpi.myid == 2){
        cout << setprecision(15) << setw(19) << g2star_flux_sY[i][j] << " ";
    }}
    cout << endl;
  }*/
 
  /*for (int j=0; j<mesh.m; j++) { 
    for (int i=0; i<mesh.n_s; i++) {
     if (mpi.myid == 3) {
        cout << setprecision(15) << setw(19) << Ex_star_sX[i][j] << " ";
      }
    }
    cout << endl;
  }*/

}



double IonTransportEqns2D::frand(double fMin, double fMax) {
  double f = (double)rand() / RAND_MAX;
  return fMin + f *(fMax-fMin);
}

void IonTransportEqns2D::perturbOneConcentration(array2<double> &data) {
  srand(1);
  double onePercentOfLocalValue;
  array1<double>::opt netPerturbation(mesh.n);
  array2<double> localPerturbation(mesh.n,mesh.m);
  for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      onePercentOfLocalValue = 0.00001*data[i][j];
      localPerturbation[i][j] = frand(-onePercentOfLocalValue, onePercentOfLocalValue);
      netPerturbation[i] = netPerturbation[i] + localPerturbation[i][j];
    }
  }
  for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      data[i][j] = data[i][j] + localPerturbation[i][j] - netPerturbation[i]/mesh.M; 
    }
  }

}

void IonTransportEqns2D::printOneConcentration(string type) {
  if (type == "C1") {
    cout << "C1: " << endl;
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
        cout << setprecision(15) << setw(19) << C1_n[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "C2") {
    cout << "C2: " << endl;
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
        cout << setprecision(15) << setw(19) << C2_n[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "C1_star") {
    cout << "C1_star: " << endl;
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
        cout << setprecision(15) << setw(19) << C1_star[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "C2_star") {
    cout << "C2_star: " << endl;
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
        cout << setprecision(15) << setw(19) << C2_star[i][j] << " ";
      }
      cout << endl;
    }
  }
  else if (type == "phi") {
    cout << "phi: " << endl;
    for (int j=mpi.hsize; j<mesh.m+mpi.hsize; j++) {
      for (int i=mpi.hsize; i<mesh.n+mpi.hsize; i++) {
        cout << setprecision(15) << setw(19) << phi[i][j] << " ";
      }
      cout << endl;
    }
  }
  else { 
    cout << "Invalid type" << endl;
  }
}

void IonTransportEqns2D::setCstarValuesfrmCn(void) {
  for(int i = istartc; i <= iendc; i++){
    for(int j = jstartc; j <= jendc; j++){
      C1_star[i][j] = C1_n[i][j];
      C2_star[i][j] = C2_n[i][j];
    }
  }
}

void IonTransportEqns2D::updateConvergedValues(void) {
  for(int i = mpi.hsize; i < mesh.n+mpi.hsize; i++){
    for(int j = mpi.hsize; j < mesh.m+mpi.hsize; j++){
      C1_nMinus1[i][j] = C1_n[i][j];
      C2_nMinus1[i][j] = C2_n[i][j];
      C1_n[i][j] = C1_star[i][j];
      C2_n[i][j] = C2_star[i][j];
    }
  }
}

