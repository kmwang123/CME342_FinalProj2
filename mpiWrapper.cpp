#include <mpi.h>
#include "mpiWrapper.hpp"

#define LEFT 0
#define RIGHT 1
#define SOUTH 2
#define NORTH 3

using namespace std;

MPI_Wrapper::MPI_Wrapper(void) {
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  comm = MPI_COMM_WORLD;
  Allocate(neighbor,numNeighbors);
}

void MPI_Wrapper::getNeighbors(void) {

  //set all to -1 first
  for (int i=0; i<numNeighbors; i++) {
    neighbor[i] = -1;
  }

  //figure out which neighbors are to the right
  if (pi != p1-1) {
    neighbor[RIGHT] = myid+1;
  }
  //figure out which neighbors are to the left
  if (pi != 0) {
    neighbor[LEFT] = myid-1;
  }
  //figure out which neighbors are to the north
  if (pj != p2-1) {
    neighbor[SOUTH] = myid+p1;
  }
  //figure out which neighbors are to the south
  if (pj != 0) {
    neighbor[NORTH] = myid-p1;
  }
  /*if (myid==2) {
    for (int i=0; i<numNeighbors; i++) {
      cout << "neighbors: " << neighbor[i] << endl;
    }
  }*/

}

void MPI_Wrapper::isend(array1<double>::opt temp, int count, string type, int dir, int req_num) {

  if (type == "double") {
   /* if (myid == 1) {
      for (int j=0; j<count; j++) {
        cout << "sending: " << temp[j] << endl;
      }
    }*/
    MPI_Isend(temp,count,MPI_DOUBLE,neighbor[dir],0,comm,&send_request[req_num]); 
    
  }
  else {
    cout << "invalid MPI_Isend type" << endl;
  }
}

void MPI_Wrapper::irecv(array1<double>::opt temp, int count, string type, int dir, int req_num) {

  if (type == "double") {
    MPI_Irecv(temp,count,MPI_DOUBLE,neighbor[dir],0,comm,&recv_request[req_num]);
  }
  else {
    cout << "invalid MPI_IRecv type" << endl;
  }
}

void MPI_Wrapper::waitall(string type) {
  if (type == "recv") {
    MPI_Wait(&recv_request[0],&status[0]);
    //MPI_Waitall(2,recv_request,status);
  }
  else if (type == "send") {
    MPI_Wait(&send_request[0],&status[0]);
    //MPI_Waitall(2,send_request,status);
  }
}
