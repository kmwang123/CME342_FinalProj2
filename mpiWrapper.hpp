#ifndef MPIWRAPPER_HPP
#define MPIWRAPPER_HPP

#include <mpi.h>
#include "Array.hpp"

using namespace Array;
using namespace std;
class MPI_Wrapper {

  public:
    int myid;
    int num_procs;
    int p1;
    int p2;
    int pi;
    int pj;
    MPI_Request send_request[2];
    MPI_Request recv_request[2];
    MPI_Status status[2];

    int numNeighbors = 4;
    int hsize = 1;
    array1<int>::opt neighbor;
    MPI_Comm comm;
    MPI_Wrapper();

    void getNeighbors(void);
    void isend(array1<double>::opt temp, int count, string type, int dir, int req_num);
    void irecv(array1<double>::opt temp, int count, string type, int dir, int req_num); 
    void waitall(string type);

};

#endif /*MPIWRAPPER_HPP*/
