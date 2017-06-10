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
    int numNeighbors = 4;
    int hsize = 1;
    array1<int>::opt neighbor;
    MPI_Comm comm;
    MPI_Wrapper();

    void getNeighbors(void);

};

#endif /*MPIWRAPPER_HPP*/
