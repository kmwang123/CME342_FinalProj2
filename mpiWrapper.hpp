#ifndef MPIWRAPPER_HPP
#define MPIWRAPPER_HPP

#include <mpi.h>

using namespace std;
class MPI_Wrapper {

  public:
    int myid;
    int num_procs;
    int p1;
    int p2;
    MPI_Comm comm;
    MPI_Wrapper();

};

#endif /*MPIWRAPPER_HPP*/
