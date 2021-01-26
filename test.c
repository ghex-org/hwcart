#include "hwcart.h"
#include <stdio.h>
#include <stdlib.h>

#define NSPLITS 5
int main(int argc, char *argv[])
{
    int domain[NSPLITS] = {
                               HWCART_MD_CORE,
                               HWCART_MD_L3CACHE,
                               HWCART_MD_NUMA,
                               HWCART_MD_SOCKET,
                               HWCART_MD_NODE
    };

    int topo[3*NSPLITS] = {
                           2, 1, 1, // core grid
                           2, 2, 1, // l3cache grid
                           1, 2, 2, // numa grid
                           1, 1, 2, // socket grid
                           2, 1, 1  // compute node grid
    };

    int ierr, comm_rank, comm_size, new_rank;
    int order = HWCartOrderXYZ;
    MPI_Comm hwcart_comm;

    if(argc>1){
      order = atoi(argv[1]);
    }

    MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    
    if(!hwcart_create(MPI_COMM_WORLD, NSPLITS, domain, topo, order, &hwcart_comm)){
        hwcart_print_rank_topology(hwcart_comm, NSPLITS, domain, topo, order);
        hwcart_free(&hwcart_comm);
    }

    MPI_Finalize();
}
