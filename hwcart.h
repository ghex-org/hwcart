#ifndef _HWCART_H
#define _HWCART_H

#include <mpi.h>

#define HWCART_MD_HWTHREAD 1
#define HWCART_MD_CORE 2
#define HWCART_MD_L1CACHE 3
#define HWCART_MD_L2CACHE 4
#define HWCART_MD_L3CACHE 5
#define HWCART_MD_SOCKET 6
#define HWCART_MD_NUMA 7
#define HWCART_MD_BOARD 8
#define HWCART_MD_HOST 9
#define HWCART_MD_NODE 10
#define HWCART_MD_CU 11
#define HWCART_MD_CLUSTER 12

#define HWCartOrderXYZ  1
#define HWCartOrderXZY  2
#define HWCartOrderZYX  3
#define HWCartOrderYZX  4
#define HWCartOrderZXY  5
#define HWCartOrderYXZ  6

/* 
   for now:
    - only 3D cartesian space (ndim=3)
    - domain and topo are arrays of size nsplits*ndim
*/

int hwcart_create(MPI_Comm mpi_comm, int nsplits, int *domain, int *topo, int cart_order, MPI_Comm *hwcart_comm_out);
int hwcart_free(MPI_Comm *hwcart_comm);
int hwcart_rank2coord(MPI_Comm hwcart_comm, int *dims, int rank, int cart_order, int *coord_out);
int hwcart_coord2rank(MPI_Comm hwcart_comm, int *dims, int *periodic, int *coord, int cart_order, int *rank_out);

int hwcart_print_rank_topology(MPI_Comm comm, int nlevels, int *domain, int *topo, int order);

#endif /* _HWCART_H */
