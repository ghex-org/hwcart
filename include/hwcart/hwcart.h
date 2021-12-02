#ifndef _HWCART_H
#define _HWCART_H

#include <mpi.h>

typedef enum {
              HWCART_MD_HWTHREAD,
              HWCART_MD_CORE,
              HWCART_MD_L1CACHE,
              HWCART_MD_L2CACHE,
              HWCART_MD_L3CACHE,
              HWCART_MD_SOCKET,
              HWCART_MD_NUMA,
              HWCART_MD_NODE
} hwcart_split_t;

typedef enum {
              HWCartOrderXYZ,
              HWCartOrderXZY,
              HWCartOrderZYX,
              HWCartOrderYZX,
              HWCartOrderZXY,
              HWCartOrderYXZ
} hwcart_order_t;

/* 
   for now:
   - only 3D cartesian space (ndim=3)
   - domain and topo are arrays of size nlevels*ndim
*/

int hwcart_init();
int hwcart_finalize();

int hwcart_create(MPI_Comm mpi_comm,
                  int nlevels,
                  hwcart_split_t *domain,
                  int *topo,
                  int *periodic,
                  hwcart_order_t cart_order,
                  MPI_Comm *hwcart_comm_out);

int hwcart_comm_free(MPI_Comm *hwcart_comm);

int hwcart2mpicart(MPI_Comm hwcart_comm, MPI_Comm *mpicart_comm_out);

int hwcart_sub(MPI_Comm comm, int rank, int *belongs, MPI_Comm *hwcart_comm_out);

int hwcart_rank2coord(MPI_Comm hwcart_comm, int rank, int *coord_out);
int hwcart_coord2rank(MPI_Comm hwcart_comm, int *coord, int *rank_out);

int hwcart_print_rank_topology(MPI_Comm comm);

#endif /* _HWCART_H */
