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

typedef struct hwcart_topo_struct_t* hwcart_topo_t;

int hwcart_init(hwcart_topo_t *hwtopo_out);
int hwcart_topo_free(hwcart_topo_t *hwtopo);

int hwcart_create(hwcart_topo_t hwtopo, MPI_Comm mpi_comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t cart_order, MPI_Comm *hwcart_comm_out);
int  hwcart2mpicart(MPI_Comm hwcart_comm, int nlevels, int *topo, int *periodic, hwcart_order_t cart_order, MPI_Comm *mpicart_comm_out);

int hwcart_sub(MPI_Comm comm, int *dims, int rank, hwcart_order_t order, int *belongs, MPI_Comm *hwcart_comm_out);
int hwcart_comm_free(MPI_Comm *hwcart_comm);

int hwcart_rank2coord(MPI_Comm hwcart_comm, int *dims, int rank, hwcart_order_t cart_order, int *coord_out);
int hwcart_coord2rank(MPI_Comm hwcart_comm, int *dims, int *periodic, int *coord, hwcart_order_t cart_order, int *rank_out);

int hwcart_print_rank_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t order);

#endif /* _HWCART_H */
