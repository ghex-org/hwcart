#include <hwcart/hwcart.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{

#define NLEVELS 5
    hwcart_split_t domain[NLEVELS] = {
                               HWCART_MD_CORE,
                               HWCART_MD_L3CACHE,
                               HWCART_MD_NUMA,
                               HWCART_MD_SOCKET,
                               HWCART_MD_NODE
    };
    int topo[3*NLEVELS] = {
                           2, 2, 1, // core grid
                           1, 1, 1, // l3cache grid
                           1, 1, 1, // numa grid
                           1, 1, 1, // socket grid
                           1, 1, 1, // node grid
    };
    
    int periodic[3] = {1, 1, 1};
    
    int comm_rank, comm_size;
    hwcart_order_t cart_order = HWCartOrderXYZ;
    hwcart_topo_t hwtopo;
    MPI_Comm hwcart_comm;
    MPI_Comm mpicart_comm;

    int ii, gdim[3] = {1,1,1};
    for(ii=0; ii<NLEVELS; ii++){
      gdim[0] *= topo[ii*3+0];
      gdim[1] *= topo[ii*3+1];
      gdim[2] *= topo[ii*3+2];
    }

    if(argc>1){
	cart_order = (hwcart_order_t)atoi(argv[1]);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if(hwcart_init(&hwtopo)){
        MPI_Finalize();
        exit(1);
    }
    
    if(!hwcart_create(hwtopo, MPI_COMM_WORLD, NLEVELS, domain, topo, periodic, cart_order, &hwcart_comm)){
        hwcart_print_rank_topology(hwtopo, hwcart_comm);

	// Convert a HWCART communicator to a native MPI_Cart communicator.
	// Dimension order will be set to MPI order (ZYX). All standard MPI_Cart functions
	// (MPI_Cart_shift, MPI_Cart_sub, MPI_Cart_rank, etc.) can be used on mpicart_com.
	hwcart2mpicart(hwcart_comm, &mpicart_comm);

	int  hwcart_coord[3], hwcart_rank;
	int mpicart_coord[3], mpicart_rank;
	
	// hwcart cartesian coordinates
	MPI_Comm_rank(hwcart_comm, &hwcart_rank);
	hwcart_rank2coord(hwcart_comm, hwcart_rank, hwcart_coord);
	
	// mpicart cartesian coordinates
	MPI_Comm_rank(mpicart_comm, &mpicart_rank);
	MPI_Cart_coords(mpicart_comm, mpicart_rank, 3, mpicart_coord);

	if(0 != memcmp(mpicart_coord, hwcart_coord, sizeof(int)*3)){
	  fprintf(stderr, "ERROR: inconsistent MPI_Cart and hwcart coordinates\n");
	  exit(1);
	}
	
	MPI_Comm_free(&mpicart_comm);
	hwcart_comm_free(&hwcart_comm);
        hwcart_topo_free(&hwtopo);
    } else exit(1);

    MPI_Finalize();
}
