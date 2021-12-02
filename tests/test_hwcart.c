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
    hwcart_order_t cart_order = HWCartOrderZYX;
    MPI_Comm hwcart_comm;

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

    if(hwcart_init()){
        MPI_Finalize();
        exit(1);
    }
    
    if(!hwcart_create(MPI_COMM_WORLD, NLEVELS, domain, topo, periodic, cart_order, &hwcart_comm)){
        hwcart_print_rank_topology(hwcart_comm);
	MPI_Barrier(hwcart_comm);

	// get our cartesian coordinates
	int hwcart_rank, coord[3], nbleft, nbright;
	MPI_Comm_rank(hwcart_comm, &hwcart_rank);
	hwcart_rank2coord(hwcart_comm, hwcart_rank, coord);

	// query left and right neighbors (in X dimension)
	coord[0]-=1;
	hwcart_coord2rank(hwcart_comm, coord, &nbleft);
	coord[0]+=2;
	hwcart_coord2rank(hwcart_comm, coord, &nbright);
	coord[0]-=1;

	{
	    int belongs[3];
	    int col_size, plane_size;
	    MPI_Comm hwcart_comm_col, hwcart_comm_plane;

	    // create column communicator
	    belongs[0] = 0; belongs[1] = 0; belongs[2] = 1;
	    hwcart_sub(hwcart_comm, hwcart_rank, belongs, &hwcart_comm_col);
	    MPI_Comm_size(hwcart_comm_col, &col_size);
	    if (col_size != gdim[2]){
		fprintf(stderr, "ERROR: wrong size of the column communicator: %d /= %d", col_size, gdim[2]);
		exit(1);
	    }

	    // create plane communicator
	    belongs[0] = 1; belongs[1] = 1; belongs[2] = 0;
	    hwcart_sub(hwcart_comm, hwcart_rank, belongs, &hwcart_comm_plane);
	    MPI_Comm_size(hwcart_comm_plane, &plane_size);
	    if (plane_size != gdim[0]*gdim[1]){
		fprintf(stderr, "ERROR: wrong size of the plane communicator: %d /= %d", col_size, gdim[0]*gdim[1]);
		exit(1);
	    }

	    hwcart_comm_free(&hwcart_comm_col);
	    hwcart_comm_free(&hwcart_comm_plane);
	}

	hwcart_comm_free(&hwcart_comm);
    } else exit(1);

    hwcart_finalize();
    MPI_Finalize();
}
