#include <hwcart/hwcart.h>
#include <hwcart/hwcart_utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>


int hwcart_remap_ranks(MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, int *level_rank, hwcart_order_t order, MPI_Comm *hwcart_comm_out)
{
    int topo_coord[3] = {0}, gdim[3] = {1,1,1}, periodic[3] = {0};
    int ii, comm_rank, newrank;
    int *cartXYZ;
    int *dims;

    HWCART_MPI_CALL( MPI_Comm_rank(comm, &comm_rank) );

    cartXYZ = calloc(3*nlevels, sizeof(int));
    dims    = calloc(3*nlevels, sizeof(int));

    // compute rank cartesian coordinates for each topological level
    for(ii=0; ii<nlevels; ii++){
        hwcart_rank2coord(MPI_COMM_NULL, topo+3*ii, level_rank[ii], order, cartXYZ+3*ii);
    }

    // compute rank grid resulution for each level
    dims[0] = 1;
    dims[1] = 1;
    dims[2] = 1;
    for(ii=1; ii<nlevels; ii++){
        dims[ii*3+0] = dims[(ii-1)*3+0]*topo[(ii-1)*3+0];
        dims[ii*3+1] = dims[(ii-1)*3+1]*topo[(ii-1)*3+1];
        dims[ii*3+2] = dims[(ii-1)*3+2]*topo[(ii-1)*3+2];
    }

    // compute resulting global cartesian index
    for(ii=0; ii<nlevels; ii++){
        topo_coord[0] = topo_coord[0] + cartXYZ[ii*3+0]*dims[ii*3+0];
        topo_coord[1] = topo_coord[1] + cartXYZ[ii*3+1]*dims[ii*3+1];
        topo_coord[2] = topo_coord[2] + cartXYZ[ii*3+2]*dims[ii*3+2];
    }

    // compute global grid dimensions
    for(ii=0; ii<nlevels; ii++){
        gdim[0] *= topo[ii*3+0];
        gdim[1] *= topo[ii*3+1];
        gdim[2] *= topo[ii*3+2];
    }

    // compute rank id in global rank space
    hwcart_coord2rank(comm, gdim, periodic, topo_coord, order, &newrank);

    // create the new communicator with remapped ranks
    HWCART_MPI_CALL( MPI_Comm_split(comm, 0, newrank, hwcart_comm_out) );

    // cleanup
    free(dims);
    free(cartXYZ);

    return 0;
}


int  hwcart_create(hwcart_topo_t hwtopo, MPI_Comm mpi_comm, int nlevels, hwcart_split_t *domain, int *topo, int *periodic, hwcart_order_t order, MPI_Comm *hwcart_comm_out)
{
    int retval;
    int *level_rank;
    int ii, gdim[3] = {1,1,1};
    int comm_size;

    if(domain[nlevels-1] != HWCART_MD_NODE){
        fprintf(stderr, "top memory domain must be HWCART_MD_NODE\n");
        return -1;
    }

    // verify number of ranks vs. the config
    for(ii=0; ii<nlevels; ii++){
        gdim[0] *= topo[ii*3+0];
        gdim[1] *= topo[ii*3+1];
        gdim[2] *= topo[ii*3+2];
    }
    HWCART_MPI_CALL( MPI_Comm_size(mpi_comm, &comm_size) );
    if(comm_size != gdim[0]*gdim[1]*gdim[2]) {
        fprintf(stderr, "number of ranks (%d) is different than the specified rank grid size (%d)\n",
		comm_size, gdim[0]*gdim[1]*gdim[2]);
        return -1;
    }
    
    level_rank = calloc(nlevels, sizeof(int));
    retval = hwcart_topology(hwtopo, mpi_comm, nlevels, domain, topo, level_rank, nlevels-1);
    if(retval<0) {
        free(level_rank);
        return retval;
    }

    if(0){
        MPI_Comm temp_comm;
        retval = hwcart_remap_ranks(mpi_comm, nlevels, domain, topo, level_rank, order, &temp_comm);
        HWCART_MPI_CALL( MPI_Cart_create(temp_comm, 3, gdim, periodic, 0, hwcart_comm_out) );
        HWCART_MPI_CALL( MPI_Comm_free(&temp_comm) );
    } else {
        retval = hwcart_remap_ranks(mpi_comm, nlevels, domain, topo, level_rank, order, hwcart_comm_out);
    }

#ifdef DEBUG
    int comm_rank, new_rank;
    MPI_Comm_rank(mpi_comm, &comm_rank);
    MPI_Comm_rank(*hwcart_comm_out, &new_rank);
    printf("%d -> %d: level_rank ", comm_rank, new_rank);
    for(int i=0; i<nlevels; i++) printf("%d ", level_rank[i]); printf("\n");
#endif

    free(level_rank);

    return retval;
}


int hwcart_comm_free(MPI_Comm *hwcart_comm)
{
    if(hwcart_comm){
      HWCART_MPI_CALL( MPI_Comm_free(hwcart_comm) );
    }
    return 0;
}


int  hwcart_sub(MPI_Comm comm, int *dims, int rank, hwcart_order_t order, int *belongs, MPI_Comm *hwcart_comm_out)
{
  int topo = -1, color;
  int coord[3], periodic[3] = {0};

  // check if this is a cartesian communicator
  if (comm != MPI_COMM_NULL) {
    HWCART_MPI_CALL( MPI_Topo_test(comm, &topo) );
  }

  if (topo != MPI_CART) {

    // Find ranks belonging to the new communicator based on each rank's cartesian coordinates.
    // color is computed by zeroing out 'belongs' in the cartesian coordinate of each rank
    // and then computing the 'collapsed' rank by coord2rank
    hwcart_rank2coord(comm, dims, rank, order, coord);
    coord[0] = coord[0]*(belongs[0]==0);
    coord[1] = coord[1]*(belongs[1]==0);
    coord[2] = coord[2]*(belongs[2]==0);
    hwcart_coord2rank(comm, dims, periodic, coord, order, &color);

    // create the new communicator
    HWCART_MPI_CALL( MPI_Comm_split(comm, color, 0, hwcart_comm_out) );
  } else {
    HWCART_MPI_CALL( MPI_Cart_sub(comm, belongs, hwcart_comm_out) );
  }
  return 0;
}


int  hwcart_rank2coord(MPI_Comm comm, int *dims, int rank, hwcart_order_t order, int *coord)
{
    int topo, tmp;

    // check if this is a cartesian communicator
    if (comm == MPI_COMM_NULL) {
        topo = -1;
    } else {
        HWCART_MPI_CALL( MPI_Topo_test(comm, &topo) );
    }

    if (topo!=MPI_CART){
        switch (order){
        case (HWCartOrderXYZ):
            tmp = rank;        coord[0] = tmp%dims[0];
            tmp = tmp/dims[0]; coord[1] = tmp%dims[1];
            tmp = tmp/dims[1]; coord[2] = tmp;
            break;
        case (HWCartOrderXZY):
            tmp = rank;        coord[0] = tmp%dims[0];
            tmp = tmp/dims[0]; coord[2] = tmp%dims[2];
            tmp = tmp/dims[2]; coord[1] = tmp;
            break;
        case (HWCartOrderZYX):
            tmp = rank;        coord[2] = tmp%dims[2];
            tmp = tmp/dims[2]; coord[1] = tmp%dims[1];
            tmp = tmp/dims[1]; coord[0] = tmp;
            break;
        case (HWCartOrderYZX):
            tmp = rank;        coord[1] = tmp%dims[1];
            tmp = tmp/dims[1]; coord[2] = tmp%dims[2];
            tmp = tmp/dims[2]; coord[0] = tmp;
            break;
        case (HWCartOrderZXY):
            tmp = rank;        coord[2] = tmp%dims[2];
            tmp = tmp/dims[2]; coord[0] = tmp%dims[0];
            tmp = tmp/dims[0]; coord[1] = tmp;
            break;
        case (HWCartOrderYXZ):
            tmp = rank;        coord[1] = tmp%dims[1];
            tmp = tmp/dims[1]; coord[0] = tmp%dims[0];
            tmp = tmp/dims[0]; coord[2] = tmp;
            break;
        default:
            fprintf(stderr, "unknown value of argument 'order': %d\n", order);
            return -1;
        }
    } else {
        HWCART_MPI_CALL( MPI_Cart_coords(comm, rank, 3, coord) );
    }
    return 0;
}


int  hwcart_coord2rank(MPI_Comm comm, int *dims, int *periodic, int *coord, hwcart_order_t order, int *rank_out)
{
    int tcoord[3], ii, topo;

    // apply periodicity inside a temporary array
    tcoord[0] = coord[0];
    tcoord[1] = coord[1];
    tcoord[2] = coord[2];

    // wrap-around negative cartesian indices.
    // TODO: currently only correctly handles -1 and dim
    for(ii=0; ii<3; ii++){
        if(tcoord[ii] < 0){
            if(!periodic[ii] || tcoord[ii]!=-1){
                *rank_out = MPI_PROC_NULL;
                return 0;
            }
            tcoord[ii] = dims[ii]-1;
        }
        if(tcoord[ii] >= dims[ii]){
            if(!periodic[ii] || tcoord[ii]>dims[ii]){
                *rank_out = MPI_PROC_NULL;
                return 0;
            }
            tcoord[ii] = 0;
        }
    }

    // check if this is a cartesian communicator
    if (comm == MPI_COMM_NULL){
        topo = -1;
    } else {
        HWCART_MPI_CALL( MPI_Topo_test(comm, &topo) );
    }

    if (topo!=MPI_CART) {

        switch (order) {
        case (HWCartOrderXYZ):
            *rank_out = (tcoord[2]*dims[1] + tcoord[1])*dims[0] + tcoord[0];
            break;
        case (HWCartOrderXZY):
            *rank_out = (tcoord[1]*dims[2] + tcoord[2])*dims[0] + tcoord[0];
            break;
        case (HWCartOrderZYX):
            *rank_out = (tcoord[0]*dims[1] + tcoord[1])*dims[2] + tcoord[2];
            break;
        case (HWCartOrderYZX):
            *rank_out = (tcoord[0]*dims[2] + tcoord[2])*dims[1] + tcoord[1];
            break;
        case (HWCartOrderZXY):
            *rank_out = (tcoord[1]*dims[0] + tcoord[0])*dims[2] + tcoord[2];
            break;
        case (HWCartOrderYXZ):
            *rank_out = (tcoord[2]*dims[0] + tcoord[0])*dims[1] + tcoord[1];
            break;
        default:
            fprintf(stderr, "unknown value of argument 'order': %d\n", order);
            return -1;
        }
    } else {
        HWCART_MPI_CALL( MPI_Cart_rank(comm, tcoord, rank_out) );
    }
    return 0;
}


int hwcart_print_rank_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t order)
{
    int *buff;
    int ii, gdim[3] = {1,1,1};
    int comm_rank, comm_size, orank, ncpus;
    int *sbuff;
    int *level_id;
    char split_name[256];

    // compute global grid dimensions
    for(ii=0; ii<nlevels; ii++){
        gdim[0] *= topo[ii*3+0];
        gdim[1] *= topo[ii*3+1];
        gdim[2] *= topo[ii*3+2];
    }

    level_id = calloc(nlevels, sizeof(int));
    for(ii=0; ii<nlevels; ii++){
        hwcart_get_noderank(hwtopo, comm, domain[ii], level_id+ii);
    }

    sbuff = calloc(nlevels+3, sizeof(int));

    // obtain all values at master
    HWCART_MPI_CALL( MPI_Comm_rank(MPI_COMM_WORLD, &orank) );
    HWCART_MPI_CALL( MPI_Comm_rank(comm, &comm_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(comm, &comm_size) );

    buff = calloc((nlevels+3)*comm_size, sizeof(int));
    sbuff[0] = comm_rank;
    sbuff[1] = orank;
    sbuff[2] = sched_getcpu(); // -D_GNU_SOURCE
    for(ii=0; ii<nlevels; ii++) sbuff[3+ii] = level_id[ii];

    // assuming HT is enabled, use the lower core ID
    ncpus = get_nprocs_conf()/2;
    if (sbuff[2] >= ncpus) sbuff[2] = sbuff[2] - ncpus;
    HWCART_MPI_CALL( MPI_Gather(sbuff, nlevels+3, MPI_INT, buff, nlevels+3, MPI_INT, 0, comm) );

    if (comm_rank==0) {
        printf("\n");
        switch(order){
        case(HWCartOrderXYZ):
            printf("cart order XYZ\n"); break;
        case(HWCartOrderXZY):
            printf("cart order XZY\n"); break;
        case(HWCartOrderZYX):
            printf("cart order ZYX\n"); break;
        case(HWCartOrderYZX):
            printf("cart order YZX\n"); break;
        case(HWCartOrderZXY):
            printf("cart order ZXY\n"); break;
        case(HWCartOrderYXZ):
            printf("cart order YXZ\n"); break;
        default:
            fprintf(stderr, "Unknown cart order %d\n", order);
            return -1;
        }

        printf("\n");
        printf("Rank to node mapping\n");
        printf("\n");

        for(ii=0; ii<nlevels; ii++){
            hwcart_split_type_to_name(domain[ii], split_name);
            printf("\n");
            printf("Level %d %s\n", ii, split_name);
            printf("\n");

            hwcart_print_cube(comm, gdim, 3+ii, buff, nlevels+3, order);
        }

        printf("\n");
        printf("Rank layout in HWCART communicator\n");
        printf("\n");
        hwcart_print_cube(comm, gdim, 0, buff, nlevels+3, order);

        printf("\n");
        printf("MPI_COMM_WORLD layout\n");
        printf("\n");
        hwcart_print_cube(comm, gdim, 1, buff, nlevels+3, order);

        printf("\n");
        printf("CPU binding\n");
        printf("\n");
        hwcart_print_cube(comm, gdim, 2, buff, nlevels+3, order);

        printf("\n");
        printf("Z |    / Y   \n");
        printf("  |   /      \n");
        printf("  |  /       \n");
        printf("  | /        \n");
        printf("  |_______ X \n");
        printf("\n");

    }

    free(level_id);
    free(sbuff);
    
    return 0;
}


int hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, hwcart_order_t order){
    int k, j, i, kk, n;
    int comm_size;
    int periodic[3] = {0}, idx[3];

    HWCART_MPI_CALL( MPI_Comm_size(comm, &comm_size) );

    for(k=gdim[2]-1; k>=0; k--){
        for(j=gdim[1]-1; j>=0; j--){
            for(kk=0; kk<j*2+5; kk++) printf(" ");
            for(i=0; i<gdim[0]; i++){
                idx[0] = i;
                idx[1] = j;
                idx[2] = k;
                hwcart_coord2rank(comm, gdim, periodic, idx, order, &n);
                if(comm_size < 1000)
                    printf("%4d", buff[n*line_size + id]);
                else
                    printf("%5d", buff[n*line_size + id]);
            }
            printf("\n");
        }
        printf("\n");
    }
    return 0;
}


void hwcart_split_type_to_name(int split_type, char *name) {

    switch (split_type) {
    case (HWCART_MD_HWTHREAD):
        sprintf(name, "HWCART_MD_HWTHREAD");
        break;
    case (HWCART_MD_CORE):
        sprintf(name, "HWCART_MD_CORE");
        break;
    case (HWCART_MD_L1CACHE):
        sprintf(name, "HWCART_MD_L1CACHE");
        break;
    case (HWCART_MD_L2CACHE):
        sprintf(name, "HWCART_MD_L2CACHE");
        break;
    case (HWCART_MD_L3CACHE):
        sprintf(name, "HWCART_MD_L3CACHE");
        break;
    case (HWCART_MD_SOCKET):
        sprintf(name, "HWCART_MD_SOCKET");
        break;
    case (HWCART_MD_NUMA):
        sprintf(name, "HWCART_MD_NUMA");
        break;
    case (HWCART_MD_NODE):
        sprintf(name, "HWCART_MD_NODE");
        break;
    }
}


#ifdef HWCART_USE_FORTRAN

/* Fortran versions */
int hwcart_create_f(hwcart_topo_t hwtopo, int mpi_comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t order, int *hwcart_comm_out)
{
    /* MPI_Comm in_comm = MPI_Comm_f2c(mpi_comm), out_comm; */
    /* int retval = hwcart_create(hwtopo, in_comm, nlevels, domain, topo, order, &out_comm); */
    /* *hwcart_comm_out = MPI_Comm_c2f(out_comm); */
    /* return retval; */
}


int  hwcart_sub_f(int mpi_comm, int *dims, int rank, hwcart_order_t order, int *belongs, int *hwcart_comm_out)
{
    MPI_Comm in_comm = MPI_Comm_f2c(mpi_comm), out_comm;
    int retval = hwcart_sub(in_comm, dims, rank, order, belongs, &out_comm);
    *hwcart_comm_out = MPI_Comm_c2f(out_comm);
    return retval;
}

int hwcart_comm_free_f(int *hwcart_comm)
{
    if(hwcart_comm){
        MPI_Comm in_comm = MPI_Comm_f2c(*hwcart_comm);
        *hwcart_comm = MPI_Comm_c2f(MPI_COMM_NULL);
	if(in_comm != MPI_COMM_NULL){
	  HWCART_MPI_CALL( MPI_Comm_free(&in_comm) );
	}
    }
    return 0;
}

int hwcart_print_rank_topology_f(hwcart_topo_t hwtopo, int comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t order)
{
    MPI_Comm in_comm = MPI_Comm_f2c(comm);
    return hwcart_print_rank_topology(hwtopo, in_comm, nlevels, domain, topo, order);
}

int hwcart_rank2coord_f(int comm, int *dims, int rank, hwcart_order_t order, int *coord_out)
{
    MPI_Comm in_comm = MPI_Comm_f2c(comm);
    return hwcart_rank2coord(in_comm, dims, rank, order, coord_out);
}

int hwcart_coord2rank_f(int comm, int *dims, int *periodic, int *coord, hwcart_order_t order, int *rank_out)
{
    MPI_Comm in_comm = MPI_Comm_f2c(comm);
    return hwcart_coord2rank(in_comm, dims, periodic, coord, order, rank_out);
}

#endif
