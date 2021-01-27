#include "hwcart.h"
#include "hwcart_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>


#ifndef USE_HWLOC
int hwcart_init(hwcart_topo_t *topo_out)
{
    topo_out->ptopo = NULL;
    return 0;
}


int  hwcart_free_hwtopo(hwcart_topo_t *hwtopo)
{
    return 0;
}


int hwcart_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank_out, int level)
{
    int *sbuff, *rbuff;
    int comm_rank, comm_size, nodeid, noderank, color;
    MPI_Comm split_comm;
    int split_type, split_rank, split_size;
    MPI_Comm master_comm;
    int master_size;
    int retval = 0;

    // parent communicator
    HWCART_MPI_CALL( MPI_Comm_rank(comm, &comm_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(comm, &comm_size) );

    if (level == 0) {

        // we've reached the bottom of the topology
        // verify topology validity on this level
        if (comm_size != topo[0]*topo[1]*topo[2]){
            fprintf(stderr, "ERROR: wrong topology on level 0: expected %d domains (config), but found %d (hardware)\n",
                    topo[0]*topo[1]*topo[2], comm_size);
            return -1;
        }
        level_rank_out[level] = comm_rank;
        return 0;
    }

    split_type = hwcart_split_type(domain[level]);
    if(split_type < 0){
        fprintf(stderr, "unknown memory domain %d, level %d\n", domain[level], level);
        return -1;
    }

    // create communicator for this topology level
    HWCART_MPI_CALL( MPI_Comm_split_type(comm, split_type, 0, MPI_INFO_NULL, &split_comm) );
    HWCART_MPI_CALL( MPI_Comm_rank(split_comm, &split_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(split_comm, &split_size) );

    // no split on this topology level
    if (split_size == comm_size) {
        if (1 != topo[level*3+0]*topo[level*3+1]*topo[level*3+2]){
            fprintf(stderr, "ERROR (2): wrong topology on level %d: expected %d domains (config), but found %d (hardware)\n",
                    level, topo[level*3+0]*topo[level*3+1]*topo[level*3+2], 1);
            return -1;
        }
        level_rank_out[level] = 0;
        return hwcart_topology(hwtopo, split_comm, nsplits, domain, topo, level_rank_out, level-1);
    }

    // make a master-rank communicator: masters from each split comm join
    color = 0;
    if (split_rank != 0){

        // non-masters
        HWCART_MPI_CALL( MPI_Comm_split(comm, color, 0, &master_comm) );
    } else {

        // masters
        // temporary nodeid identifier: rank of the split master, +1 needed
        nodeid = comm_rank+1;
        color = 1;
        HWCART_MPI_CALL( MPI_Comm_split(comm, color, 0, &master_comm) );
        HWCART_MPI_CALL( MPI_Comm_size(master_comm, &master_size) );

        // verify topology validity on this level
        if (master_size != topo[level*3+0]*topo[level*3+1]*topo[level*3+2]){
            fprintf(stderr, "ERROR (3): wrong topology on level %d: expected %d domains (config), but found %d (hardware)\n",
                    level, topo[level*3+0]*topo[level*3+1]*topo[level*3+2], master_size);
            return -1;
        }

        // comm buffers to establish unique node id's for each master
        sbuff = (int*)calloc(sizeof(int),comm_size);
        rbuff = (int*)calloc(sizeof(int),comm_size);

        // find node rank based on unique node id from above
        for(int ii=0; ii<comm_size; ii++) sbuff[ii] = nodeid;
        HWCART_MPI_CALL( MPI_Alltoall(sbuff, 1, MPI_INT, rbuff, 1, MPI_INT, master_comm) );

        // mark each unique node id with a 1
        for(int ii=0; ii<comm_size; ii++) sbuff[ii] = 0;
        for(int ii=0; ii<comm_size; ii++) if(rbuff[ii]) sbuff[rbuff[ii]-1] = 1;

        // cumsum: finds node rank for each unique node id
        for(int ii=1; ii<comm_size; ii++) sbuff[ii] = sbuff[ii] + sbuff[ii-1];
        noderank = sbuff[nodeid] - 1;

        // cleanup
        free(sbuff);
        free(rbuff);
    }

    HWCART_MPI_CALL( MPI_Comm_disconnect(&master_comm) );

    // distribute the node id to all split ranks
    HWCART_MPI_CALL( MPI_Bcast(&noderank, 1, MPI_INT, 0, split_comm) );

    // save our level rank
    level_rank_out[level] = noderank;

    // sub-divide lower levels
    retval = hwcart_topology(hwtopo, split_comm, nsplits, domain, topo, level_rank_out, level-1);

    // cleanup
    HWCART_MPI_CALL( MPI_Comm_disconnect(&split_comm) );

    return retval;
}

// obtain level node ID of the calling rank
int hwcart_get_noderank(hwcart_topo_t hwtopo, MPI_Comm comm, int split_type, int *noderank_out)
{
    int *sbuff, *rbuff;
    int nodeid, rank, size, ii;
    MPI_Comm shmem_comm;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int resultlen;

    // old rank
    HWCART_MPI_CALL( MPI_Comm_rank(comm, &rank) );
    HWCART_MPI_CALL( MPI_Comm_size(comm, &size) );

    // communication buffers
    sbuff = calloc(sizeof(int), size);
    rbuff = calloc(sizeof(int), size);

    //check for OpenMPI to use fine-grain hwloc domains
    HWCART_MPI_CALL( MPI_Get_library_version(version, &resultlen) );
    if(0 == strncmp(version, "Open MPI", strlen("Open MPI"))){

        // create local communicator
        split_type = hwcart_split_type(split_type);
        if(split_type < 0){
            fprintf(stderr, "unknown memory domain %d\n", split_type);
            return -1;
        }
        HWCART_MPI_CALL( MPI_Comm_split_type(comm, split_type, 0, MPI_INFO_NULL, &shmem_comm) );
    } else {

        // create local communicator
        HWCART_MPI_CALL( MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm) );
    }

    // figure out unique compute node id = max rank id on each compute node
    sbuff[0] = rank+1;
    HWCART_MPI_CALL( MPI_Allreduce(sbuff, rbuff, 1, MPI_INT, MPI_MAX, shmem_comm) );
    nodeid = rbuff[0];

    // find node rank based on unique node id from above
    for(ii=0; ii<size; ii++) sbuff[ii] = nodeid;
    HWCART_MPI_CALL( MPI_Alltoall(sbuff, 1, MPI_INT, rbuff, 1, MPI_INT, comm) );

    // mark each unique node id with a 1
    for(ii=0; ii<size; ii++) sbuff[ii] = 0;
    for(ii=0; ii<size; ii++) sbuff[rbuff[ii]] = 1;

    // cumsum: finds node rank for each unique node id
    for(ii=1; ii<size; ii++){
        sbuff[ii] = sbuff[ii-1] + sbuff[ii];
    }
    *noderank_out = sbuff[nodeid-1];

    // cleanup
    free(sbuff);
    free(rbuff);
    HWCART_MPI_CALL( MPI_Comm_disconnect(&shmem_comm) );

    return 0;
}

int hwcart_split_type(int split_type)
{
    switch (split_type) {
    case (HWCART_MD_HWTHREAD):
        return OMPI_COMM_TYPE_HWTHREAD;
    case (HWCART_MD_CORE):
        return OMPI_COMM_TYPE_CORE;
    case (HWCART_MD_L1CACHE):
        return OMPI_COMM_TYPE_L1CACHE;
    case (HWCART_MD_L2CACHE):
        return OMPI_COMM_TYPE_L2CACHE;
    case (HWCART_MD_L3CACHE):
        return OMPI_COMM_TYPE_L3CACHE;
    case (HWCART_MD_SOCKET):
        return OMPI_COMM_TYPE_SOCKET;
    case (HWCART_MD_NUMA):
        return OMPI_COMM_TYPE_NUMA;
    case (HWCART_MD_NODE):
        return OMPI_COMM_TYPE_NODE;
    default:
        return -1;
    }
}

#endif
