#include "hwcart.h"
#include "hwcart_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>

#ifdef USE_HWLOC
#include <hwloc.h>
int hwcart_topology(MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank_out, int level)
{
    int *sbuff, *rbuff;
    int ierr, comm_rank, comm_size, nodeid, noderank, color, ii, split_type;
    MPI_Comm shmem_comm, split_comm;
    int shmem_rank, shmem_size, split_rank, split_size;
    MPI_Comm master_comm;
    int master_size;
    int retval = 0;
    int nlevel_nodes;
    for(int i=0; i<nsplits; i++) level_rank_out[i] = -1;

    hwloc_topology_t hwloctopo;
    int res = hwloc_topology_init (&hwloctopo);
    res = hwloc_topology_load(hwloctopo);

    // parent communicator
    ierr = MPI_Comm_rank(comm, &comm_rank);
    ierr = MPI_Comm_size(comm, &comm_size);

    if(domain[level] != HWCART_MD_NODE){
        fprintf(stderr, "top-1 memory domain must be HWCART_MD_NODE\n");
        hwloc_topology_destroy (hwloctopo);
        return -1;
    }

    // create communicator for this topology level
    ierr = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm);
    ierr = MPI_Comm_rank(shmem_comm, &shmem_rank);
    ierr = MPI_Comm_size(shmem_comm, &shmem_size);

    // make a master-rank communicator: masters from each split comm join
    color = 0;
    if (shmem_rank != 0){

        // non-masters
        ierr = MPI_Comm_split(comm, color, 0, &master_comm);
    } else {

        // masters
        color = 1;
        ierr = MPI_Comm_split(comm, color, 0, &master_comm);
        ierr = MPI_Comm_rank(master_comm, level_rank_out+level);
    }
    
    // cleanup
    ierr = MPI_Comm_disconnect(&master_comm);
    
    // distribute the node id to all split ranks
    ierr = MPI_Bcast(level_rank_out+level, 1, MPI_INT, 0, shmem_comm);

    hwloc_bitmap_t m_cpuset = hwloc_bitmap_alloc();
    hwloc_obj_t obj;
    color = 0;
    for(int i=level-1; i>=0; i--){
        split_type = hwcart_split_type(domain[i]);
        if(split_type < 0){
            fprintf(stderr, "unknown memory domain %d, level %d\n", domain[i], i+1);
            hwloc_bitmap_free(m_cpuset);
            hwloc_topology_destroy (hwloctopo);
            return -1;
        }

        int n = hwloc_get_nbobjs_by_type(hwloctopo, split_type);
        for(int j=0; j<n; j++){
            // figure out on which object we reside
            obj = hwloc_get_obj_by_type (hwloctopo, split_type, j);
            hwloc_get_cpubind(hwloctopo, m_cpuset, 0);
            if(hwloc_bitmap_isincluded (m_cpuset, obj->cpuset)){
                // include myself on this node
                level_rank_out[i] = j;
                break;
            }
        }

        // find number of nodes on this level
        ierr = MPI_Allreduce(level_rank_out+i, &nlevel_nodes, 1, MPI_INT, MPI_MAX, shmem_comm);
        nlevel_nodes++;

        if (nlevel_nodes != topo[i*3+0]*topo[i*3+1]*topo[i*3+2]){
            
            // those either have to match, or there is no split on this level
            
            // no split
            if((i != 0) && (1 == topo[i*3+0]*topo[i*3+1]*topo[i*3+2])) {
                level_rank_out[i] = 0;
                continue;
            }

            fprintf(stderr, "ERROR: wrong topology on level %d: expected %d nodes on this level, instead found %d\n",
                    i, topo[i*3+0]*topo[i*3+1]*topo[i*3+2], nlevel_nodes);
            return -1;
        } else {
            // no split
            if(1 == topo[i*3+0]*topo[i*3+1]*topo[i*3+2]) {
                level_rank_out[i] = 0;
                continue;
            }
        }

        // create a sub-topology
        hwloc_topology_restrict(hwloctopo, obj->cpuset, 0);
    }

    // cleanup
    hwloc_bitmap_free(m_cpuset);
    hwloc_topology_destroy (hwloctopo);

    return 0;
}

// obtain level node ID of the calling rank
int hwcart_get_noderank(MPI_Comm comm, int split_type, int *noderank_out)
{
    int *sbuff, *rbuff;
    int nodeid, rank, size, ierr, ii;
    MPI_Comm shmem_comm;
    int resultlen;

    if(split_type == HWCART_MD_NODE){
    
        // old rank
        ierr = MPI_Comm_rank(comm, &rank);
        ierr = MPI_Comm_size(comm, &size);

        sbuff = calloc(sizeof(int), size);
        rbuff = calloc(sizeof(int), size);

        // create local communicator
        ierr = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm);

        // figure out unique compute node id = max rank id on each compute node
        sbuff[0] = rank+1;
        ierr = MPI_Allreduce(sbuff, rbuff, 1, MPI_INT, MPI_MAX, shmem_comm);
        nodeid = rbuff[0];

        // find node rank based on unique node id from above
        for(ii=0; ii<size; ii++) sbuff[ii] = nodeid;
        ierr = MPI_Alltoall(sbuff, 1, MPI_INT, rbuff, 1, MPI_INT, comm);

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
        ierr = MPI_Comm_disconnect(&shmem_comm);
        return 0;
    }
  
    hwloc_topology_t hwloctopo;
    int res = hwloc_topology_init (&hwloctopo);
    res = hwloc_topology_load(hwloctopo);
  
    hwloc_bitmap_t m_cpuset = hwloc_bitmap_alloc();
    hwloc_obj_t obj;
    int hwloc_split_type = hwcart_split_type(split_type);
    if(hwloc_split_type < 0){
        fprintf(stderr, "unknown memory domain %d\n", split_type);
        hwloc_bitmap_free(m_cpuset);
        hwloc_topology_destroy (hwloctopo);
        return -1;
    }

    int n = hwloc_get_nbobjs_by_type(hwloctopo, hwloc_split_type);
    for(int j=0; j<n; j++){
        // figure out on which object we reside
        obj = hwloc_get_obj_by_type (hwloctopo, hwloc_split_type, j);
        hwloc_get_cpubind(hwloctopo, m_cpuset, 0);
        if(hwloc_bitmap_isincluded (m_cpuset, obj->cpuset)){
            // include myself on this node
            *noderank_out = j;
            break;
        }
    }
    hwloc_bitmap_free(m_cpuset);
    hwloc_topology_destroy (hwloctopo);
}

int hwcart_split_type(int split_type)
{
    switch (split_type) {
    case (HWCART_MD_HWTHREAD):
        return HWLOC_OBJ_PU;
    case (HWCART_MD_CORE):
        return HWLOC_OBJ_CORE;
    case (HWCART_MD_L1CACHE):
        return HWLOC_OBJ_L1CACHE;
    case (HWCART_MD_L2CACHE):
        return HWLOC_OBJ_L2CACHE;
    case (HWCART_MD_L3CACHE):
        return HWLOC_OBJ_L3CACHE;
    case (HWCART_MD_SOCKET):
        return HWLOC_OBJ_SOCKET;
    case (HWCART_MD_NUMA):
        return HWLOC_OBJ_NUMANODE;
    default:
        return -1;
    }
}

#endif
