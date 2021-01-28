#include "hwcart.h"
#include "hwcart_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>

#ifdef USE_HWLOC
#include <hwloc.h>

hwloc_obj_type_t hwcart_split_type(hwcart_split_t split_type);

struct hwcart_topo_struct_t {
    hwloc_topology_t topo;
};


int hwcart_init(hwcart_topo_t *hwtopo_out)
{
    int res;
    *hwtopo_out = malloc(sizeof(struct hwcart_topo_struct_t));
    res = hwloc_topology_init(&(*hwtopo_out)->topo);
    if(res) {
        free(*hwtopo_out);
        *hwtopo_out = NULL;
        return res;
    }
    res = hwloc_topology_load((*hwtopo_out)->topo);
    if(res) {
        hwloc_topology_destroy ((*hwtopo_out)->topo);
        free(*hwtopo_out);
        *hwtopo_out = NULL;        
        return res;
    }
    return 0;
}


int  hwcart_free_hwtopo(hwcart_topo_t *hwtopo)
{
    hwloc_topology_destroy((*hwtopo)->topo);
    free(*hwtopo);
    *hwtopo = NULL;
    return 0;
}


int hwcart_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, int *level_rank_out, int level)
{
    int comm_rank, comm_size;
    hwloc_obj_type_t split_type;
    MPI_Comm shmem_comm;
    int shmem_rank, shmem_size;
    MPI_Comm master_comm;
    int nlevel_nodes;

    // init
    for(int i=0; i<nlevels; i++) level_rank_out[i] = -1;

    // parent communicator
    HWCART_MPI_CALL( MPI_Comm_rank(comm, &comm_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(comm, &comm_size) );

    // create communicator for shared memory nodes
    HWCART_MPI_CALL( MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm) );
    HWCART_MPI_CALL( MPI_Comm_rank(shmem_comm, &shmem_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(shmem_comm, &shmem_size) );

    // make a master-rank communicator: masters from each split comm join
    if (shmem_rank != 0){

        // non-masters
        HWCART_MPI_CALL( MPI_Comm_split(comm, 0, 0, &master_comm) );
    } else {

        // masters
        HWCART_MPI_CALL( MPI_Comm_split(comm, 1, 0, &master_comm) );
        HWCART_MPI_CALL( MPI_Comm_rank(master_comm, level_rank_out+level) );
    }
    
    // cleanup
    HWCART_MPI_CALL( MPI_Comm_disconnect(&master_comm) );
    
    // distribute the node id to all split ranks
    HWCART_MPI_CALL( MPI_Bcast(level_rank_out+level, 1, MPI_INT, 0, shmem_comm) );

    // we need a copy of the topology: we modify it
    hwloc_topology_t hwloctopo;
    hwloc_topology_dup(&hwloctopo, hwtopo->topo);
    
    hwloc_bitmap_t m_obj_cpuset = hwloc_bitmap_alloc();
    hwloc_bitmap_t m_cpuset = hwloc_bitmap_alloc();
    hwloc_obj_t obj;
    hwloc_get_cpubind(hwloctopo, m_cpuset, 0);
    
    for(int i=level-1; i>=0; i--){
        split_type = hwcart_split_type(domain[i]);
        if(split_type < 0){
            fprintf(stderr, "unknown memory domain %d, level %d\n", domain[i], i+1);
            hwloc_bitmap_free(m_cpuset);
            hwloc_bitmap_free(m_obj_cpuset);
            hwloc_topology_destroy (hwloctopo);
            return -1;
        }

        // clear the resource mask
        hwloc_bitmap_xor (m_obj_cpuset, m_obj_cpuset, m_obj_cpuset);
        
        int n = hwloc_get_nbobjs_by_type(hwloctopo, split_type);
        int ncomponents = 0;
        for(int j=0; j<n; j++){
          
            // figure out on which object we reside
            obj = hwloc_get_obj_by_type (hwloctopo, split_type, j);
            
            if(!obj){
                char name[256];
                hwcart_split_type_to_name(domain[i], name);
                fprintf(stderr, "no objects of type %s found on level %d\n", name, i);
                hwloc_bitmap_free(m_cpuset);
                hwloc_bitmap_free(m_obj_cpuset);
                hwloc_topology_destroy (hwloctopo);
                return -1;
            }
            
            if(hwloc_bitmap_isincluded (m_cpuset, obj->cpuset)){
                // include myself on this node
                hwloc_bitmap_or (m_obj_cpuset, m_obj_cpuset, obj->cpuset);
                level_rank_out[i] = j;
                break;
            }

            // m_cpuset might be larger than just one node: merge
            if(hwloc_bitmap_isincluded (obj->cpuset, m_cpuset)){
                ncomponents++;
                hwloc_bitmap_or (m_obj_cpuset, m_obj_cpuset, obj->cpuset);
            }

            if(hwloc_bitmap_isincluded (m_cpuset, m_obj_cpuset)){
                level_rank_out[i] = j/ncomponents;
                break;
            }
        }

        // find number of nodes on this level
        HWCART_MPI_CALL( MPI_Allreduce(level_rank_out+i, &nlevel_nodes, 1, MPI_INT, MPI_MAX, shmem_comm) );
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
        hwloc_topology_restrict(hwloctopo, m_obj_cpuset, 0);
    }

    // cleanup
    hwloc_bitmap_free(m_cpuset);
    hwloc_bitmap_free(m_obj_cpuset);
    hwloc_topology_destroy (hwloctopo);

    return 0;
}


// obtain level node ID of the calling rank
int hwcart_get_noderank(hwcart_topo_t hwtopo, MPI_Comm comm, hwcart_split_t in_split_type, int *noderank_out)
{
    int *sbuff, *rbuff;
    int nodeid, rank, size, ii;
    hwloc_obj_type_t split_type;    
    MPI_Comm shmem_comm;

    if(in_split_type == HWCART_MD_NODE){
    
        // old rank
        HWCART_MPI_CALL( MPI_Comm_rank(comm, &rank) );
        HWCART_MPI_CALL( MPI_Comm_size(comm, &size) );

        sbuff = calloc(size, sizeof(int));
        rbuff = calloc(size, sizeof(int));

        // create local communicator
        HWCART_MPI_CALL( MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm) );

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
  
    split_type = hwcart_split_type(in_split_type);
    if(split_type < 0){
        fprintf(stderr, "unknown memory domain %d\n", in_split_type);
        return -1;
    }

    hwloc_topology_t *phwloctopo = &hwtopo->topo;
    hwloc_bitmap_t m_obj_cpuset = hwloc_bitmap_alloc();
    hwloc_bitmap_t m_cpuset = hwloc_bitmap_alloc();
    hwloc_obj_t obj;
    hwloc_get_cpubind(*phwloctopo, m_cpuset, 0);

    int n = hwloc_get_nbobjs_by_type(*phwloctopo, split_type);
    int ncomponents = 0;
    for(int j=0; j<n; j++){

        // figure out on which object we reside
        obj = hwloc_get_obj_by_type (*phwloctopo, split_type, j);
        if(!obj){
            char name[256];
            hwcart_split_type_to_name(in_split_type, name);
            fprintf(stderr, "no objects of type %s found\n", name);
            hwloc_bitmap_free(m_cpuset);
            hwloc_bitmap_free(m_obj_cpuset);
            return -1;
        }
        
        if(hwloc_bitmap_isincluded (m_cpuset, obj->cpuset)){
            // include myself on this node
            hwloc_bitmap_or (m_obj_cpuset, m_obj_cpuset, obj->cpuset);
            *noderank_out = j;
            break;
        }

        // m_cpuset might be larger than just one node: merge
        if(hwloc_bitmap_isincluded (obj->cpuset, m_cpuset)){
            ncomponents++;
            hwloc_bitmap_or (m_obj_cpuset, m_obj_cpuset, obj->cpuset);
        }

        if(hwloc_bitmap_isincluded (m_cpuset, m_obj_cpuset)){
            *noderank_out = j/ncomponents;
            break;
        }
    }

    hwloc_bitmap_free(m_cpuset);
    hwloc_bitmap_free(m_obj_cpuset);
    return 0;
}


hwloc_obj_type_t hwcart_split_type(hwcart_split_t split_type)
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
