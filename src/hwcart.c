#include <hwcart/hwcart.h>
#include <hwcart/hwcart_utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#ifdef HWCART_USE_HWLOC
#include <hwloc.h>
#include <hwloc/shmem.h>

hwloc_obj_type_t hwcart_split_type(hwcart_split_t split_type);

struct hwcart_topo_struct_t {
    hwloc_topology_t topo;
};


int load_hwtopo(hwcart_topo_t hwtopo)
{
    int res;
    res = hwloc_topology_init(&hwtopo->topo);
    if(res) {
        fprintf(stderr, "failed to initialize topology data structure\n");
        return res;
    }
    res = hwloc_topology_load(hwtopo->topo);
    if(res) {
        fprintf(stderr, "failed to load topology\n");
        hwloc_topology_destroy (hwtopo->topo);
        return res;
    }
    return 0;
}


int hwcart_init(hwcart_topo_t *hwtopo_out)
{
    int res;
    pid_t master_pid;
    size_t shmem_size;
    void *shmem_addr;
    char shmem_filename[256];
    MPI_Comm shmem_comm;   
    int rank;

    *hwtopo_out = malloc(sizeof(struct hwcart_topo_struct_t));

    // a single process reads the topo for each compute node
    // the topo is then shared through shared memory
    HWCART_MPI_CALL( MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmem_comm) );
    HWCART_MPI_CALL( MPI_Comm_rank(shmem_comm, &rank) );

    if(rank==0){

        res = load_hwtopo(*hwtopo_out);
        if(0 != res){
            free(*hwtopo_out);
            *hwtopo_out = NULL;
            return res;
        }

        hwloc_shmem_topology_get_length((*hwtopo_out)->topo, &shmem_size, 0);

        master_pid = getpid();
        snprintf(shmem_filename, 256, "/dev/shm/hwcart_topo.sm.%i", master_pid);
        int shmem_fd = open(shmem_filename, O_CREAT | O_RDWR, 0600);
        if(-1 == shmem_fd){
            fprintf(stderr, "failed to create backing file in shared memory.\n");

            // each process loads the topology itself
            master_pid = 0;
            HWCART_MPI_CALL( MPI_Bcast(&master_pid, sizeof(master_pid), MPI_BYTE, 0, shmem_comm) );
            return 0;
        }

        // use mmap to get a default address
        shmem_addr = mmap(NULL, shmem_size, PROT_READ | PROT_WRITE, MAP_SHARED, shmem_fd, 0);
        if(NULL == shmem_addr){
            fprintf(stderr, "failed to mmap shared memory\n");
            close(shmem_fd);
            unlink(shmem_filename);

            // each process loads the topology itself
            master_pid = 0;
            HWCART_MPI_CALL( MPI_Bcast(&master_pid, sizeof(master_pid), MPI_BYTE, 0, shmem_comm) );
            return 0;
        }

        // free the address
        munmap(shmem_addr, shmem_size);

        // store the topology
        res = hwloc_shmem_topology_write((*hwtopo_out)->topo, shmem_fd, 0, shmem_addr, shmem_size, 0);
        if(0 != res){
            fprintf(stderr, "failed to write topology info to shared memory\n");
            close(shmem_fd);
            unlink(shmem_filename);

            // each process loads the topology itself
            master_pid = 0;
            HWCART_MPI_CALL( MPI_Bcast(&master_pid, sizeof(master_pid), MPI_BYTE, 0, shmem_comm) );
            return 0;
        }
        close(shmem_fd);
        
        // broadcast topology info
        HWCART_MPI_CALL( MPI_Bcast(&master_pid, sizeof(master_pid), MPI_BYTE, 0, shmem_comm) );
        HWCART_MPI_CALL( MPI_Bcast(&shmem_size, sizeof(shmem_size), MPI_BYTE, 0, shmem_comm) );
        HWCART_MPI_CALL( MPI_Bcast(&shmem_addr, sizeof(shmem_addr), MPI_BYTE, 0, shmem_comm) );

        MPI_Barrier(shmem_comm);
        unlink(shmem_filename);

    } else {

        // do we have shared topology?
        HWCART_MPI_CALL( MPI_Bcast(&master_pid, sizeof(master_pid), MPI_BYTE, 0, shmem_comm) );
        if(0 != master_pid){

            HWCART_MPI_CALL( MPI_Bcast(&shmem_size, sizeof(shmem_size), MPI_BYTE, 0, shmem_comm) );
            HWCART_MPI_CALL( MPI_Bcast(&shmem_addr, sizeof(shmem_addr), MPI_BYTE, 0, shmem_comm) );

            snprintf(shmem_filename, 256, "/dev/shm/hwcart_topo.sm.%i", master_pid);
            int shmem_fd = open(shmem_filename, O_RDONLY, 0600);

            // we can clean up the shm now
            MPI_Barrier(shmem_comm);
            
            if(-1 == shmem_fd){
                fprintf(stderr, "failed to open backing file in shared memory.\n");

                // load the topology
                res = load_hwtopo(*hwtopo_out);
                if(0 != res){
                    free(*hwtopo_out);
                    *hwtopo_out = NULL;
                    return res;
                }
                return 0;
            }
            
            res = hwloc_shmem_topology_adopt(&(*hwtopo_out)->topo, shmem_fd, 0, shmem_addr, shmem_size, 0);
            if(0 != res){
                // fprintf(stderr, "failed to adopt topology, reading topology directly\n");
                // load the topology
                res = load_hwtopo(*hwtopo_out);
                if(0 != res){
                    free(*hwtopo_out);
                    *hwtopo_out = NULL;
                    return res;
                }
                return 0;
            }
            close(shmem_fd);
	    return 0;
        }
        
        // load the topology
        res = load_hwtopo(*hwtopo_out);
        if(0 != res){
            free(*hwtopo_out);
            *hwtopo_out = NULL;
            return res;
        }
    }
    return 0;
}


int  hwcart_topo_free(hwcart_topo_t *hwtopo)
{
    if(hwtopo && (*hwtopo)){
        hwloc_topology_destroy((*hwtopo)->topo);
        free(*hwtopo);
        *hwtopo = NULL;
    }
    return 0;
}


int hwcart_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, int *level_rank_out, int level)
{
    int retval = 0;
    int comm_rank, comm_size;
    hwloc_obj_type_t split_type;
    MPI_Comm parent_comm;
    int shmem_rank;
    MPI_Comm level_comm;
    int level_rank;
    MPI_Comm master_comm;
    int nlevel_nodes;

    // init
    for(int i=0; i<nlevels; i++) level_rank_out[i] = -1;

    // parent communicator
    HWCART_MPI_CALL( MPI_Comm_rank(comm, &comm_rank) );
    HWCART_MPI_CALL( MPI_Comm_size(comm, &comm_size) );

    // create top-level parent communicator for shared memory nodes
    HWCART_MPI_CALL( MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &parent_comm) );
    HWCART_MPI_CALL( MPI_Comm_rank(parent_comm, &shmem_rank) );

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
    HWCART_MPI_CALL( MPI_Bcast(level_rank_out+level, 1, MPI_INT, 0, parent_comm) );

    // we need a copy of the topology: we modify it
    hwloc_topology_t hwloctopo;
    hwloc_topology_dup(&hwloctopo, hwtopo->topo);
    
    hwloc_bitmap_t m_obj_cpuset = hwloc_bitmap_alloc();
    hwloc_bitmap_t m_cpuset = hwloc_bitmap_alloc();
    hwloc_get_cpubind(hwloctopo, m_cpuset, 0);

    for(int i=level-1; i>=0; i--){
        split_type = hwcart_split_type(domain[i]);
        if(split_type < 0){
            fprintf(stderr, "unknown memory domain %d, level %d\n", domain[i], i+1);
	    retval = -1;
	    goto cleanup;
        }

        // clear the resource mask
        hwloc_bitmap_xor (m_obj_cpuset, m_obj_cpuset, m_obj_cpuset);
        
        int n = hwloc_get_nbobjs_by_type(hwloctopo, split_type);
        int ncomponents = 1;

        for(int j=0; j<n; j++){
	    
            // figure out on which object we reside
            hwloc_obj_t obj = hwloc_get_obj_by_type (hwloctopo, split_type, j);
          
            if(!obj){
                char name[256];
                hwcart_split_type_to_name(domain[i], name);
                fprintf(stderr, "no objects of type %s found on level %d\n", name, i);
		retval = -1;
		goto cleanup;
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

            if(hwloc_bitmap_isequal (m_cpuset, m_obj_cpuset)){
                
		// include myself on this node
                level_rank_out[i] = j;
		break;
            }
        }

	if(level_rank_out[i] == -1){
	    fprintf(stderr, "rank %d was not included in any memory domains on level %d\n", comm_rank, i);
	    retval = -1;
	    goto cleanup;
	}

	// create communicator for ranks sharing the same node on this level
	HWCART_MPI_CALL( MPI_Comm_split(parent_comm, level_rank_out[i], 0, &level_comm) );
	HWCART_MPI_CALL( MPI_Comm_rank(level_comm, &level_rank) );

	// make a master-rank communicator: masters from each split comm join
	nlevel_nodes=-1;
	if (level_rank != 0){

	    // non-masters
	    HWCART_MPI_CALL( MPI_Comm_split(parent_comm, 0, 0, &master_comm) );
	} else {

	    // masters
	    HWCART_MPI_CALL( MPI_Comm_split(parent_comm, 1, 0, &master_comm) );
	    HWCART_MPI_CALL( MPI_Comm_rank(master_comm, level_rank_out+i) );
	    HWCART_MPI_CALL( MPI_Comm_size(master_comm, &nlevel_nodes) );
	}
    
	// cleanup
	HWCART_MPI_CALL( MPI_Comm_disconnect(&master_comm) );
    
	// distribute the node id to all split ranks
	HWCART_MPI_CALL( MPI_Bcast(level_rank_out+i, 1, MPI_INT, 0, level_comm) );
	HWCART_MPI_CALL( MPI_Bcast(&nlevel_nodes, 1, MPI_INT, 0, level_comm) );
	HWCART_MPI_CALL( MPI_Comm_disconnect(&parent_comm) );
	parent_comm = level_comm;
	
        if (nlevel_nodes != topo[i*3+0]*topo[i*3+1]*topo[i*3+2]){
            
            // those either have to match, or there is no split on this level
            
            // no split
            if((i != 0) && (1 == topo[i*3+0]*topo[i*3+1]*topo[i*3+2])) {
                level_rank_out[i] = 0;
                continue;
            }
	    
	    if(level_rank==0)
		fprintf(stderr, "ERROR: wrong topology on level %d: expected %d nodes on this level, instead found %d\n",
			i, topo[i*3+0]*topo[i*3+1]*topo[i*3+2], nlevel_nodes);
	    retval = -1;
	    goto cleanup;
        }

        // create a sub-topology
        hwloc_topology_restrict(hwloctopo, m_obj_cpuset, 0);
    }

    // cleanup
 cleanup:
    HWCART_MPI_CALL( MPI_Comm_disconnect(&parent_comm) );
    hwloc_bitmap_free(m_cpuset);
    hwloc_bitmap_free(m_obj_cpuset);
    hwloc_topology_destroy (hwloctopo);

    return retval;
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
            // NOTE: this works under the assumption that ranks are assigned 
            // to subsequently numbered resources
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
