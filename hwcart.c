#include "hwcart.h"
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <string.h>


int  hwcart_ompi_split_type(int split_type);
void hwcart_split_type_to_name(int split_type, char *name);
void hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, int order);


int hwcart_topology(MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank_out, int level)
{
    int *sbuff, *rbuff;
    int ierr, comm_rank, comm_size, nodeid, noderank, color, ii;
    MPI_Comm split_comm;
    int split_type, split_rank, split_size;
    MPI_Comm master_comm;
    int master_size;
    int retval = 0;

    // parent communicator
    ierr = MPI_Comm_rank(comm, &comm_rank);
    ierr = MPI_Comm_size(comm, &comm_size);

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

    split_type = hwcart_ompi_split_type(domain[level]);
    if(split_type < 0){
      fprintf(stderr, "unknown memory domain %d, level %d\n", domain[level], level);
      return -1;
    }

    // create communicator for this topology level
    ierr = MPI_Comm_split_type(comm, split_type, 0, MPI_INFO_NULL, &split_comm);
    ierr = MPI_Comm_rank(split_comm, &split_rank);
    ierr = MPI_Comm_size(split_comm, &split_size);

    // no split on this topology level
    if (split_size == comm_size) {
        if (1 != topo[level*3+0]*topo[level*3+1]*topo[level*3+2]){
            fprintf(stderr, "ERROR (2): wrong topology on level %d: expected %d domains (config), but found %d (hardware)\n",
                    level, topo[level*3+0]*topo[level*3+1]*topo[level*3+2], 1);
            return -1;
        }
        level_rank_out[level] = 0;
        return hwcart_topology(split_comm, nsplits, domain, topo, level_rank_out, level-1);
    }

    // make a master-rank communicator: masters from each split comm join
    color = 0;
    if (split_rank != 0){

        // non-masters
        ierr = MPI_Comm_split(comm, color, 0, &master_comm);
    } else {

        // masters
        // temporary nodeid identifier: rank of the split master, +1 needed
        nodeid = comm_rank+1;
        color = 1;
        ierr = MPI_Comm_split(comm, color, 0, &master_comm);
        ierr = MPI_Comm_size(master_comm, &master_size);

        // verify topology validity on this level
        if (master_size != topo[level*3+0]*topo[level*3+1]*topo[level*3+2]){
            fprintf(stderr, "ERROR (3): wrong topology on level %d: expected %d domains (config), but found %d (hardware)\n",
                    level, topo[level*3+0]*topo[level*3+1]*topo[level*3+2], master_size);
            return -1;
        }

        // comm buffers to establish unique node id's for each master
        sbuff = (int*)calloc(comm_size,sizeof(int));
        rbuff = (int*)calloc(comm_size,sizeof(int));

        // find node rank based on unique node id from above
        for(int ii=0; ii<comm_size; ii++) sbuff[ii] = nodeid;
        ierr = MPI_Alltoall(sbuff, 1, MPI_INT, rbuff, 1, MPI_INT, master_comm);

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

    ierr = MPI_Comm_disconnect(&master_comm);

    // distribute the node id to all split ranks
    ierr = MPI_Bcast(&noderank, 1, MPI_INT, 0, split_comm);

    // save our level rank
    level_rank_out[level] = noderank;

    // sub-divide lower levels
    retval = hwcart_topology(split_comm, nsplits, domain, topo, level_rank_out, level-1);

    // cleanup
    ierr = MPI_Comm_disconnect(&split_comm);

    return retval;
}


int hwcart_remap_ranks(MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank, int order, MPI_Comm *hwcart_comm_out)
{
  int topo_coord[3] = {0}, gdim[3] = {1,1,1}, periodic[3] = {0};
  int ierr, ii, comm_rank, newrank;
  int *cartXYZ;
  int *dims;

  ierr = MPI_Comm_rank(comm, &comm_rank);

  cartXYZ = calloc(3*nsplits, sizeof(int));
  dims    = calloc(3*nsplits, sizeof(int));

  // compute rank cartesian coordinates for each topological level
  for(ii=0; ii<nsplits; ii++){
    hwcart_rank2coord(MPI_COMM_NULL, topo+3*ii, level_rank[ii], order, cartXYZ+3*ii);
  }

  // compute rank grid resulution for each level
  dims[0] = 1;
  dims[1] = 1;
  dims[2] = 1;
  for(ii=1; ii<nsplits; ii++){
    dims[ii*3+0] = dims[(ii-1)*3+0]*topo[(ii-1)*3+0];
    dims[ii*3+1] = dims[(ii-1)*3+1]*topo[(ii-1)*3+1];
    dims[ii*3+2] = dims[(ii-1)*3+2]*topo[(ii-1)*3+2];
  }

  // compute resulting global cartesian index
  for(ii=0; ii<nsplits; ii++){
    topo_coord[0] = topo_coord[0] + cartXYZ[ii*3+0]*dims[ii*3+0];
    topo_coord[1] = topo_coord[1] + cartXYZ[ii*3+1]*dims[ii*3+1];
    topo_coord[2] = topo_coord[2] + cartXYZ[ii*3+2]*dims[ii*3+2];
  }

  // compute global grid dimensions
  for(ii=0; ii<nsplits; ii++){
    gdim[0] *= topo[ii*3+0];
    gdim[1] *= topo[ii*3+1];
    gdim[2] *= topo[ii*3+2];
  }

  // compute rank id in global rank space
  hwcart_coord2rank(comm, gdim, periodic, topo_coord, order, &newrank);

  // create the new communicator with remapped ranks
  ierr = MPI_Comm_split(comm, 0, newrank, hwcart_comm_out);

  // cleanup
  free(dims);
  free(cartXYZ);

  return 0;
}


int  hwcart_create(MPI_Comm mpi_comm, int nsplits, int *domain, int *topo, int order, MPI_Comm *hwcart_comm_out)
{
  int retval, cleanup = 0;
  int *level_rank;

  level_rank = calloc(nsplits, sizeof(int));
  retval = hwcart_topology(mpi_comm, nsplits, domain, topo, level_rank, nsplits-1);
  if(retval<0) {
    free(level_rank);
    return retval;
  }

  retval = hwcart_remap_ranks(mpi_comm, nsplits, domain, topo, level_rank, order, hwcart_comm_out);
  free(level_rank);
  return retval;
}


int hwcart_free(MPI_Comm *hwcart_comm)
{
  return MPI_Comm_free(hwcart_comm);
}


int  hwcart_rank2coord(MPI_Comm comm, int *dims, int rank, int order, int *coord)
{
  int topo, ierr, tmp;

  // check if this is a cartesian communicator
  if (comm == MPI_COMM_NULL) {
    topo = -1;
  } else {
    ierr =  MPI_Topo_test(comm, &topo);
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
    ierr = MPI_Cart_coords(comm, rank, 3, coord);
  }
  return 0;
}


int  hwcart_coord2rank(MPI_Comm comm, int *dims, int *periodic, int *coord, int order, int *rank_out)
{
  int tcoord[3], ii, topo, ierr;

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
    ierr =  MPI_Topo_test(comm, &topo);
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
    ierr = MPI_Cart_rank(comm, tcoord, rank_out);
  }
  return 0;
}


// obtain level node ID of the calling rank
int ghex_get_noderank(MPI_Comm comm, int split_type, int *noderank_out)
{
  int *sbuff, *rbuff;
  int nodeid, rank, size, ierr, ii;
  MPI_Comm shmcomm;
  char version[MPI_MAX_LIBRARY_VERSION_STRING];
  int resultlen;

  // old rank
  ierr = MPI_Comm_rank(comm, &rank);
  ierr = MPI_Comm_size(comm, &size);

  // communication buffers
  sbuff = calloc(size, sizeof(int));
  rbuff = calloc(size, sizeof(int));

  //check for OpenMPI to use fine-grain hwloc domains
  ierr = MPI_Get_library_version(version, &resultlen);
  if(0 == strncmp(version, "Open MPI", strlen("Open MPI"))){

    // create local communicator
    split_type = hwcart_ompi_split_type(split_type);
    if(split_type < 0){
      fprintf(stderr, "unknown memory domain %d\n", split_type);
      return -1;
    }
    ierr = MPI_Comm_split_type(comm, split_type, 0, MPI_INFO_NULL, &shmcomm);
  } else {

    // create local communicator
    ierr = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
  }

  // figure out unique compute node id = max rank id on each compute node
  sbuff[0] = rank+1;
  ierr = MPI_Allreduce(sbuff, rbuff, 1, MPI_INT, MPI_MAX, shmcomm);
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
  ierr = MPI_Comm_disconnect(&shmcomm);

  return 0;
}


int hwcart_print_rank_topology(MPI_Comm comm, int nlevels, int *domain, int *topo, int order)
{
  int *buff;
  int ierr, ii, li, gdim[3] = {1,1,1};
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
    ghex_get_noderank(comm, domain[ii], level_id+ii);
  }

  sbuff = calloc(nlevels+3, sizeof(int));

  // obtain all values at master
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &orank);
  ierr = MPI_Comm_rank(comm, &comm_rank);
  ierr = MPI_Comm_size(comm, &comm_size);

  buff = calloc((nlevels+3)*comm_size, sizeof(int));
  sbuff[0] = comm_rank;
  sbuff[1] = orank;
  sbuff[2] = sched_getcpu(); // -D_GNU_SOURCE
  for(ii=0; ii<nlevels; ii++) sbuff[3+ii] = level_id[ii];

  // assuming HT is enabled, use the lower core ID
  ncpus = get_nprocs_conf()/2;
  if (sbuff[2] >= ncpus) sbuff[2] = sbuff[2] - ncpus;
  ierr = MPI_Gather(sbuff, nlevels+3, MPI_INT, buff, nlevels+3, MPI_INT, 0, comm);

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
    printf("Rank layout\n");
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


void hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, int order){
  int ierr, k, j, i, kk, n;
  int comm_size;
  int periodic[3] = {0}, idx[3];

  ierr = MPI_Comm_size(comm, &comm_size);

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
	  printf("%3d", buff[n*line_size + id]);
      }
      printf("\n");
    }
    printf("\n");
  }
}


void hwcart_split_type_to_name(int split_type, char *name) {

  switch (split_type) {
  case (MPI_COMM_TYPE_SHARED):
    sprintf(name, "MPI_COMM_TYPE_SHARED");
    break;
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


int hwcart_ompi_split_type(int split_type)
{
  switch (split_type) {
  case (MPI_COMM_TYPE_SHARED):
    return MPI_COMM_TYPE_SHARED;
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
