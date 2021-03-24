# hwcart

This library implements a hardware-aware Cartesian MPI communicator. 
A standard MPI Cartesian communicator created with `MPI_Cart_create`
allows the user to specify rank grid dimensions (`[nx, ny, nz]`, where
`np=nx*ny*nz` is the total number of ranks) and periodicity. How the ranks
are then distributed among the compute nodes, and how they are bound to
compute cores depends on the MPI runtime.

`hwcart` gives the user much finer control over how the ranks are
placed on various levels of the memory hierarchy. Consider 2-socket compute nodes
with AMD EPYC 7742 64-Core CPUs. 

```
$ hwloc-info
depth 0:            1 Machine (type #0)
 depth 1:           2 Package (type #1)
  depth 2:          8 Group0 (type #12)
   depth 3:         32 L3Cache (type #6)
    depth 4:        128 L2Cache (type #5)
     depth 5:       128 L1dCache (type #4)
      depth 6:      128 L1iCache (type #9)
       depth 7:     128 Core (type #2)
        depth 8:    256 PU (type #3)
Special depth -3:   8 NUMANode (type #13)

```
All cores have their private L1 and L2 caches, but there are 32 shared L3 caches. 
This means that 4 cores share an L3 cache, 16 cores (4 L3 cache groups) share 
the same NUMA node, and 4 NUMA-nodes are located on a single socket.

`hwcart` uses [`hwloc`](https://www.open-mpi.org/projects/hwloc/) to understand
the memory locality of the hardware and creates an equivalent of the MPI Cartesian
communicator, in which the neigoboring ranks are placed close to each other in
the complex memory hierarchy. This is done level by level, starting from the lowest
desired level. 

For installation instructions see [Installation](INSTALL.md).

## Example
Let's assume we want to start 128 ranks per compute node, on a total of
1024 compute nodes. First, you define the split granularity, i.e., which memory domain
levels are of interest. On the Epyc system this could be:

```
#define NLEVELS 5
    hwcart_split_t domain[NLEVELS] = {
                               HWCART_MD_CORE,
                               HWCART_MD_L3CACHE,
                               HWCART_MD_NUMA,
                               HWCART_MD_SOCKET,
                               HWCART_MD_NODE
    };
```
On the lowest level 4 cores share an L3 cache. Hence, we can define the rank
grid at this level to be `nx=4, ny=1, nz=1`, `nx=2, ny=2, nz=1`, `nx=2, ny=1, nz=2` or `nx=1, ny=2, nz=2`. 
The 4 L3 cache groups that belong to the same NUMA node should now be arranged themselves 
into a 3D grid. Again, the dimensions can be e.g., `nx=1, ny=2, nz=2`. On each higher level
the user specifies the grid dimension of the lower-level blocks. The final rak topology
may look as follows

```
    int topo[3*NLEVELS] = {
                           4, 1, 1, // core grid
                           1, 2, 2, // l3cache grid
                           1, 2, 2, // numa grid
                           1, 2, 1, // socket grid
                          16, 8, 8, // node grid
    };
```
In terms of a standard MPI Cartesian communicator this would correspond to a `[64,64,32]` grid
of ranks. Contrary to the MPI communicator, here it is possible to specify exactly how the
in-node rank grid should look to achieve best speed and minimize off-node communication.


## Example: rank placement
For the topology described in [Example](#example), here is how the rank grid looks on a single compute node, 
on the socket level:
```
Level 3 HWCART_MD_SOCKET

                      1   1   1   1
                    1   1   1   1
                  1   1   1   1
                1   1   1   1
              0   0   0   0
            0   0   0   0
          0   0   0   0
        0   0   0   0

                      1   1   1   1
                    1   1   1   1
                  1   1   1   1
                1   1   1   1
              0   0   0   0
            0   0   0   0
          0   0   0   0
        0   0   0   0

                      1   1   1   1
                    1   1   1   1
                  1   1   1   1
                1   1   1   1
              0   0   0   0
            0   0   0   0
          0   0   0   0
        0   0   0   0

                      1   1   1   1
                    1   1   1   1
                  1   1   1   1
                1   1   1   1
              0   0   0   0
            0   0   0   0
          0   0   0   0
        0   0   0   0

```
In the above rank grid, all ranks with the same value are executed on the same memory domain - in this case socket.
Here is how the distribution looks among NUMA nodes:
```
Level 2 HWCART_MD_NUMA

                      7   7   7   7
                    7   7   7   7
                  6   6   6   6
                6   6   6   6
              3   3   3   3
            3   3   3   3
          2   2   2   2
        2   2   2   2

                      7   7   7   7
                    7   7   7   7
                  6   6   6   6
                6   6   6   6
              3   3   3   3
            3   3   3   3
          2   2   2   2
        2   2   2   2

                      5   5   5   5
                    5   5   5   5
                  4   4   4   4
                4   4   4   4
              1   1   1   1
            1   1   1   1
          0   0   0   0
        0   0   0   0

                      5   5   5   5
                    5   5   5   5
                  4   4   4   4
                4   4   4   4
              1   1   1   1
            1   1   1   1
          0   0   0   0
        0   0   0   0
```
There are 8 NUMA nodes on this 2-socket system, and the printout above shows which ranks
are bound to which NUMA node.


## Rank numbering
Consider 4 ranks running within the same L3 cache domain, arranged into a grid of size `[2,2,1]`. 
These are two possible rank numberings:

```
2 3
0 1
```
or 
```
1 3
0 2
```
`hwcart` converts a 3-dimensional Cartesian coordinate of a rank into a linear rank id dimension by dimension.
The order in which it is done is defined by the user at initialization time. In the first example above 
the order is XY, in the second - YX.

For 3D grids 6 different order types are defined (see `hwcart_order_t` in `hwcart.h`), e.g., `HWCartOrderXYZ`, or `HWCartOrderZYX`. 
Note that while the order will only affect rank numbering and not placement, for real-world applications 
this can have a measurable impact on performance. For example, it may affect the order in which applications
perform certain computations, or communication.


## Cartesian communicator API
`hwcart` is initialized using
```
int hwcart_init(hwcart_topo_t *hwtopo_out);
```
This call returns a pointer to hardware information. When no longer needed, the hw topology information
can be freed by 
```
int hwcart_topo_free(hwcart_topo_t *hwtopo);
```

A `hwcart` communicator is created based on an existing MPI communicator (e.g., `MPI_COMM_WORLD`):
```
int hwcart_create(hwcart_topo_t hwtopo, MPI_Comm mpi_comm, int nlevels, hwcart_split_t *domain, int *topo, hwcart_order_t cart_order, MPI_Comm *hwcart_comm_out);
```
It can be freed when no longer needed:
```
int hwcart_comm_free(MPI_Comm *hwcart_comm);
```

Most importantly, `hwcart` provides functions to obtain a Cartesian index from the MPI rank:

```
hwcart_rank2coord(MPI_Comm hwcart_comm, int *dims, int rank, hwcart_order_t cart_order, int *coord_out);
```
and the other way round - to obtain MPI rank id from a set of Cartesian coordinates:

```
int hwcart_coord2rank(MPI_Comm hwcart_comm, int *dims, int *periodic, int *coord, hwcart_order_t cart_order, int *rank_out);
```
