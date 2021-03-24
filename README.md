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
the complex memory hierarchy. This is done level by level, starring from the lowest
desired level. Let's assume we want to start 128 ranks per compute node, on a total of
1024 compute nodes. First, you define the split granularity, i.e., which memory domain
levels are of interest:

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
On the lowest level 4 cores share an L3 cache. Hence, we can define the core
grid at the lowest level to be `nx=4, ny=1, nz=1`, or `nx=2, ny=1, nz=1`. 
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

## Resulting rank placement
Looking at the above example, here is how the rank grid looks on a single compute node, on the socket level:
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

## Cartesian communicator API
`hwcart` provides functions to obtain a Cartesian index from the MPI rank:

```
hwcart_rank2coord(MPI_Comm hwcart_comm, int *dims, int rank, hwcart_order_t cart_order, int *coord_out);
```
and the other way round - to obtain MPI rank ID from a set of Cartesian coordinates:

```
int hwcart_coord2rank(MPI_Comm hwcart_comm, int *dims, int *periodic, int *coord, hwcart_order_t cart_order, int *rank_out);
```
In the above, `hwcart_order_t` defines the significance of the three dimensions when computing the cartesian indices, e.g., 
`HWCartOrderXYZ`, or `HWCartOrderZYX`. The order will affect the rank numbering only, not placement. However, 
our experiments show that for real-world applications this can have a measurable impact on performance,
as it may affect the order in which applications perform certain computations, or communication.
