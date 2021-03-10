PROGRAM test_halo_exchange
  use mpi
  use hwcart_mod

  implicit none

  type(hwcart_topo_t) :: hwcart_topo
  integer :: hwcart_comm
  integer :: ierr

  ! integer, parameter :: NLEVELS=5
  ! integer, dimension(:), target :: domain(NLEVELS) = [   &
  !   HWCART_MD_CORE,   &
  !   HWCART_MD_L3CACHE,&
  !   HWCART_MD_NUMA,   &
  !   HWCART_MD_SOCKET, &
  !   HWCART_MD_NODE]

  ! integer, dimension(:), target :: topo(3*NLEVELS) = [&
  !   4, 1, 1, & ! core grid 
  !   1, 2, 2, & ! l3cache grid
  !   1, 2, 2, & ! numa grid
  !   1, 2, 1, & ! socket grid
  !   1, 1, 1]   ! node grid

  integer, parameter :: NLEVELS=2
  integer, dimension(:), target :: domain(NLEVELS) = [   &
    HWCART_MD_CORE,   &
    HWCART_MD_NODE]

  integer, dimension(:), target :: topo(3*NLEVELS) = [&
    2, 2, 1, & ! core grid 
    1, 1, 1]   ! node grid

  integer :: order = HWCartOrderZYX
  integer :: hwcart_rank, nbleft, nbright
  integer, dimension(:) :: gdim(3)
  integer, dimension(:) :: periodic(3) = [1, 1, 1]
  integer, dimension(:) :: coord(3)

  ! global rank dimensions
  gdim = product(reshape(topo, [3,NLEVELS]), 2)

  call MPI_Init(ierr);

  ierr = hwcart_init(hwcart_topo)
  if (ierr /= 0) then
    call MPI_Finalize(ierr);
    call exit(1)
  end if
  
  ierr = hwcart_create(hwcart_topo, MPI_COMM_WORLD, NLEVELS, domain, topo, order, hwcart_comm)
  if (ierr == 0) then
    ierr = hwcart_print_rank_topology(hwcart_topo, hwcart_comm, NLEVELS, domain, topo, order);
    call MPI_Barrier(hwcart_comm, ierr)

    ! get our cartesian coordinates
    call MPI_Comm_rank(hwcart_comm, hwcart_rank, ierr)
    ierr = hwcart_rank2coord(hwcart_comm, gdim, hwcart_rank, order, coord)

    ! query left and right neighbors (in X dimension)
    coord(1) = coord(1)-1
    ierr = hwcart_coord2rank(hwcart_comm, gdim, periodic, coord, order, nbleft)
    coord(1) = coord(1)+2
    ierr = hwcart_coord2rank(hwcart_comm, gdim, periodic, coord, order, nbright)
    coord(1) = coord(1)-1

    print *, "rank", hwcart_rank, "coord", coord, "left", nbleft, "right", nbright
    call hwcart_free(hwcart_topo, hwcart_comm);
  end if
  
  call MPI_Finalize(ierr);
END PROGRAM test_halo_exchange
