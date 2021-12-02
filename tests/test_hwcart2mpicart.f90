PROGRAM hwcart_test_hwcart2mpicart_f
  use mpi
  use hwcart_mod

  implicit none

  type(hwcart_topo_t) :: hwcart_topo
  integer :: hwcart_comm, mpicart_comm
  integer :: ierr

  integer, parameter :: NLEVELS=5
  integer, dimension(:), target :: domain(NLEVELS) = [   &
    HWCART_MD_CORE,   &
    HWCART_MD_L3CACHE,&
    HWCART_MD_NUMA,   &
    HWCART_MD_SOCKET, &
    HWCART_MD_NODE]

  integer, dimension(:,:), target :: topo(3,NLEVELS) = reshape([&
    2, 2, 1, & ! core grid 
    1, 1, 1, & ! l3cache grid
    1, 1, 1, & ! numa grid
    1, 1, 1, & ! socket grid
    1, 1, 1],& ! node grid
    shape(topo))

  integer :: cart_order = HWCartOrderXYZ
  integer, dimension(:) :: gdim(3)
  logical, dimension(:) :: periodic(3) = [.true.,.true.,.true.]

  integer :: hwcart_rank, mpicart_rank
  integer, dimension(:) :: hwcart_coord(3), mpicart_coord(3)

  ! global rank dimensions
  gdim = product(topo, 2)

  call MPI_Init(ierr);

  ierr = hwcart_init(hwcart_topo)
  if (ierr /= 0) then
    call MPI_Finalize(ierr);
    call exit(1)
  end if

  ierr = hwcart_create(hwcart_topo, MPI_COMM_WORLD, domain, topo, periodic, cart_order, hwcart_comm)
  if (ierr == 0) then
    ierr = hwcart_print_rank_topology(hwcart_topo, hwcart_comm);

    ! Convert a HWCART communicator to a native MPI_Cart communicator.
    ! Dimension order will be set to MPI order (ZYX). All standard MPI_Cart functions
    ! (MPI_Cart_shift, MPI_Cart_sub, MPI_Cart_rank, etc.) can be used on mpicart_com.
    ierr = hwcart2mpicart(hwcart_comm, mpicart_comm)

    ! hwcart cartesian coordinates
    call MPI_Comm_rank(hwcart_comm, hwcart_rank, ierr)
    ierr = hwcart_rank2coord(hwcart_comm, hwcart_rank, hwcart_coord)
    
    ! mpicart cartesian coordinates
    call MPI_Comm_rank(mpicart_comm, mpicart_rank, ierr);
    call MPI_Cart_coords(mpicart_comm, mpicart_rank, 3, mpicart_coord, ierr);

    if (any(hwcart_coord /= mpicart_coord)) then
       write (*,*) "ERROR: inconsistent MPI_Cart and hwcart coordinates"
       call exit(1)       
    end if

    call MPI_Comm_free(mpicart_comm, ierr)
    call hwcart_comm_free(hwcart_comm)
    call hwcart_topo_free(hwcart_topo)
  else
     call exit(1)
  end if
  
  call MPI_Finalize(ierr);
END PROGRAM hwcart_test_hwcart2mpicart_f
