PROGRAM hwcart_test_f
  use mpi
  use hwcart_mod

  implicit none

  type(hwcart_topo_t) :: hwcart_topo
  integer :: hwcart_comm, hwcart_comm_col, hwcart_comm_plane
  integer :: ierr

  integer, parameter :: NLEVELS=5
  integer, dimension(:), target :: domain(NLEVELS) = [   &
    HWCART_MD_CORE,   &
    HWCART_MD_L3CACHE,&
    HWCART_MD_NUMA,   &
    HWCART_MD_SOCKET, &
    HWCART_MD_NODE]

  integer, dimension(:,:), target :: topo(3,NLEVELS) = reshape([&
    4, 1, 1, & ! core grid 
    1, 2, 2, & ! l3cache grid
    1, 2, 2, & ! numa grid
    1, 2, 1, & ! socket grid
    1, 1, 1],& ! node grid
    shape(topo))

  integer :: cart_order = HWCartOrderZYX
  integer :: hwcart_rank, nbleft, nbright
  integer, dimension(:) :: gdim(3)
  logical, dimension(:) :: periodic(3) = [.true.,.true.,.true.]
  integer, dimension(:) :: coord(3)

  ! global rank dimensions
  gdim = product(topo, 2)

  call MPI_Init(ierr);

  ierr = hwcart_init(hwcart_topo)
  if (ierr /= 0) then
    call MPI_Finalize(ierr);
    call exit(1)
  end if

  ierr = hwcart_create(hwcart_topo, MPI_COMM_WORLD, domain, topo, cart_order, hwcart_comm)
  if (ierr == 0) then
    ierr = hwcart_print_rank_topology(hwcart_topo, hwcart_comm, domain, topo, cart_order);
    call MPI_Barrier(hwcart_comm, ierr)

    ! get our cartesian coordinates
    call MPI_Comm_rank(hwcart_comm, hwcart_rank, ierr)
    ierr = hwcart_rank2coord(hwcart_comm, gdim, hwcart_rank, cart_order, coord)

    ! query left and right neighbors (in X dimension)
    coord(1) = coord(1)-1
    ierr = hwcart_coord2rank(hwcart_comm, gdim, periodic, coord, cart_order, nbleft)
    coord(1) = coord(1)+2
    ierr = hwcart_coord2rank(hwcart_comm, gdim, periodic, coord, cart_order, nbright)
    coord(1) = coord(1)-1
    ! print *, "rank", hwcart_rank, "coord", coord, "left", nbleft, "right", nbright

    block
      logical, dimension(:) :: belongs(3)
      integer :: col_size, plane_size
      
      ! create column communicator
      belongs(:) = (/.false.,.false.,.true./)
      ierr = hwcart_sub(hwcart_comm, gdim, hwcart_rank, cart_order, belongs, hwcart_comm_col)
      call MPI_Comm_size(hwcart_comm_col, col_size, ierr)
      if (col_size /= gdim(3)) then
         write (*,*) "ERROR: wrong size of the column communicator: ", col_size, "/=", gdim(3)
         call exit(1)
      end if

      ! create plane communicator
      belongs(:) = (/.true.,.true.,.false./)
      ierr = hwcart_sub(hwcart_comm, gdim, hwcart_rank, cart_order, belongs, hwcart_comm_plane)
      call MPI_Comm_size(hwcart_comm_plane, plane_size, ierr)
      if (plane_size /= gdim(1)*gdim(2)) then
         write (*,*) "ERROR: wrong size of the plane communicator: ", plane_size, "/=", gdim(1)*gdim(2)
         call exit(1)
      end if
      
    end block
    
    call hwcart_comm_free(hwcart_comm_col)
    call hwcart_comm_free(hwcart_comm_plane)
    call hwcart_comm_free(hwcart_comm)
    call hwcart_topo_free(hwcart_topo)
  else
     call exit(1)
  end if
  
  call MPI_Finalize(ierr);
END PROGRAM hwcart_test_f
