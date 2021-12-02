MODULE hwcart_mod
  use iso_c_binding
  use mpi
  
  implicit none

  enum, bind(C)
     enumerator :: HWCART_MD_HWTHREAD = 0
     enumerator :: HWCART_MD_CORE
     enumerator :: HWCART_MD_L1CACHE
     enumerator :: HWCART_MD_L2CACHE
     enumerator :: HWCART_MD_L3CACHE
     enumerator :: HWCART_MD_SOCKET
     enumerator :: HWCART_MD_NUMA
     enumerator :: HWCART_MD_NODE
  end enum

  enum, bind(C)
     enumerator :: HWCartOrderXYZ = 0
     enumerator :: HWCartOrderXZY
     enumerator :: HWCartOrderZYX
     enumerator :: HWCartOrderYZX
     enumerator :: HWCartOrderZXY
     enumerator :: HWCartOrderYXZ
  end enum
  
  interface
     function hwcart_init() bind(C)
       use iso_c_binding
       integer :: hwcart_init
     end function hwcart_init

     function hwcart_finalize() bind(C)
       use iso_c_binding
       integer :: hwcart_finalize
     end function hwcart_finalize

     subroutine hwcart_comm_free(hwcart_comm) bind(C, name="hwcart_comm_free_f")
       use iso_c_binding
       integer :: hwcart_comm
     end subroutine hwcart_comm_free

     function hwcart_create_f(mpi_comm, nlevels, domain, topo, periodic, cart_order, hwcart_comm_out) bind(C)
       use iso_c_binding
       integer, value :: mpi_comm
       integer, value :: nlevels
       type(c_ptr), value :: domain
       type(c_ptr), value :: topo
       type(c_ptr), value :: periodic
       integer, value :: cart_order
       integer :: hwcart_comm_out
       integer :: hwcart_create_f
     end function hwcart_create_f

     function hwcart2mpicart_f(hwcart_comm, mpicart_comm_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       integer :: mpicart_comm_out
       integer :: hwcart2mpicart_f
     end function hwcart2mpicart_f
     
     function hwcart_sub_f(hwcart_comm, rank, belongs, hwcart_comm_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       integer, value :: rank
       type(c_ptr), value :: belongs
       integer :: hwcart_comm_out
       integer :: hwcart_sub_f
     end function hwcart_sub_f
     
     function hwcart_print_rank_topology_f(mpi_comm) bind(C)
       use iso_c_binding
       integer, value :: mpi_comm
       integer :: hwcart_print_rank_topology_f
     end function hwcart_print_rank_topology_f
     
     function hwcart_rank2coord_f(hwcart_comm, rank, coord_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       integer, value :: rank
       type(c_ptr), value :: coord_out
     end function hwcart_rank2coord_f
     
     function hwcart_coord2rank_f(hwcart_comm, coord, rank_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       type(c_ptr), value :: coord
       integer :: rank_out
     end function hwcart_coord2rank_f
     
  end interface
CONTAINS

  function hwcart_create(mpi_comm, domain, topo, periodic, cart_order, hwcart_comm_out)
    integer, value :: mpi_comm
    integer, dimension(:), target :: domain
    integer, dimension(:,:), target :: topo
    logical, dimension(:), target :: periodic
    integer, value :: cart_order
    integer :: hwcart_comm_out
    integer :: hwcart_create
    integer :: nlevels

    nlevels = size(domain)
    hwcart_create = &
      hwcart_create_f(mpi_comm, nlevels, c_loc(domain), c_loc(topo), c_loc(periodic), cart_order, hwcart_comm_out)
  end function hwcart_create
 
  function hwcart2mpicart(hwcart_comm, mpicart_comm_out)
    integer, value :: hwcart_comm
    integer :: mpicart_comm_out
    integer :: hwcart2mpicart

    hwcart2mpicart = &
      hwcart2mpicart_f(hwcart_comm, mpicart_comm_out)    
  end function hwcart2mpicart

  function hwcart_sub(hwcart_comm, rank, belongs, hwcart_comm_out)
    integer, value :: hwcart_comm
    integer, value :: rank
    logical, dimension(:), target :: belongs(3)
    integer :: hwcart_comm_out
    integer :: hwcart_sub

    hwcart_sub = &
      hwcart_sub_f(hwcart_comm, rank, c_loc(belongs), hwcart_comm_out)
  end function hwcart_sub

  function hwcart_print_rank_topology(hwcart_comm)
    integer, value :: hwcart_comm
    integer :: hwcart_print_rank_topology
    
    hwcart_print_rank_topology = &
      hwcart_print_rank_topology_f(hwcart_comm)
  end function hwcart_print_rank_topology

  function hwcart_rank2coord(hwcart_comm, rank, coord_out)
    use iso_c_binding
    integer, value :: hwcart_comm
    integer, value :: rank
    integer, dimension(:), target :: coord_out(3)
    integer :: hwcart_rank2coord

    hwcart_rank2coord = &
      hwcart_rank2coord_f(hwcart_comm, rank, c_loc(coord_out))
  end function hwcart_rank2coord
 
  function hwcart_coord2rank(hwcart_comm, coord, rank_out)
    use iso_c_binding
    integer, value :: hwcart_comm
    integer, dimension(:), target :: coord(3)
    integer :: rank_out
    integer :: hwcart_coord2rank

    hwcart_coord2rank = &
      hwcart_coord2rank_f(hwcart_comm, c_loc(coord), rank_out)
  end function hwcart_coord2rank
 
END MODULE hwcart_mod
