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
  
  type, bind(c) :: hwcart_topo_t
     type(c_ptr) :: ptr = c_null_ptr
  end type hwcart_topo_t
  
  interface
     function hwcart_init(hwcart_topo) bind(C)
       use iso_c_binding
       import hwcart_topo_t
       integer :: hwcart_init
       type(hwcart_topo_t), intent(out) :: hwcart_topo
     end function hwcart_init

     subroutine hwcart_free(hwcart_topo, hwcart_comm) bind(C, name="hwcart_free_f")
       use iso_c_binding
       import hwcart_topo_t
       type(hwcart_topo_t) :: hwcart_topo
       integer :: hwcart_comm
     end subroutine hwcart_free

     function hwcart_create_f(hwtopo, mpi_comm, nlevels, domain, topo, cart_order, hwcart_comm_out) bind(C)
       use iso_c_binding
       import hwcart_topo_t
       type(hwcart_topo_t), value :: hwtopo
       integer, value :: mpi_comm
       integer, value :: nlevels
       type(c_ptr), value :: domain
       type(c_ptr), value :: topo
       integer, value :: cart_order
       integer :: hwcart_comm_out
       integer :: hwcart_create_f
     end function hwcart_create_f
     
     function hwcart_print_rank_topology_f(hwtopo, mpi_comm, nlevels, domain, topo, cart_order) bind(C)
       use iso_c_binding
       import hwcart_topo_t
       type(hwcart_topo_t), value :: hwtopo
       integer, value :: mpi_comm
       integer, value :: nlevels
       type(c_ptr), value :: domain
       type(c_ptr), value :: topo
       integer, value :: cart_order
       integer :: hwcart_print_rank_topology_f
     end function hwcart_print_rank_topology_f
     
     function hwcart_rank2coord_f(hwcart_comm, dims, rank, cart_order, coord_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       type(c_ptr), value :: dims
       integer, value :: rank
       integer, value :: cart_order
       type(c_ptr), value :: coord_out
     end function hwcart_rank2coord_f
     
     function hwcart_coord2rank_f(hwcart_comm, dims, periodic, coord, cart_order, rank_out) bind(C)
       use iso_c_binding
       integer, value :: hwcart_comm
       type(c_ptr), value :: dims
       type(c_ptr), value :: periodic
       type(c_ptr), value :: coord
       integer, value :: cart_order
       integer :: rank_out
     end function hwcart_coord2rank_f
     
  end interface
CONTAINS

  function hwcart_create(hwtopo, mpi_comm, domain, topo, cart_order, hwcart_comm_out)
    type(hwcart_topo_t), value :: hwtopo
    integer, value :: mpi_comm
    integer, dimension(:), target :: domain
    integer, dimension(:,:), target :: topo
    integer, value :: cart_order
    integer :: hwcart_comm_out
    integer :: hwcart_create
    integer :: nlevels

    nlevels = size(domain)
    hwcart_create = hwcart_create_f(hwtopo, mpi_comm, nlevels, c_loc(domain), c_loc(topo), cart_order, hwcart_comm_out)
  end function hwcart_create

  function hwcart_print_rank_topology(hwtopo, mpi_comm, domain, topo, cart_order)
    type(hwcart_topo_t), value :: hwtopo
    integer, value :: mpi_comm
    integer, dimension(:), target :: domain
    integer, dimension(:,:), target :: topo
    integer, value :: cart_order
    integer :: hwcart_comm_out
    integer :: hwcart_print_rank_topology
    integer :: nlevels
    
    nlevels = size(domain)
    hwcart_print_rank_topology = hwcart_print_rank_topology_f(hwtopo, mpi_comm, nlevels, c_loc(domain), c_loc(topo), cart_order)
  end function hwcart_print_rank_topology

  function hwcart_rank2coord(hwcart_comm, dims, rank, cart_order, coord_out)
    use iso_c_binding
    integer, value :: hwcart_comm
    integer, dimension(:), target :: dims(3)
    integer, value :: rank
    integer, value :: cart_order
    integer, dimension(:), target :: coord_out(3)
    integer :: hwcart_rank2coord

    hwcart_rank2coord = hwcart_rank2coord_f(hwcart_comm, c_loc(dims), rank, cart_order, c_loc(coord_out))
  end function hwcart_rank2coord
 
  function hwcart_coord2rank(hwcart_comm, dims, periodic, coord, cart_order, rank_out)
    use iso_c_binding
    integer, value :: hwcart_comm
    integer, dimension(:), target :: dims(3)
    integer, dimension(:), target :: periodic(3)
    integer, dimension(:), target :: coord(3)
    integer, value :: cart_order
    integer :: rank_out
    integer :: hwcart_coord2rank

    hwcart_coord2rank = hwcart_coord2rank_f(hwcart_comm, c_loc(dims), c_loc(periodic), c_loc(coord), cart_order, rank_out)
  end function hwcart_coord2rank
 
END MODULE hwcart_mod
