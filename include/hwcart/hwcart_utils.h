#ifndef _HWCART_UTILS_H

#include "hwcart.h"

#define HWCART_MPI_CALL(a) {                    \
        int ierr = a;                           \
        if(ierr) return ierr;                   \
    }

typedef struct hwcart_topo_struct_t* hwcart_topo_t;

void hwcart_split_type_to_name(int split_type, char *name);
int  hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, hwcart_order_t order);
int  hwcart_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nlevels, hwcart_split_t *domain, int *topo, int *level_rank_out, int level);
int  hwcart_get_noderank(hwcart_topo_t hwtopo, MPI_Comm comm, hwcart_split_t hwcart_split_type, int *noderank_out);

#endif /* _HWCART_UTILS_H */
