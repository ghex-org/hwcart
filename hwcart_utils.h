#ifndef _HWCART_UTILS_H

#include "hwcart.h"

#define HWCART_MPI_CALL(a) {                    \
        int ierr = a;                           \
        if(ierr) return ierr;                   \
    }

int  hwcart_split_type(int split_type);
void hwcart_split_type_to_name(int split_type, char *name);
int  hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, int order);
int  hwcart_topology(hwcart_topo_t hwtopo, MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank_out, int level);
int  hwcart_get_noderank(hwcart_topo_t hwtopo, MPI_Comm comm, int split_type, int *noderank_out);
int  hwcart_free_hwtopo(hwcart_topo_t *hwtopo);

#endif /* _HWCART_UTILS_H */
