#ifndef _HWCART_UTILS_H

int  hwcart_split_type(int split_type);
void hwcart_split_type_to_name(int split_type, char *name);
void hwcart_print_cube(MPI_Comm comm, int *gdim, int id, int *buff, int line_size, int order);
int  hwcart_topology(MPI_Comm comm, int nsplits, int *domain, int *topo, int *level_rank_out, int level);
int  hwcart_get_noderank(MPI_Comm comm, int split_type, int *noderank_out);

#endif /* _HWCART_UTILS_H */
