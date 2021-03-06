#ifndef SPARSETOOLS_HEADER
#define SPARSETOOLS_HEADER

// Sparse matrix format conversion.
// List of lists to compressed row storage.
void lil_to_csr(int nvertices, int nedges,
    int *lil_in, int *row_ptr_out, int *col_ind_out);

void csr_degree_min_max(
    const int *row_ptr, int nvertices,
    int *min_degree, int *max_degree);

int csr_count_loops(const int *row_ptr, const int *col_ind, int n);

void print_csr_graph(const int *row_ptr, const int *col_ind, int n);

#endif
