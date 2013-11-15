#ifndef SPARSETOOLS_HEADER
#define SPARSETOOLS_HEADER

// Sparse matrix format conversion.
// List of lists to compressed row storage.
void lil_to_csr(int nvertices, int nedges,
    int *lil_in, int *row_ptr_out, int *col_ind_out);

#endif
