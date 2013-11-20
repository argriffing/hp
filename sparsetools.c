// Sparse matrix format conversion.
// List of lists to compressed row storage.
void lil_to_csr(int nvertices, int nedges,
    int *lil_in, int *row_ptr_out, int *col_ind_out)
{
  row_ptr_out[0] = 0;

  int v, e;
  int va, vb;
  int nfilled = 0;
  for (v=0; v<nvertices; ++v) {
    for (e=0; e<nedges; ++e) {
      int va = lil_in[e*2 + 0];
      int vb = lil_in[e*2 + 1];
      if (v == va) {
        col_ind_out[nfilled++] = vb;
      }
      if (v == vb) {
        col_ind_out[nfilled++] = va;
      }
    }
    row_ptr_out[v+1] = nfilled;
  }
}

void csr_degree_min_max(
    const int *row_ptr, int nvertices,
    int *min_degree, int *max_degree)
{
  *min_degree = -1;
  *max_degree = -1;
  int i, v;
  for (v=0; v<nvertices; ++v) {
    int degree = row_ptr[v+1] - row_ptr[v];
    if (*min_degree < 0 || degree < *min_degree) {
      *min_degree = degree;
    }
    if (*max_degree < 0 || degree > *max_degree) {
      *max_degree = degree;
    }
  }
}

