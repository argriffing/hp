// Conventions:
// Component label -1 means unlabeled.
// Parent -2 means root.
// Parent -1 means not visited by depth first search.
//
// Connected components in csr format.


// The idea of this structure is to avoid malloc/free in an inner loop.
// This structure should not own the memory referenced by any of its pointers.
typedef struct tagSUBGRAPH {

  // This is the number of local vertices.
  int nvertices;

  // Points somewhere in ccgraph.compo_row_ptr.
  int *row_ptr;

  // Points somewhere in ccgraph.local_to_global.
  // Maps local vertex indices to global vertex indices.
  int *local_to_global;

} SUBGRAPH;


typedef struct tagCCGRAPH {
  int nvertices;
  int ncomponents;
  int nedges;
  LOCALGRAPH *subgraph; // length is nvertices (as upper bound on ncomponents)
  int *global_to_local; // length is nvertices
  int *local_to_global; // length is nvertices
  int *compo_row_ptr; // length is nvertices + 1
  int *compo_ind_ptr; // length is nedges
} CCGRAPH;


// This function needs to be called only once per application,
// if the max number of vertices and max number of edges is known.
// You do not need to call this function every time you want
// to compute the connected components.
void ccgraph_init(CCGRAPH *p, int max_nvertices, int max_nedges)
{
  p->subgraph = (LOCALGRAPH *) malloc(max_nvertices * sizeof(CCGRAPH));
  p->global_to_local = (int *) malloc(max_nvertices * sizeof(int));
  p->local_to_global = (int *) malloc(max_nvertices * sizeof(int));
  p->compo_row_ptr = (int *) malloc((max_nvertices + 1) * sizeof(int));
  p->compo_ind_ptr = (int *) malloc(max_nedges * sizeof(int));
}


// This function assumes that the ccgraph has already been initialized.
void ccgraph_compute(CCGRAPH *p,
    const int *row_ptr, const int *col_ind, int nvertices)
{
  p->ncomponents = 0;
  p->nvertices = 0;
  p->nedges = 0;
  p->compo_row_ptr[p->nvertices] = 0;
}


// This function is not specific to ccgraph.
void compute_component(
    const int *row_ptr, const int *col_ind,
    int *parent_ws, int *deck_ws, int *next_ws,
    int root, int label, int *component_labels,
    int *vertices_out, int *nvertices_out
    )
{
  // The component at the root should be unknown (-1).
  // Set it to the current label.
  assert(component_labels[root] == -1);
  component_labels[root] = label;

  // The root parent should be unknown (-1).
  // Set it to -2 to indicate that it is the root.
  assert(parent_ws[root] == -1);
  parent_ws[root] = -2;

  // The root vertex is on deck for breadth first search.
  deck_ws[0] = root;
  int ndeck = 1;

  // Do the breadth first search.
  // During the search update the local/global vertex index maps.
  // After the search update the local csr graph.
  while (ndeck)
  {
    int nnext = 0;
    int ideck;
    for (ideck=0; ideck<ndeck; ++ideck) {
      int v = deck_ws[ideck];

      // Update the map from the local to the global vertex index.
      // Increment the vertex count to include the added vertex.
      vertices_out[(*nvertices_out)++] = v;

      for (i=row_ptr[v]; i<row_ptr[v+1]; ++i) {
        int w = col_ind[i];
        if (w != parent_ws[v]) {
          if (parent_ws[w] >= 0) {
            assert(component_labels[w] == label);
          } else {
            assert(component_labels[w] == -1);
            component_labels[w] = label;
            parent_ws[w] = v;
            next_ws[nnext++] = w;
          }
        }
      }
    }

    // Put the next array on deck.
    int *tmp = deck_ws;
    deck_ws = next_ws;
    next_ws = tmp;
    ndeck = nnext;
  }
}


// This is a helper function for flood filling a single component.
// The function name underscore prefix indicates that this function is not
// part of the public interface of this module.
// C may have keywords to enforce this somehow but I'm not doing that.
//
// The parent vertices are assumed to have been set to -1,
// and all vertices in the connected component of the root are assumed to be -1
// representing membership in a currently unknown component.
// The set of vertices in the component is also tracked.
// The local to global map and the global to local map should each
// have -1 entries to represent unknown values.
//
void _ccgraph_compute_component(CCGRAPH *p,
    const int *row_ptr, const int *col_ind,
    int *parent_ws, int *deck_ws, int *next_ws,
    int root, int *component_labels
    )
{
  // The current label is the number of components in the ccgraph.
  int label = p->ncomponents;

  // Pick the subgraph that corresponds to this component.
  SUBGRAPH *subgraph = p->subgraph[p->ncomponents];

  // Set up the subgraph indirection relative to the ccgraph object.
  subgraph->local_to_global = p->local_to_global + p->nvertices;
  subgraph->row_ptr = p->compo_row_ptr + p->nvertices;

  // Initialize the subgraph to have no vertices.
  subgraph->nvertices = 0;

  // Compute the connected component.
  compute_component(row_ptr, col_ind,
      parent_ws, deck_ws, next_ws,
      root, label, component_labels,
      subgraph->local_to_global, &subgraph->nvertices)

  // Do some simple post-processing for each vertex in the component.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {

    // Get the global index of the vertex.
    v_global = subgraph->local_to_global[v_local];

    // Reset the vertex parent.
    parent_ws[v_global] = -1;

    // Fill the global to local map.
    p->global_to_local[v_global] = v_local;
  }

  // Update the subgraph csr.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {

    // Get the global index of the vertex.
    v_global = subgraph->local_to_global[v_local];

    // For each edge, add the local sink vertex to the compo col ind array.
    for (i=row_ptr[v_global]; i<row_ptr[v_global+1]; ++i) {
      int w_global = col_ind[i];
      int w_local = p->global_to_local[w_global];
      p->compo_col_ind[p->nedges++] = w_local;
    }

    // Add to the compo row ptr list.
    p->compo_row_ptr[p->nvertices + 1] = p->nedges;

    // Increment the number of vertices for which the edges are tracked.
    p->nvertices++;
  }

  // Increment the number of components for which the local graph is available.
  p->ncomponents++;
}


void ccgraph_destroy(CCGRAPH *p)
{
  free(p->subgraph);
  free(p->global_to_local);
  free(p->local_to_global);
  free(p->compo_row_ptr);
  free(p->compo_ind_ptr);
}





// Get the connected component associated with a single vertex.
// This is a helper for the function that gets all connected
// components of a graph.
// The parent vertices are assumed to have been set to -1,
// and all vertices in the connected component of the root are assumed to be -1
// representing membership in a currently unknown component.
// The set of vertices in the component is also tracked.
void fill_connected_component(
    const int *row_ptr, const int *col_ind,
    int *parent_ws, int *deck_ws, int *next_ws,
    int root_in, int label_in, int *component_labels_out,
    int *vertices_out, int *nvertices_out
    )
{

  // The component label array should have entries that are -1
  // when the label is unknown, and it has entries that are non-negative
  // if the labels have been determined.
  assert(component_labels_out[root_in] == -1);
  assert(parent_ws[r] == -2);

  // Set the root label and parent vertex,
  // and add the vertex to the list of vertices in the connected component.
  parent_ws[r] = -2;
  component_labels_out[root_in] = label_in;
  vertices_out[(*nvertices_out)++] = root_in;

  // The root vertex is on deck for breadth first search.
  deck_ws[0] = r;
  int ndeck = 1;

  // Do the breadth first search.
  int depth = 0;
  while (ndeck)
  {
    int nnext = 0;
    int ideck;
    for (ideck=0; ideck<ndeck; ++ideck) {
      int v = deck_ws[ideck];
      for (i=row_ptr[v]; i<row_ptr[v+1]; ++i) {
        int w = col_ind[i];
        if (w != parent_ws[v]) {

          // If we have already reached this node earlier
          // (possible if the graph has cycles) then skip the vertex.
          // Otherwise if we have not reached the vertex earlier
          // then set the label and the parent and add it to the shell,
          // and also add it to the list of vertices in the component,
          // and increment the number of vertices in the component.
          if (parent_ws[w] >= 0) {
            assert(component_labels_out[w] == label_in);
          } else {
            assert(component_labels_out[w] == -1);
            component_labels_out[w] = label_in;
            vertices_out[(*nvertices_out)++] = w;
            parent_ws[w] = v;
            next_ws[nnext++] = w;
          }
        }
      }
    }

    // Put the next array on deck.
    int *tmp = deck_ws;
    deck_ws = next_ws;
    next_ws = tmp;
    ndeck = nnext;
  }
}


// Get the connected components of a graph in csr format.
// The first argument group defines the graph.
// The second argument group provides workspace for the depth first search.
// The final argument groups provide actual outputs of interest to the caller.
void fill_connected_components(
    const int *row_ptr, const int *col_ind, int nvertices,
    int *parent_ws, int *deck_ws, int *next_ws,
    int *component_labels_out, int *ncomponents_out,
    int *compo_row_ptr_out, int *compo_col_ind_out
    )
{
  int root, v;
  *ncomponents_out = 0;

  // Set all vertex labels and all parent vertices to -1.
  for (v=0; v<nvertices; ++v) {
    component_labels_out[v] = -1;
    parent_ws[v] = -1;
  }

  // For each vertex that has not been assigned a component,
  // flood fill a component rooted at that vertex.
  for (root=0; root<nvertices; ++root) {
    if (component_labels_out[root] < 0) {

      // The component label is the current number components.
      label = *ncomponents_out;
      fill_connected_component(
          row_ptr, col_ind,
          parent_ws, deck_ws, next_ws,
          root, label, component_labels_out,
          component_ws, &nvertices)

      // For each vertex in the newly defined component,
      // set the parent vertex to -1.
      for (v=0; v<

      // Increment the number of components.
      ++(*ncomponents_out);
    }
  }
}

