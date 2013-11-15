// Conventions:
// Component label -1 means unlabeled.
// Parent -2 means root.
// Parent -1 means not visited by depth first search.
//
// Connected components in csr format.


typedef struct tagLOCALGRAPH {

  // This is the number of local vertices.
  int nvertices;

  // Points somewhere in ccgraph.compo_row_ptr.
  int *row_ptr;

  // Points somewhere in ccgraph.compo_ind_ptr.
  int *ind_ptr;

  // Points somewhere in ccgraph.local_to_global.
  // Maps local vertex indices to global vertex indices.
  int *local_to_global;

} LOCALGRAPH;


typedef struct tagCCGRAPH {
  int nvertices;
  int ncomponents;
  LOCALGRAPH *subgraphs;
  int *global_to_local; // length is nvertices
  int *local_to_global; // length is nvertices
  int *compo_row_ptr; // length is nvertices + 1
  int *compo_ind_ptr; // length is nedges
} CCGRAPH;

void ccgraph_init(CCGRAPH *p,
    const int *row_ptr, const int *col_ind, int nvertices)
{
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

