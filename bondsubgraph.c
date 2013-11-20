// Functions downstream of the 2d hp bond graph.
//
// The main function in this module computes the densest k-subgraph
// of a 2d hp bond graph.
// This can be done efficiently because of the particular
// constraints on the graph (all of its vertices have degree less than three,
// except for at most two vertices which may be of degree three).
//
// Although this function is not of general interest,
// it depends on functions of more general interest,
// such as segmentation into connected components,
// calculation of graph girth, and the solutions of subset sum problems.
//
// The bsg_ function name prefix stands for "bond subgraph".


#include "assert.h"
#include "stdlib.h"

#include "connectedcomponents.h"
#include "subsetsum.h"
#include "graphgirth.h"
#include "breadthfirst.h"
#include "bondsubgraph.h"

// This is only for allocating memory.
void solver_init(SOLVER *solver, int max_nvertices, int max_ncomponents)
{
  solver->k = -1;
  solver->nsolution = -1;
  solver->score = -1;
  solver->ncomponents = -1;
  solver->solution = (int *) malloc(max_nvertices * sizeof(int));
  solver->data = (BSG_COMPONENT *) malloc(
      max_ncomponents * sizeof(BSG_COMPONENT));
  solver->components = (BSG_COMPONENT **) malloc(
      max_ncomponents * sizeof(BSG_COMPONENT **));
}

// This is only for freeing memory.
void solver_destroy(SOLVER *solver)
{
  free(solver->solution);
  free(solver->data);
  free(solver->components);
}

void solver_swap_components(SOLVER *solver, int ca, int cb)
{
  BSG_COMPONENT *tmp = solver->components[ca];
  solver->components[ca] = solver->components[cb];
  solver->components[cb] = tmp;
}

// Swap the component pointer at the given index
// with the component pointer at the last index.
// Decrement the number of components.
// Note that the component does not index into the ccgraph
// array of connected components, but rather into the solver array
// of pointers to connected components.
void solver_remove_component(SOLVER *solver, int component)
{
  solver_swap_components(solver, component, solver->ncomponents - 1);
  solver->ncomponents--;
}

// Call this before each solving attempt.
void solver_prepare(SOLVER *solver,
    const int *row_ptr, const int *col_ind,
    CCGRAPH *ccgraph, int k,
    int *va_trace, int *vb_trace,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  solver->k = k;
  solver->nsolution = 0;
  solver->score = 0;
  solver->ncomponents = ccgraph->ncomponents;

  // For each component, move the smallest cycle to the front of the array
  // and order all vertices so that each new vertex is adjacent
  // to at least one earlier vertex.
  // Also set some component attributes.
  int c;
  BSG_COMPONENT *bsg;
  for (c=0; c<solver->ncomponents; ++c) {
    int computed_girth = _move_smallest_cycle_to_front(
        row_ptr, col_ind,
        ccgraph, c,
        va_trace, vb_trace,
        bfs_ws, depth_ws);
    bsg = solver->components[c];
    bsg->nvertices = ccgraph_get_component_nvertices(ccgraph, c);
    bsg->nedges = ccgraph_get_component_nedges(ccgraph, c);
    bsg->ell = bsg->nedges - bsg->nvertices;
    bsg->girth = computed_girth;
    bsg->component = c;
  }
}

// Return -1 if no special component is found.
int solver_get_special_component_index(SOLVER *solver) {
  int c;
  for (c=0; c<solver->ncomponents; ++c) {
    if (solver->components[c]->ell > 1) {
      return c;
    }
  }
  return -1;
}

// Process the special component if it exists.
void solver_do_special(SOLVER *solver,
    const int *row_ptr, const int *col_ind,
    CCGRAPH *ccgraph
    )
{
  // Get the special component index.
  // If none was found, then return.
  int special_idx = solver_get_special_component_index(solver);
  if (special_idx < 0) {
    return;
  }

  // If the number of vertices in the special component is at most k,
  // then add the special component into the solution
  // and remove it from the list of available components.
  
  // Add the vertices of the special component into the solution array.
  int v_local, v_global;
  BSG_COMPONENT *bsg = solver->components[special_idx];
  for (v_local=0; v_local<bsg->nvertices; ++v_local) {
    v_global = ccgraph_local_to_global(ccgraph, bsg->component, v_local);
    solver->solution[solver->nsolution++] = v_global;
  }

  // Add the number of edges to the score.
  solver->score += bsg->nedges;

  // Swap the special component with the component at the back
  // of the component array, and decrement the number of components.
  solver_remove_component(solver, special_idx);
}

// Helper function for _qsort_components().
// If the return value is negative
// then the first element goes before the second element.
int _cmp_components(const void *a, const void *b)
{
  BSG_COMPONENT *pa = *(BSG_COMPONENT **) a;
  BSG_COMPONENT *pb = *(BSG_COMPONENT **) b;

  // If one of the components has a cycle
  // but the other does not, then the one with the cycle goes in front.
  // If both components have a cycle then order does not matter.
  // If neither component has a cycle then we defer to another criterion.
  int pa_has_cycle = (pa->girth > -1);
  int pb_has_cycle = (pb->girth > -1);
  if (pa_has_cycle && !pb_has_cycle) {
    return -1;
  } else if (!pa_has_cycle && pb_has_cycle) {
    return 1;
  } else if (pa_has_cycle && pb_has_cycle) {
    return 0;
  }

  // If neither component has a cycle,
  // then sort by decreasing number of vertices.
  return pb->nvertices - pa->nvertices;
}

// Sort an array of component pointers in decreasing order.
void _qsort_components(BSG_COMPONENT **components, int n)
{
  qsort(components, n, sizeof(BSG_COMPONENT *), _cmp_components);
}

// Sort all of the remaining components.
// This is not efficient, because it could be broken into three steps.
// The cycle-containing components could be pushed to the front,
// the isolated vertices could be pushed to the back,
// and the remaining tree-like shapes in the middle could be sorted
// by decreasing length.
// But instead, I will write my code inefficiently
// and just sort the whole thing.
void solver_sort_components(SOLVER *solver)
{
  _qsort_components(solver->components, solver->ncomponents);
}


// The special component has been dealt with.
// All of the usable components that have a cycle have been moved
// to the front of the list of components.
// If the sum of their vertices is at least as large as k,
// then we will try to solve a subset sum problem
// to get k edges using k vertices.
// If k is greater than the sum of the vertices,
// then the subset sum problem has no solution so we skip it.
//
// The s3table object and the low, high, and s3solution arrays are workspaces
// for solving the subset sum problem.
//
void solver_subset_sum(SOLVER *solver, CCGRAPH *ccgraph,
    S3TABLE *s3table, int *low, int *high, int *s3solution
    )
{
  // Get the number of available components that have cycles.
  // Get the sum of vertices of these components.
  int cyclic_component_count = 0;
  int acyclic_component_count = 0;
  int cyclic_component_vertex_count = 0;
  int c;
  for (c=0; c<solver->ncomponents; ++c)
  {
    if (solver->components[c]->girth > 0) {

      // Components with cycles are not allowed to appear
      // after components without cycles according to the sorting criterion.
      assert(!cyclic_component_count);

      // Take note of the cycle-containing component.
      cyclic_component_count++;
      cyclic_component_vertex_count += solver->components[c]->nvertices;

    } else {
      acyclic_component_count++;
    }
  }

  // If the cycle-containing components do not collectively have enough
  // vertices, then no solution to the subset sum problem is possible.
  if (solver->k > cyclic_component_vertex_count) {
    return;
  }

  // Construct the low array and the high array for the subset sum solver.
  // These values are related to the girth and the number of vertices
  // of the components that contain cycles.
  for (c=0; c<cyclic_component_count; ++c) {
    low[c] = solver->components[c]->girth;
    high[c] = solver->components[c]->nvertices;
  }

  // Attempt to solve the subset sum problem.
  // Clean up the subset sum solving workspace.
  int target = solver->k;
  int ncomponents = cyclic_component_count;
  s3table_forward(s3table, low, high, ncomponents, target);
  s3table_backward(s3table, ncomponents, target, s3solution);
  int success = s3table_attainable(s3table, ncomponents-1, target);
  s3table_clear(s3table, ncomponents, target);

  // If the attempt failed, then return.
  if (!success) {
    return;
  }

  // If we have successfully found a solution to the subset sum problem,
  // then update the vertex set for the densest k-subgraph solution.
  int v_local, v_global;
  int nvertices_added = 0;
  for (c=0; c<cyclic_component_count; ++c) {
    BSG_COMPONENT *bsg = solver->components[c];
    for (v_local=0; v_local<s3solution[c]; ++v_local) {
      v_global = ccgraph_local_to_global(ccgraph, bsg->component, v_local);
      solver->solution[solver->nsolution++] = v_global;
      nvertices_added++;
    }
  }

  // Check that the number of vertices added is exactly k.
  // Then set k to zero, showing that we have finished the search.
  assert(solver->k == nvertices_added);
  solver->k = 0;

  // Because of the graph constraints and the vertex ordering,
  // the number of edges added is exactly equal to the number
  // of vertices added, when the subset sum problem has been solved.
  solver->score += nvertices_added;
}


// Components are sorted with cycle-containing components appearing first,
// followed by tree-like components.
// The components are secondarily sorted according to decreasing
// number of vertices.
// Because we have failed the subset sum search by this point,
// there is no way to add k or more edges by adding k vertices.
// We are penalized by one point if we include a partial cycle-containing
// component -- note that we are assuming that we include at most
// one such component and that because we have failed the subset sum
// problem, this partial cycle-containing component cannot contain
// a whole cycle.
// We are also penalized one point
// for each whole and each partial cycle-less component.
void solver_greedy(SOLVER *solver, CCGRAPH *ccgraph)
{
  int v, v_global;
  int npenalties = 0;
  int nvertices_added = 0;
  int c;
  for (c=0; c<solver->ncomponents && nvertices_added < solver->k; ++c) {
    BSG_COMPONENT *bsg = solver->components[c];
    for (v=0; v<bsg->nvertices && nvertices_added < solver->k; ++v) {
      v_global = ccgraph_local_to_global(ccgraph, bsg->component, v);
      solver->solution[solver->nsolution++] = v_global;
      nvertices_added++;
      
      // We incur a penalty if we include a partial cycle-containing component.
      if (nvertices_added == solver->k && v < bsg->nvertices - 1) {
        npenalties++;
      }

      // We incur a penalty for each cycle-free component regardless
      // of whether or not the whole component is included.
      if (v == 0 && bsg->girth < 0) {
        npenalties++;
      }
    }
  }

  // Assert that we have incurred at least one penalty
  // and that we have added a number of vertices equal to k.
  assert(npenalties > 0);
  assert(nvertices_added == solver->k);

  // Add the score, which is k minus the penalties, and set k to zero.
  solver->score += solver->k - npenalties;
  solver->k = 0;
}


// Reorder vertices within a connected component of an undirected graph
// so that the vertices in the smallest cycle are before the other vertices,
// and so that each vertex is adjacent to a vertex somewhere
// before it in the vertex list.
// The length of the smallest cycle, also called the girth, is returned.
// If the component has no cycle, then return -1 to represent an infinite girth
// and do not change the vertex order.
int _move_smallest_cycle_to_front(
    const int *row_ptr, const int *col_ind,
    CCGRAPH *ccgraph, int component,
    int *va_trace, int *vb_trace,
    BFS_WS *bfs_ws, int *depth_ws
    )
{
  int nedges = ccgraph_get_component_nedges(ccgraph, component);
  int nvertices = ccgraph_get_component_nvertices(ccgraph, component);
  int ell = nedges - nvertices;
  int girth = -1;
  SUBGRAPH *subgraph = ccgraph->subgraph + component;

  // If the number of vertices is much greater than the number
  // of edges, then we do not have a connected component.
  assert(ell >= -1);

  // If the number of vertices is greater than the number of edges
  // then we have something like a point or path or tree.
  // In any case we have no cycles, so do not change the order,
  // and return -1 indicating infinite girth.
  if (ell == -1) {
    girth = -1;
    return girth;
  }

  // If the number of vertices is equal to the number of edges
  // then we have a pure cycle.
  // Because there are no other vertices in this component,
  // this cycle is already at the front of the vertex list,
  // so we do not have to reorder any vertices.
  // The girth is the number of vertices in this cycle.
  if (ell == 0) {
    girth = nvertices;
    return girth;
  }

  // For more interesting connected components,
  // each cycle includes at least one vertex of degree at least three.
  // So we can use a modified topo sort to check each of these vertices
  // to see which one is part of the smallest cycle in the component.
  int *local_row_ptr = ccgraph_get_component_row_ptr(ccgraph, component);
  int *local_col_ind = ccgraph_get_component_col_ind(ccgraph, component);
  int local_witness = -1;
  get_girth_and_vertex_conn(local_row_ptr, local_col_ind, nvertices,
      bfs_ws, depth_ws, &girth, &local_witness);
  int global_witness = subgraph->local_to_global[local_witness];

  // Put the entries of the smallest cycle into the local to global map.
  int ncycle = -1;
  get_smallest_cycle_ub(
      row_ptr, col_ind, global_witness,
      bfs_ws, depth_ws,
      va_trace, vb_trace,
      subgraph->local_to_global, &ncycle);
  assert(ncycle == girth);

  // Use breadth first search to fill the rest of the connected component.
  int nfilled = bfs_fill(row_ptr, col_ind,
      bfs_ws->parent, subgraph->local_to_global, ncycle);
  assert(nfilled == nvertices);
  bfs_clear(bfs_ws->parent, subgraph->local_to_global, nfilled);

  int i;
  int v_local, v_global;
  int w_local, w_global;

  // Translate this global vertex ordering back into local coordinates.
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {
    v_global = subgraph->local_to_global[v_local];
    ccgraph->global_to_local[v_global] = v_local;
  }

  // Update the subgraph csr.
  int nedges_updated = local_row_ptr[0];
  for (v_local=0; v_local<subgraph->nvertices; ++v_local) {

    // Get the global index of the vertex.
    v_global = subgraph->local_to_global[v_local];

    // For each edge, add the local sink vertex to the compo col ind array.
    for (i=row_ptr[v_global]; i<row_ptr[v_global+1]; ++i) {
      w_global = col_ind[i];
      w_local = ccgraph->global_to_local[w_global];
      ccgraph->compo_col_ind[nedges_updated++] = w_local;
    }

    // Add to the compo row ptr list.
    local_row_ptr[v_local + 1] = nedges_updated;
  }

  return girth;
}


// The connected component decomposition is assumed to have been performed.
int solver_toplevel(SOLVER *solver,
    const int *row_ptr, const int *col_ind, CCGRAPH *ccgraph, int k,
    int *va_trace, int *vb_trace,
    BFS_WS *bfs_ws, int *depth_ws,
    S3TABLE *s3table, int *low, int *high, int *s3solution
    )
{
  // Check that enough vertices are available.
  assert(solver->k <= ccgraph->nvertices);

  // Prepare the solver.
  // This reorders the vertices within the components that contain cycles.
  solver_prepare(solver,
      row_ptr, col_ind, ccgraph, k,
      va_trace, vb_trace,
      bfs_ws, depth_ws);

  // Deal with the special component if it exists.
  solver_do_special(solver, row_ptr, col_ind, ccgraph);
  if (!solver->k) return;

  // Components with cycles have priority;
  // secondarily, large components are preferred.
  solver_sort_components(solver);

  // Attempt to solve the subset sum problem, if appropriate.
  solver_subset_sum(solver, ccgraph, s3table, low, high, s3solution);
  if (!solver->k) return;

  // Add components according to their order within the array,
  // and add vertices according to their order with their component.
  solver_greedy(solver, ccgraph);
  assert(!solver->k);
}


