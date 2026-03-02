#define LG_FREE_WORK          \
  {                           \
    GrB_free(&bottom_scalar); \
    GrB_free(&IPI_set);       \
    GrB_free(&IPI_mult);      \
    GrB_free(&PI_set);        \
    GrB_free(&PI_mult);       \
    GrB_free(&PI_add);        \
    GrB_free(&PI_semiring);   \
    GrB_free(&Theta);         \
    GrB_free(&PI_monoid);     \
  }

#include "LG_internal.h"
#include <LAGraphX.h>

void add_path_index(PathIndex *z, const PathIndex *x, const PathIndex *y)
{
  // If one is ⊥, then we take another.
  if (x->height == 0 && x->middle == 0)
  {
    *z = *y;
  }
  else if (y->height == 0 && y->middle == 0)
  {
    *z = *x;
  }
  // Take the path with the minimum height
  else if (x->height < y->height || (x->height == y->height && x->middle <= y->middle))
  {
    *z = *x;
  }
  else
  {
    *z = *y;
  }
}

void mult_path_index(PathIndex *z,
                     const PathIndex *x, GrB_Index ix, GrB_Index jx,
                     const PathIndex
                         *y,
                     GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
  // Сoncatenation with ⊥ gives ⊥
  if ((x->height == 0 && x->middle == 0) || (y->height == 0 && y->middle == 0))
  {
    z->height = 0;
    z->middle = 0;
  }
  else
  {
    z->middle = jx;
    z->height = (x->height > y->height ? x->height : y->height) + 1;
  }
}

void set_path_index(PathIndex *z,
                    const PathIndex *x, GrB_Index ix, GrB_Index jx,
                    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,
                    const void *theta)
{
  if (*edge_exist)
  {
    z->middle = ix;
    z->height = 1;
  }
  else
  {
    z->middle = 0;
    z->height = 0;
  }
}

#define SET_PATH_INDEX_DEFN                                                   \
  "void set_path_index(                                                   \n" \
  "    PathIndex *z,                                                      \n" \
  "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n" \
  "    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,        \n" \
  "    const void *theta)                                                 \n" \
  "{                                                                      \n" \
  "    if (*edge_exist)                                                   \n" \
  "    {                                                                  \n" \
  "        z->middle = ix;                                                \n" \
  "        z->height = 1;                                                 \n" \
  "    }                                                                  \n" \
  "    else                                                               \n" \
  "    {                                                                  \n" \
  "        z->middle = 0;                                                 \n" \
  "        z->height = 0;                                                 \n" \
  "    }                                                                  \n" \
  "}"

#define MULT_PATH_INDEX_DEFN                                                   \
  "void mult_path_index(                                                  \n"  \
  "    PathIndex *z,                                                      \n"  \
  "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n"  \
  "    const PathIndex *y, GrB_Index iy, GrB_Index jy,                    \n"  \
  "    const void *theta)                                                 \n"  \
  "{                                                                      \n"  \
  "    if ((x->height == 0 && x->middle == 0) ||                           \n" \
  "        (y->height == 0 && y->middle == 0))                             \n" \
  "    {                                                                  \n"  \
  "        z->height = 0;                                                 \n"  \
  "        z->middle = 0;                                                 \n"  \
  "    }                                                                  \n"  \
  "    else                                                               \n"  \
  "    {                                                                  \n"  \
  "        z->middle = jx;                                                \n"  \
  "        z->height = (x->height > y->height ?                           \n"  \
  "                     x->height : y->height) + 1;                       \n"  \
  "    }                                                                  \n"  \
  "}"

GrB_Info LAGraph_CFL_single_path(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) contains a PathIndex structure if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the non-terminal 'k' of the specified CFG.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, terms_count - 1]).
    int64_t terms_count,            // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
)
{
  // Semiring components
  GrB_Type PI_type = NULL; // Type PathIndex
  GrB_BinaryOp PI_add = NULL;
  GrB_Monoid PI_monoid = NULL;
  GxB_IndexBinaryOp IPI_mult = NULL;
  GrB_BinaryOp PI_mult = NULL;
  GrB_Semiring PI_semiring = NULL;
  GxB_IndexBinaryOp IPI_set = NULL;
  GrB_BinaryOp PI_set = NULL;
  GrB_Scalar Theta = NULL;
  GrB_Scalar bottom_scalar = NULL;

  GRB_TRY(GrB_Type_new(&PI_type, sizeof(PathIndex))); // the memory is not being freed yet

  // Theta cannot be NULL
  GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, false));

  PathIndex bottom = {0, 0};
  GRB_TRY(GrB_Scalar_new(&bottom_scalar, PI_type));
  GRB_TRY(GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom)));

  // Create semiring
  GRB_TRY(GrB_BinaryOp_new(
      &PI_add,
      (void *)add_path_index,
      PI_type,
      PI_type,
      PI_type));

  GRB_TRY(GrB_Monoid_new(
      &PI_monoid,
      PI_add,
      (void *)(&bottom))); // ⊥ - neutral element for the addition operation

  GRB_TRY(GxB_IndexBinaryOp_new(
      &IPI_mult,
      (void *)mult_path_index,
      PI_type,
      PI_type,
      PI_type,
      GrB_BOOL,
      "mult_path_index",
      MULT_PATH_INDEX_DEFN));

  GRB_TRY(GxB_BinaryOp_new_IndexOp(
      &PI_mult,
      IPI_mult,
      Theta));

  GRB_TRY(GrB_Semiring_new(
      &PI_semiring,
      PI_monoid,
      PI_mult));

  GRB_TRY(GxB_IndexBinaryOp_new(
      &IPI_set,
      (void *)set_path_index,
      PI_type,
      PI_type,
      GrB_BOOL,
      GrB_BOOL,
      "set_path_index",
      SET_PATH_INDEX_DEFN));

  GRB_TRY(GxB_BinaryOp_new_IndexOp(
      &PI_set,
      IPI_set,
      Theta));

  CFL_Semiring semiring = {.type = PI_type,
                           .semiring = PI_semiring,
                           .add = PI_add,
                           .mult = PI_mult,
                           .init_path = PI_set,
                           .bottom_scalar = bottom_scalar};
  LG_TRY(LAGraph_CFPQ_core(outputs, adj_matrices, terms_count, nonterms_count, rules, rules_count, &semiring, msg));
  LG_FREE_WORK;
  return GrB_SUCCESS;
}
