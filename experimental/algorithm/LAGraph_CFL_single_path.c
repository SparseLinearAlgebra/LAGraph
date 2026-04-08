#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        GrB_free(&bottom_scalar);                                                        \
        GrB_free(&IPathIndex_set);                                                       \
        GrB_free(&IPathIndex_mult);                                                      \
        GrB_free(&PathIndex_set);                                                        \
        GrB_free(&PathIndex_mult);                                                       \
        GrB_free(&PathIndex_add);                                                        \
        GrB_free(&PathIndex_semiring);                                                   \
        GrB_free(&Theta);                                                                \
        GrB_free(&PathIndex_monoid);                                                     \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

void add_path_index(PathIndex *z, const PathIndex *x, const PathIndex *y) {
    // If one is ⊥, then we take another.
    if (x->height == 0 && x->middle == 0) {
        *z = *y;
    } else if (y->height == 0 && y->middle == 0) {
        *z = *x;
    }
    // Take the path with the minimum height
    else if (x->height < y->height ||
             (x->height == y->height && x->middle <= y->middle)) {
        *z = *x;
    } else {
        *z = *y;
    }
}

void mult_path_index(PathIndex *z, const PathIndex *x, GrB_Index ix, GrB_Index jx,
                     const PathIndex *y, GrB_Index iy, GrB_Index jy, const void *theta) {
    // Сoncatenation with ⊥ gives ⊥
    if ((x->height == 0 && x->middle == 0) || (y->height == 0 && y->middle == 0)) {
        z->height = 0;
        z->middle = 0;
    } else {
        z->middle = jx;
        z->height = (x->height > y->height ? x->height : y->height) + 1;
    }
}

void set_path_index(PathIndex *z, const PathIndex *x, GrB_Index ix, GrB_Index jx,
                    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,
                    const void *theta) {
    if (*edge_exist) {
        z->middle = ix;
        z->height = 1;
    } else {
        z->middle = 0;
        z->height = 0;
    }
}

#define SET_PATH_INDEX_DEFN                                                              \
    "void set_path_index(                                                   \n"          \
    "    PathIndex *z,                                                      \n"          \
    "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n"          \
    "    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,        \n"          \
    "    const void *theta)                                                 \n"          \
    "{                                                                      \n"          \
    "    if (*edge_exist)                                                   \n"          \
    "    {                                                                  \n"          \
    "        z->middle = ix;                                                \n"          \
    "        z->height = 1;                                                 \n"          \
    "    }                                                                  \n"          \
    "    else                                                               \n"          \
    "    {                                                                  \n"          \
    "        z->middle = 0;                                                 \n"          \
    "        z->height = 0;                                                 \n"          \
    "    }                                                                  \n"          \
    "}"

#define MULT_PATH_INDEX_DEFN                                                             \
    "void mult_path_index(                                                  \n"          \
    "    PathIndex *z,                                                      \n"          \
    "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n"          \
    "    const PathIndex *y, GrB_Index iy, GrB_Index jy,                    \n"          \
    "    const void *theta)                                                 \n"          \
    "{                                                                      \n"          \
    "    if ((x->height == 0 && x->middle == 0) ||                           \n"         \
    "        (y->height == 0 && y->middle == 0))                             \n"         \
    "    {                                                                  \n"          \
    "        z->height = 0;                                                 \n"          \
    "        z->middle = 0;                                                 \n"          \
    "    }                                                                  \n"          \
    "    else                                                               \n"          \
    "    {                                                                  \n"          \
    "        z->middle = jx;                                                \n"          \
    "        z->height = (x->height > y->height ?                           \n"          \
    "                     x->height : y->height) + 1;                       \n"          \
    "    }                                                                  \n"          \
    "}"

// LAGraph_CFL_single_path: Context-Free Language Single Path Matrix-Based Algorithm
//
// This function computes, for each non-terminal and each pair of vertices (u, v) in the
// graph (represented by adjacency matrices), one path from u to v where the edge labels
// form a word from the language generated by the context-free grammar (represented by
// `rules`) with the minimum derivation tree height.
//
// Terminals and non-terminals are enumerated by integers starting from zero.
// The start non-terminal is the non-terminal with index 0.
//
// Note: This function does not compute the actual path, only information about it
// (represented by the PathIndex structure). To reconstruct the full path as PathArray and
// Path structures, use the LAGraph_CFL_extract_single_path function.
//
// Important: Do not free path_index_t until all work with the output matrices is
// complete. Accessing matrices after freeing their type is undefined behavior.
//
// Example:
//
// Graph:
// ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
// │ 0 ├───► 1 ├───► 2 ├───► 3 ├───► 4 │
// └───┘ a └─┬─┘ a └─▲─┘ b └───┘ b └───┘
//           │       │
//           │ ┌───┐ │
//          a└─► 5 ├─┘b
//             └───┘
//
// Grammar: S -> aSb | ab
//
// There are paths from node [1] to node [3] and from node [1] to node [2] that form the
// word "ab" ([1]-a->[2]-b->[3] and [1]-a->[5]-b->[2]). The word "ab" is in the language
// generated by our context-free grammar, so the pairs (1, 3) and (1, 2) will contain the
// corresponding path information.
//
// In contrast, the path from node [1] to node [4] forms the word "abb"
// ([1]-a->[2]-b->[3]-b->[4]) and the word "abbb" ([1]-a->[5]-b->[2]-b->[3]-b->[4]). The
// words "abb" and "abbb" are not in the language, so the pair (1, 4) will not contain any
// path information.
//
// With this graph and grammar, information about the found paths from the start
// non-terminal will be stored in the following matrix entries:
// (0, 4) - because there exists a path (0-1-2-3-4) that forms the word "aabb"
// (1, 3) - because there exists a path (1-2-3) that forms "ab"
// (1, 2) - because there exists a path (1-5-2) that forms the word "ab"
// (0, 3) - because there exists a path (0-1-5-2-3) that forms the word "aabb"

GrB_Info LAGraph_CFL_single_path(
    // Output
    GrB_Matrix
        *outputs, // Array of matrices containing results.
                  // The size of the array must be equal to nonterms_count.
                  //
                  // outputs[k]: (i, j) contains a PathIndex structure if and only if
                  // there is a path from node i to node j whose edge labels form a word
                  // derivable from the non-terminal 'k' of the specified CFG.
    GrB_Type *path_index_t, // PathIndex type - elements of the output matrices
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
    int64_t nonterms_count, // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
) {

#if GxB_IMPLEMENTATION < GxB_VERSION(9, 4, 5)
    return (GrB_NOT_IMPLEMENTED);
#else
    // Semiring components
    GrB_BinaryOp PathIndex_add = NULL;
    GrB_Monoid PathIndex_monoid = NULL;
    GxB_IndexBinaryOp IPathIndex_mult = NULL;
    GrB_BinaryOp PathIndex_mult = NULL;
    GrB_Semiring PathIndex_semiring = NULL;
    GxB_IndexBinaryOp IPathIndex_set = NULL;
    GrB_BinaryOp PathIndex_set = NULL;
    GrB_Scalar Theta = NULL;
    GrB_Scalar bottom_scalar = NULL;

    GRB_TRY(GrB_Type_new(path_index_t, sizeof(PathIndex)));
    GrB_Type PathIndex_type = *path_index_t;

    // Theta cannot be NULL
    GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, false));

    PathIndex bottom = {0, 0};
    GRB_TRY(GrB_Scalar_new(&bottom_scalar, PathIndex_type));
    GRB_TRY(GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom)));

    // Create semiring
    GRB_TRY(GrB_BinaryOp_new(&PathIndex_add, (void *)add_path_index, PathIndex_type,
                             PathIndex_type, PathIndex_type));

    GRB_TRY(GrB_Monoid_new(
        &PathIndex_monoid, PathIndex_add,
        (void *)(&bottom))); // ⊥ - neutral element for the addition operation

    GRB_TRY(GxB_IndexBinaryOp_new(&IPathIndex_mult, (void *)mult_path_index,
                                  PathIndex_type, PathIndex_type, PathIndex_type,
                                  GrB_BOOL, "mult_path_index", MULT_PATH_INDEX_DEFN));

    GRB_TRY(GxB_BinaryOp_new_IndexOp(&PathIndex_mult, IPathIndex_mult, Theta));

    GRB_TRY(GrB_Semiring_new(&PathIndex_semiring, PathIndex_monoid, PathIndex_mult));

    GRB_TRY(GxB_IndexBinaryOp_new(&IPathIndex_set, (void *)set_path_index, PathIndex_type,
                                  PathIndex_type, GrB_BOOL, GrB_BOOL, "set_path_index",
                                  SET_PATH_INDEX_DEFN));

    GRB_TRY(GxB_BinaryOp_new_IndexOp(&PathIndex_set, IPathIndex_set, Theta));

    CFL_Semiring semiring = {.type = PathIndex_type,
                             .semiring = PathIndex_semiring,
                             .add = PathIndex_add,
                             .mult = PathIndex_mult,
                             .init_path = PathIndex_set,
                             .bottom_scalar = bottom_scalar,
                             .get_nvals = GrB_Matrix_nvals};
    LG_TRY(LAGraph_CFPQ_core(outputs, adj_matrices, terms_count, nonterms_count, rules,
                             rules_count, &semiring, msg));
    LG_FREE_WORK;
    return GrB_SUCCESS;
#endif
}