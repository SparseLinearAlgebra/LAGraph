#define LG_FREE_ALL                                  \
    {                                                \
        LAGraph_Free((void **)&output->paths, NULL); \
        output->count = 0;                           \
        output->capacity = 0;                        \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

#define ADD_TO_MSG(...)                                                   \
    {                                                                     \
        if (msg_len == 0)                                                 \
        {                                                                 \
            msg_len +=                                                    \
                snprintf(msg, LAGRAPH_MSG_LEN,                            \
                         "LAGraph failure (file %s, line %d): ",          \
                         __FILE__, __LINE__);                             \
        }                                                                 \
        if (msg_len < LAGRAPH_MSG_LEN)                                    \
        {                                                                 \
            msg_len += snprintf(msg + msg_len, LAGRAPH_MSG_LEN - msg_len, \
                                __VA_ARGS__);                             \
        }                                                                 \
    }

// LAGraph_CFL_extract_single_path: Context-Free Language Single Path Extraction Algorithm
//
// For a given pair of vertices (start, end) and a non-terminal symbol (nonterm),
// this function extracts a single path corresponding to a string with the minimum derivation tree height.
//
// The vertex parameters are optional:
//   - Passing NULL for start extracts paths from all possible start vertices
//   - Passing NULL for end extracts paths to all possible end vertices
//   - Passing NULL for both extracts one path for every pair of vertices in the graph
// In all cases, results are returned as a PathArray containing one Path per pair.
//
// Terminals and non-terminals are enumerated by integers starting from zero.
// The start non-terminal is the non-terminal with index 0.
//
// Note: This function must be called after LAGraph_CFL_single_path. The output
// from that function (path_index_matrices) is used as the basis for path extraction.
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
// Grammar: S -> aSb | ab | eps
//
// For non-terminal S (index 0) with different start/end parameters, the function extracts the following paths:
//
// 1. Single path extraction (start = 0, end = 4):
//    Returns path: 0 -a-> 1 -a-> 2 -b-> 3 -b-> 4 (word "aabb")
//
// 2. Epsilon path (start = end):
//    For any vertex v, calling with start = v, end = v returns an empty path (word eps)
//
// 3. No path exists (start = 0, end = 5):
//    A word formed in this way does not belong to the language - returns empty PathArray
//
// 4. Multiple paths (start = NULL, end = 3):
//    Returns PathArray containing all paths to vertex 3:
//    - (0, 3): 0 -a-> 1 -a-> 5 -b-> 2 -b-> 3 (word "aabb")
//    - (1, 3): 1 -a-> 2 -b-> 3 (word "ab")
//    - (3, 3): empty path (word eps)

GrB_Info LAGraph_CFL_extract_single_path(
    // Output
    PathArray *output, // Array of extracted paths.
                       // When both start and end are fixed (non-NULL),
                       // the array contains at most one path.
                       // When start or end is NULL (or both), the array may contain multiple paths — one for each valid pair.
    // Input
    GrB_Index *start,                      // Source vertex of a graph path.
                                           // Pass NULL to get paths from all vertices of the graph.
    GrB_Index *end,                        // Destination vertex of a graph path.
                                           // Pass NULL to get paths to all vertices of the graph.
    int32_t nonterm,                       // Non-terminal symbol to derive paths for.
    const GrB_Matrix *adj_matrices,        // Array of adjacency matrices representing the graph.
                                           // The length of this array is equal to the count of
                                           // terminals (terms_count).
                                           //
                                           // adj_matrices[t]: (i, j) == 1 if and only if there
                                           // is an edge between nodes i and j with the label of
                                           // the terminal corresponding to index 't' (where t is
                                           // in the range [0, terms_count - 1]).
    const GrB_Matrix *path_index_matrices, // Matrices containing information about existing paths for each non-terminal
                                           //
                                           // outputs[k]: (i, j) contains a PathIndex structure if and only if there is a path
                                           // from node i to node j whose edge labels form a word
                                           // derivable from the non-terminal 'k' of the specified CFG.
    int64_t terms_count,                   // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,                // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules,        // The rules of the CFG.
    int64_t rules_count,                   // The total number of rules in the CFG.
    char *msg                              // Message string for error reporting.
)
{
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    GrB_Info info;
    output->count = 0;
    output->paths = NULL;
    // Initial capacity 1, since most often looking for exactly 1 path with a fixed start and end
    output->capacity = 1;

    LG_CFL_CHECK_BASE_INPUTS(adj_matrices, terms_count, nonterms_count, rules_count, rules);
    LG_ASSERT_MSG(nonterm < nonterms_count, GrB_INVALID_VALUE,
                  "The start non-terminal must be no greater than the number of non-terminals.");
    LG_ASSERT_MSG(path_index_matrices != NULL, GrB_NULL_POINTER, "The path_index_matrices array cannot be null.");

    // Find null path_index_matrices matrices
    bool found_null = false;
    for (int64_t i = 0; i < nonterms_count; i++)
    {
        if (path_index_matrices[i] != NULL)
            continue;

        if (!found_null)
        {
            ADD_TO_MSG("path_index_matrices matrices with these indexes are null:");
        }
        ADD_TO_MSG(" %" PRId64, i);

        found_null = true;
    }
    if (found_null)
    {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, adj_matrices[0]));

    GrB_Index start_begin = (start == NULL) ? 0 : *start;
    GrB_Index start_end = (start == NULL) ? n - 1 : *start;
    GrB_Index end_begin = (end == NULL) ? 0 : *end;
    GrB_Index end_end = (end == NULL) ? n - 1 : *end;

    LG_TRY(LAGraph_Calloc((void **)&output->paths, rules_count, sizeof(size_t), msg));

    for (GrB_Index st = start_begin; st <= start_end; st++)
    {
        for (GrB_Index en = end_begin; en <= end_end; en++)
        {
            Path path;

            // Function that extracts one path with fixed start and end
            info = LAGraph_CFL_extract_single_path_internal(&path, st, en, nonterm, adj_matrices, path_index_matrices, terms_count, nonterms_count, rules, rules_count, msg);
            if (info == GrB_SUCCESS)
            {
                if (output->count == output->capacity)
                {
                    LG_TRY(LAGraph_Realloc((void **)&output->paths, output->capacity * 2, output->capacity,
                                           sizeof(Path), msg));
                    output->capacity *= 2;
                }
                output->paths[output->count++] = path;
            }
            else if (info != GrB_NO_VALUE) // GrB_NO_VALUE is the absence of a path, not an error
            {
                LG_FREE_ALL;
                return info;
            }
        }
    }

    // No paths found - free previously allocated memory
    if (output->count == 0)
    {
        LG_FREE_ALL;
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
