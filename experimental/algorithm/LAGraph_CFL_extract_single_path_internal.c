#define LG_FREE_WORK                              \
    {                                             \
        LAGraph_Free((void **)&eps_rules, NULL);  \
        LAGraph_Free((void **)&term_rules, NULL); \
        LAGraph_Free((void **)&bin_rules, NULL);  \
    }

#define LG_FREE_ALL                                  \
    {                                                \
        LAGraph_Free((void **)&output->edges, NULL); \
        output->len = 0;                             \
        LG_FREE_WORK;                                \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

// Internal helper for LAGraph_CFL_extract_single_path.
// Extracts a path for a fixed (start, end) vertex pair.
GrB_Info LAGraph_CFL_extract_single_path_internal(
    // Output
    Path *output, // A path extracted from a graph that forms a word produced by the rules of the grammar.
                  //
                  // If the path is not found, it is equal to the empty path: len = 0, path = NULL and GrB_NO_VALUE is returned.
                  // If the path is empty, i.e., formed by the rule: nonterminal -> eps, then the empty path: len = 0, path = NULL and GrB_SUCCESS is returned.
    // Input
    GrB_Index start,                       // Source vertex of a graph path.
    GrB_Index end,                         // Destination vertex of a graph path.
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
    output->len = 0;
    output->edges = NULL;
    // Arrays for processing rules
    size_t *eps_rules = NULL, eps_rules_count = 0;   // [Variable -> eps]
    size_t *term_rules = NULL, term_rules_count = 0; // [Variable -> term]
    size_t *bin_rules = NULL, bin_rules_count = 0;   // [Variable -> AB]

    LG_TRY(LAGraph_Calloc((void **)&eps_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&term_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&bin_rules, rules_count, sizeof(size_t), msg));

    LAGraph_CFL_classify_rules(eps_rules, &eps_rules_count, term_rules, &term_rules_count, bin_rules, &bin_rules_count, rules, rules_count);

    // Internal function: assumes all inputs have been validated by the public wrapper

    PathIndex index;
    GrB_Info info = GrB_Matrix_extractElement_UDT(&index, path_index_matrices[nonterm], start, end);
    if (info == GrB_SUCCESS) // Such a path exists
    {
        if (index.height == 1)
        {
            // Check for eps-path first (height=1, start=end, shortest possible)
            if (start == end)
            {
                for (size_t i = 0; i < eps_rules_count; i++)
                {
                    LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];
                    if (eps_rule.nonterm == nonterm)
                    {
                        LG_FREE_WORK;
                        return GrB_SUCCESS;
                    }
                }
            }
            // Height = 1 and no eps-path found - check for a term-path
            for (int64_t i = 0; i < terms_count; i++)
            {
                bool edge_exist;
                if (GrB_Matrix_extractElement_BOOL(&edge_exist, adj_matrices[i], start, end) == GrB_SUCCESS)
                {
                    for (size_t j = 0; j < term_rules_count; j++)
                    {
                        LAGraph_rule_WCNF term_rule = rules[term_rules[j]];
                        if (term_rule.nonterm == nonterm && term_rule.prod_A == i)
                        {
                            LG_TRY(LAGraph_Calloc((void **)&output->edges, 1, sizeof(Edge), msg));
                            output->len = 1;
                            output->edges[0] = (Edge){start, i, end};
                            LG_FREE_WORK;
                            return GrB_SUCCESS;
                        }
                    }
                }
            }
            // If couldn't find rules for outputting an empty or terminal path,
            // then the path were looking for doesn't match the rules
            LG_FREE_WORK;
            ADD_TO_MSG(msg_len, "The extracted path does not match the input grammar.");
            return GrB_NO_VALUE;
        }
        // Rules of the form Nonterm -> Nonterm * Nonterm are traversed recursively and merged
        for (size_t i = 0; i < bin_rules_count; i++)
        {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];
            if (bin_rule.nonterm != nonterm)
            {
                continue;
            }
            PathIndex indexB, indexC;
            if ((info = GrB_Matrix_extractElement_UDT(&indexB, path_index_matrices[bin_rule.prod_A], start, index.middle)) != GrB_SUCCESS)
            {
                // If haven't found such a piece of the path, then continue.
                if (info != GrB_NO_VALUE)
                {
                    LG_FREE_WORK;
                    return info;
                }

                continue;
            }
            if ((info = GrB_Matrix_extractElement_UDT(&indexC, path_index_matrices[bin_rule.prod_B], index.middle, end)) != GrB_SUCCESS)
            {
                if (info != GrB_NO_VALUE)
                {
                    LG_FREE_WORK;
                    return info;
                }
                continue;
            }

            // Height compliance check
            int32_t max_height = (indexB.height > indexC.height ? indexB.height : indexC.height);
            if (index.height != max_height + 1)
            {
                continue;
            }

            Path left, right;
            // If didn't find the path, try the other rules.
            if ((info = LAGraph_CFL_extract_single_path_internal(&left, start, index.middle, bin_rule.prod_A, adj_matrices, path_index_matrices, terms_count, nonterms_count, rules, rules_count, msg)) != GrB_SUCCESS)
            {
                if (info == GrB_NO_VALUE)
                {
                    continue;
                }
                LG_FREE_WORK;
                return info;
            }
            if ((info = LAGraph_CFL_extract_single_path_internal(&right, index.middle, end, bin_rule.prod_B, adj_matrices, path_index_matrices, terms_count, nonterms_count, rules, rules_count, msg)) != GrB_SUCCESS)
            {
                if (info == GrB_NO_VALUE)
                {
                    LG_TRY(LAGraph_Free((void **)&left.edges, msg));
                    continue;
                }
                LG_TRY(LAGraph_Free((void **)&left.edges, msg));
                LG_FREE_WORK;
                return info;
            }

            output->len = left.len + right.len;

            LG_TRY(LAGraph_Calloc((void **)&output->edges, output->len, sizeof(Edge), msg));

            memcpy(output->edges, left.edges, left.len * sizeof(Edge));
            memcpy(output->edges + left.len, right.edges, right.len * sizeof(Edge));
            LG_TRY(LAGraph_Free((void **)&left.edges, msg));
            LG_TRY(LAGraph_Free((void **)&right.edges, msg));
            LG_FREE_WORK;
            return GrB_SUCCESS;
        }

        // If couldn't find rules for outputting an path,
        // then the path were looking for doesn't match the rules
        LG_FREE_WORK;
        ADD_TO_MSG(msg_len, "The extracted path does not match the input grammar.");
        return GrB_NO_VALUE;
    }
    // Such a path doesn't exists - return an empty path and GrB_NO_VALUE
    else if (info == GrB_NO_VALUE)
    {
        LG_FREE_WORK;
        return GrB_NO_VALUE;
    }
    // Return some other error
    LG_FREE_WORK;
    return info;
}
