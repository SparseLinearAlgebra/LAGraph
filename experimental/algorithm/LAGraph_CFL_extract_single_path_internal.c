#define LG_FREE_WORK                              \
    {                                             \
        LAGraph_Free((void **)&eps_rules, NULL);  \
        LAGraph_Free((void **)&term_rules, NULL); \
        LAGraph_Free((void **)&bin_rules, NULL);  \
    }

#define LG_FREE_ALL                                 \
    {                                               \
        LAGraph_Free((void **)&output->edges, NULL); \
        output->len = 0;                            \
        LG_FREE_WORK;                               \
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

GrB_Info LAGraph_CFL_extract_single_path_internal(
    // Output
    Path *output,
    // Input
    GrB_Index start,
    GrB_Index end,
    int32_t nonterm,
    const GrB_Matrix *adj_matrices,
    const GrB_Matrix *path_index_matrices,
    int64_t terms_count,            // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
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

    // Internal function: assumes all inputs have been validated by the public wrapper

    // Process rules
    for (int64_t i = 0; i < rules_count; i++)
    {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // [Variable -> eps]
        if (is_rule_eps)
        {
            eps_rules[eps_rules_count++] = i;
            continue;
        }

        // [Variable -> term]
        if (is_rule_term)
        {
            term_rules[term_rules_count++] = i;
            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin)
        {
            bin_rules[bin_rules_count++] = i;
            continue;
        }
    }

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
            ADD_TO_MSG("The extracted path does not match the input grammar.");
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
        ADD_TO_MSG("The extracted path does not match the input grammar.");
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
