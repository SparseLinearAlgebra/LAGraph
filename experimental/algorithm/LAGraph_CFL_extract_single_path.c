#define LG_FREE_WORK                              \
    {                                             \
        LAGraph_Free((void **)&eps_rules, NULL);  \
        LAGraph_Free((void **)&term_rules, NULL); \
        LAGraph_Free((void **)&bin_rules, NULL);  \
    }

#define LG_FREE_ALL                                 \
    {                                               \
        LAGraph_Free((void **)&output->path, NULL); \
        output->len = 0;                            \
        LG_FREE_WORK;                               \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

#define ERROR_RULE(msg, i)                                                  \
    {                                                                       \
        LG_ASSERT_MSGF(false, GrB_INVALID_VALUE,                            \
                       "Rule with index %" PRId64 " is invalid. ", msg, i); \
    }

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

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                     \
    {                                                        \
        rule.len_indexes_str += snprintf(                    \
            rule.indexes_str + rule.len_indexes_str,         \
            LAGRAPH_MSG_LEN - rule.len_indexes_str,          \
            rule.count == 0 ? "%" PRId64 : ", %" PRId64, i); \
        rule.count++;                                        \
    }

GrB_Info LAGraph_CFL_extract_single_path(
    // Output
    Path *output,
    // Input
    GrB_Index start,
    GrB_Index end,
    int32_t nonterm,
    const GrB_Matrix *adj_matrices,
    const GrB_Matrix *T,
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
    output->path = NULL;
    // Arrays for processing rules
    size_t *eps_rules = NULL, eps_rules_count = 0;   // [Variable -> eps]
    size_t *term_rules = NULL, term_rules_count = 0; // [Variable -> term]
    size_t *bin_rules = NULL, bin_rules_count = 0;   // [Variable -> AB]
    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "The number of terminals must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "The number of non-terminals must be greater than zero.");
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "The number of rules must be greater than zero.");
    LG_ASSERT_MSG(nonterm < nonterms_count, GrB_INVALID_VALUE,
                  "The start non-terminal must be no greater than the number of non-terminals.");
    LG_ASSERT_MSG(T != NULL, GrB_NULL_POINTER, "The T array cannot be null.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "The rules array cannot be null.");
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "The adjacency matrices array cannot be null.");

    // Find null adjacency matrices
    bool found_null = false;
    for (int64_t i = 0; i < terms_count; i++)
    {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null)
        {
            ADD_TO_MSG("Adjacency matrices with these indexes are null: ");
            ADD_TO_MSG("%" PRId64, i);
        }
        else
        {
            ADD_TO_MSG("%" PRId64, i);
        }

        found_null = true;
    }

    if (found_null)
    {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    // Find null T matrices
    found_null = false;
    for (int64_t i = 0; i < nonterms_count; i++)
    {
        if (T[i] != NULL)
            continue;

        if (!found_null)
        {
            ADD_TO_MSG("T matrices with these indexes are null: ");
            ADD_TO_MSG("%" PRId64, i);
        }
        else
        {
            ADD_TO_MSG("%" PRId64, i);
        }

        found_null = true;
    }
    if (found_null)
    {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    LG_TRY(LAGraph_Calloc((void **)&eps_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&term_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&bin_rules, rules_count, sizeof(size_t), msg));
    // Process rules
    typedef struct
    {
        size_t count;
        size_t len_indexes_str;
        char indexes_str[LAGRAPH_MSG_LEN];
    } rule_error_s;
    rule_error_s term_err = {0};
    rule_error_s nonterm_err = {0};
    rule_error_s invalid_err = {0};
    for (int64_t i = 0; i < rules_count; i++)
    {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // Check that all rules are well-formed
        if (rule.nonterm < 0 || rule.nonterm >= nonterms_count)
        {
            ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
        }

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

            if (rule.prod_A < -1 || rule.prod_A >= terms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin)
        {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= nonterms_count || rule.prod_B < -1 ||
                rule.prod_B >= nonterms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
            }

            continue;
        }

        // [Variable -> _ B]
        ADD_INDEX_TO_ERROR_RULE(invalid_err, i);
    }

    if (term_err.count + nonterm_err.count + invalid_err.count > 0)
    {
        ADD_TO_MSG("Count of invalid rules: %" PRId64 ".\n",
                   (int64_t)(term_err.count + nonterm_err.count + invalid_err.count));

        if (nonterm_err.count > 0)
        {
            ADD_TO_MSG("Non-terminals must be in range [0, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indexes_str)
        }
        if (term_err.count > 0)
        {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indexes_str)
        }
        if (invalid_err.count > 0)
        {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indexes_str)
        }

        LG_FREE_ALL;
        return GrB_INVALID_VALUE;
    }
    PathIndex index;
    GrB_Info info = GrB_Matrix_extractElement_UDT(&index, T[nonterm], start, end);
    if (info == GrB_SUCCESS) // Such a path exists
    {
        if (index.height == 1)
        {
            if (start == end) // Height = 1 and start = end is an empty eps-path
            {
                for (size_t i = 0; i < eps_rules_count; i++)
                {
                    LAGraph_rule_WCNF term_rule = rules[eps_rules[i]];
                    if (term_rule.nonterm == nonterm)
                    {
                        LG_FREE_WORK;
                        return GrB_SUCCESS;
                    }
                }
            }
            // Height = 1 and different vertices is a term-path
            for (int64_t i = 0; i < terms_count; i++)
            {
                bool edge;
                if (GrB_Matrix_extractElement_BOOL(&edge, adj_matrices[i], start, end) == GrB_SUCCESS)
                {
                    for (size_t j = 0; j < term_rules_count; j++)
                    {
                        LAGraph_rule_WCNF term_rule = rules[term_rules[j]];
                        if (term_rule.nonterm == nonterm && term_rule.prod_A == i)
                        {
                            LG_TRY(LAGraph_Calloc((void **)&output->path, 1, sizeof(Edge), msg));
                            output->len = 1;
                            output->path[0] = (Edge){start, i, end};
                            LG_FREE_WORK;
                            return GrB_SUCCESS;
                        }
                    }
                }
            }
            // If couldn't find rules for outputting an empty or terminal path,
            // then the path were looking for doesn't match the rules
            LG_FREE_ALL;
            ADD_TO_MSG("The extracted path does not match the input grammar.");
            return GrB_NO_VALUE;
        }
        // Rules of the form Nonterm -> Nonterm * Nonterm are traversed recursively and merged
        for (size_t i = 0; i < bin_rules_count; i++)
        {
            LAGraph_rule_WCNF term_rule = rules[bin_rules[i]];
            if (term_rule.nonterm == nonterm)
            {
                PathIndex indexB, indexC;
                if ((info = GrB_Matrix_extractElement_UDT(&indexB, T[term_rule.prod_A], start, index.middle)) != GrB_SUCCESS)
                {
                    // If haven't found such a piece of the path, then continue.
                    if (info != GrB_NO_VALUE)
                    {
                        return info;
                    }

                    continue;
                }
                if ((info = GrB_Matrix_extractElement_UDT(&indexC, T[term_rule.prod_B], index.middle, end)) != GrB_SUCCESS)
                {
                    if (info != GrB_NO_VALUE)
                    {
                        return info;
                    }
                    continue;
                }

                if (!((indexB.height == 0 && indexB.middle == 0) || (indexC.height == 0 && indexC.middle == 0)))
                {
                    int32_t maxH = (indexB.height > indexC.height ? indexB.height : indexC.height);
                    if (index.height == maxH + 1)
                    {
                        Path left, right;
                        LG_TRY(LAGraph_CFL_extract_single_path(&left, start, index.middle, term_rule.prod_A, adj_matrices, T, terms_count, nonterms_count, rules, rules_count, msg));
                        LG_TRY(LAGraph_CFL_extract_single_path(&right, index.middle, end, term_rule.prod_B, adj_matrices, T, terms_count, nonterms_count, rules, rules_count, msg));

                        output->len = left.len + right.len;

                        LG_TRY(LAGraph_Calloc((void **)&output->path, output->len, sizeof(Edge), msg));

                        memcpy(output->path, left.path, left.len * sizeof(Edge));
                        memcpy(output->path + left.len, right.path, right.len * sizeof(Edge));
                        LG_TRY(LAGraph_Free((void **)&left.path, msg));
                        LG_TRY(LAGraph_Free((void **)&right.path, msg));
                        LG_FREE_WORK;
                        return GrB_SUCCESS;
                    }
                }
            }
        }
        // If couldn't find rules for outputting an path,
        // then the path were looking for doesn't match the rules
        LG_FREE_ALL;
        ADD_TO_MSG("The extracted path does not match the input grammar.");
        return GrB_NO_VALUE;
    }
    // Such a path doesn't exists - return an empty path and GrB_NO_VALUE
    else if (info == GrB_NO_VALUE)
    {
        LG_FREE_ALL;
        return GrB_NO_VALUE;
    }
    // Return some other error
    LG_FREE_ALL;
    return info;
}
