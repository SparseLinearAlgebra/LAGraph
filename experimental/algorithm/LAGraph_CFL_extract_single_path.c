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
    PathArray *output,
    // Input
    GrB_Index *start,
    GrB_Index *end,
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
    GrB_Info info;
    output->count = 0;
    output->paths = NULL;
    // Initial capacity 1, since most often looking for exactly 1 path with a fixed start and end
    output->capacity = 1;

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
            ADD_TO_MSG("Adjacency matrices with these indexes are null:");
        }
        ADD_TO_MSG(" %" PRId64, i);
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
            ADD_TO_MSG("T matrices with these indexes are null:");
        }
        ADD_TO_MSG(" %" PRId64, i);

        found_null = true;
    }
    if (found_null)
    {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    // Check the rules
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

        // [Variable -> term]
        if (is_rule_term)
        {
            if (rule.prod_A < -1 || rule.prod_A >= terms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }
            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin)
        {
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

    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, adj_matrices[0]));

    GrB_Index start_begin = (start == NULL) ? 0 : *start;
    GrB_Index start_end = (start == NULL) ? n - 1 : *start;
    GrB_Index end_begin = (end == NULL) ? 0 : *end;
    GrB_Index end_end = (end == NULL) ? n - 1 : *end;

    LG_TRY(LAGraph_Malloc((void **)&output->paths, output->capacity, sizeof(Path), msg));

    for (GrB_Index st = start_begin; st <= start_end; st++)
    {
        for (GrB_Index en = end_begin; en <= end_end; en++)
        {
            Path path;

            // Function that extracts one path with fixed start and end
            info = LAGraph_CFL_extract_single_path_internal(&path, st, en, nonterm, adj_matrices, T, terms_count, nonterms_count, rules, rules_count, msg);
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

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
