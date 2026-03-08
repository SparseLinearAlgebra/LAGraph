#include "LG_internal.h"
#include <LAGraphX.h>

// Structure for collecting rule check errors
// Used by grammar check macros to accumulate error information
typedef struct
{
    size_t count;
    size_t len_indexes_str;
    char indexes_str[LAGRAPH_MSG_LEN];
} LAGraph_rule_error_s;

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                     \
    {                                                        \
        rule.len_indexes_str += snprintf(                    \
            rule.indexes_str + rule.len_indexes_str,         \
            LAGRAPH_MSG_LEN - rule.len_indexes_str,          \
            rule.count == 0 ? "%" PRId64 : ", %" PRId64, i); \
        rule.count++;                                        \
    }

// Checks the input CFL grammar
GrB_Info LAGraph_CFL_check_grammar(int64_t terms_count, int64_t nonterms_count, int64_t rules_count, const LAGraph_rule_WCNF *rules, char *msg)
{
    LG_CLEAR_MSG;
    size_t msg_len = 0;
    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "The number of terminals must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "The number of non-terminals must be greater than zero.");
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "The number of rules must be greater than zero.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "The rules array cannot be null.");

    LAGraph_rule_error_s term_err = {0};
    LAGraph_rule_error_s nonterm_err = {0};
    LAGraph_rule_error_s invalid_err = {0};

    for (int64_t i = 0; i < rules_count; i++)
    {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = (rule.prod_A == -1 && rule.prod_B == -1);
        bool is_rule_term = (rule.prod_A != -1 && rule.prod_B == -1);
        bool is_rule_bin = (rule.prod_A != -1 && rule.prod_B != -1);

        /* Check that all rules are well-formed */
        if (rule.nonterm < 0 || rule.nonterm >= nonterms_count)
        {
            ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
        }

        /* [Variable -> eps]  */
        if (is_rule_eps)
        {
            continue;
        }

        /* [Variable -> term] */
        if (is_rule_term)
        {
            if (rule.prod_A < -1 || rule.prod_A >= terms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }
            continue;
        }

        /* [Variable -> A B] */
        if (is_rule_bin)
        {
            if (rule.prod_A < -1 || rule.prod_A >= nonterms_count ||
                rule.prod_B < -1 || rule.prod_B >= nonterms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
            }
            continue;
        }

        /* [Variable -> _ B] */
        ADD_INDEX_TO_ERROR_RULE(invalid_err, i);
    }

    if (term_err.count + nonterm_err.count + invalid_err.count > 0)
    {
        ADD_TO_MSG("Count of invalid rules: %" PRId64 ".\n",
                   (int64_t)(term_err.count + nonterm_err.count + invalid_err.count));

        if (nonterm_err.count > 0)
        {
            ADD_TO_MSG("Non-terminals must be in range [0, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indexes_str);
        }
        if (term_err.count > 0)
        {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indexes_str);
        }
        if (invalid_err.count > 0)
        {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indexes_str);
        }
        return GrB_INVALID_VALUE;
    }
    return GrB_SUCCESS;
}

// Checks the input graph for CFL algorithms
GrB_Info LAGraph_CFL_check_graph(const GrB_Matrix *adj_matrices, int64_t terms_count, char *msg)
{
    LG_CLEAR_MSG;
    size_t msg_len = 0;
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "The adjacency matrices array cannot be null.");

    /* Find null adjacency matrices */
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
            ADD_TO_MSG(" %" PRId64, i);
        }

        found_null = true;
    }

    if (found_null)
    {
        return GrB_NULL_POINTER;
    }
    return GrB_SUCCESS;
}

// LAGraph_CFL_check_base_inputs: Checks the input graph and CFL grammar for algorithms working with it
//
// Checks:
//   - The adjacency matrix array is not NULL
//   - There are no entries equal to NULL in the adjacency matrix array
//   - The rule array is not NULL
//   - Valid number of terms/non-terms/rules (> 0)
//   - Correctly formed grammar rules in WCNF format
//
// Parameters:
//   - adj_matrices: adjacency matrix array for terminals
//   - terms_count: number of terminal symbols
//   - nonterms_count: number of non-terminal symbols
//   - rules_count: number of rules in the grammar
//   - rules: array of rules
//
// If an error occurs: adds a message using ADD_TO_MSG,
// returns the corresponding GrB_Info

GrB_Info LAGraph_CFL_check_base_inputs(
    const GrB_Matrix *adj_matrices,
    int64_t terms_count,
    int64_t nonterms_count,
    int64_t rules_count,
    const LAGraph_rule_WCNF *rules,
    char *msg)
{
    LG_CLEAR_MSG;
    size_t msg_len = 0;
    LG_TRY(LAGraph_CFL_check_graph(adj_matrices, terms_count, msg));
    LG_TRY(LAGraph_CFL_check_grammar(terms_count, nonterms_count, rules_count, rules, msg));
    return GrB_SUCCESS;
}
