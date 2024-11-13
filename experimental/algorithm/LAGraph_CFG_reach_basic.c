#define LG_FREE_WORK                                                                     \
    do {                                                                                 \
        for (size_t i = 0; i < T_size; i++) {                                            \
            GrB_free(&T[i]);                                                             \
        }                                                                                \
        free(nnz);                                                                       \
    } while (0)

#define LG_FREE_ALL                                                                      \
    { LG_FREE_WORK; }

#include "LG_internal.h"
#include <LAGraphX.h>

// Determines the set of vertex pairs (u, v) in graph (which represented by adjacency
// matrices) such that there is a path from u to v whose edge labels form a word from the
// language generated by the context-free grammar (which represented by `rules`).
//
// Terminals and non-terminals enumerated by integers starting with zero.
// Start nonterminal is nonterminal with index 0
GrB_Info LAGraph_CFG_reach_basic(
    // output
    GrB_Matrix *output, // Matrix with result.
                        // (i, j) = true <=> there is a path
                        // from i to j nodes that form a word derivable
                        // from the start nonterminal of the specified CFG.
    /// input
    const GrB_Matrix
        *adj_matrices, // Array of adjacency matrices which represent graph.
                       // Len of that array equal to count of terminals `terms_count`.
                       //
                       // adj_matrices[t]: (i, j) == 1 <=> there is edge between
                       // i and j nodes with label of terminal with index 't'
                       // t - number in range [0, terms_count - 1]

    size_t terms_count,             // count of terminals in CFG
    size_t nonterms_count,          // count of non-terminals in CFG
    const LAGraph_rule_WCNF *rules, // rules of CFG
    size_t rules_count,             // count of rules
    char *msg) {
    // Declare workspace and clear the msg string, if not NULL
    GrB_Matrix T[nonterms_count];
    size_t T_size = 0; // Variable for correct free
    uint64_t *nnz = NULL;
    LG_CLEAR_MSG;

    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "Count of terms must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "Count of nonterms must be greater than zero.")
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "Count of rules must be greater than zero.");

    GrB_Index n;
    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    // Create nonterms matrices
    for (size_t i = 0; i < nonterms_count; i++) {
        GRB_TRY(GrB_Matrix_new(&T[i], GrB_BOOL, n, n));
        T_size++;
    }

    // Arrays for processing rules
    size_t eps_rules[rules_count], eps_rules_count = 0;   // [Variable -> eps]
    size_t term_rules[rules_count], term_rules_count = 0; // [Variable -> term]
    size_t bin_rules[rules_count], bin_rules_count = 0;   // [Variable -> AB]

    // Process rules
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_WCNF rule = rules[i];

        // [Variable -> eps]
        if (rule.prod_A == -1 && rule.prod_B == -1) {
            eps_rules[eps_rules_count++] = i;
            continue;
        }

        // [Variable -> term]
        if (rule.prod_A != -1 && rule.prod_B == -1) {
            term_rules[term_rules_count++] = i;
            continue;
        }

        // [Variable -> A B]
        if (rule.prod_A != -1 && rule.prod_B != -1) {
            bin_rules[bin_rules_count++] = i;
            continue;
        }

        // [Variable -> _ B]
        LG_ASSERT_MSG(true, GrB_INVALID_VALUE,
                      "Rule with index %d is invalid. [Variable -> _ B] type of rule "
                      "don't accepted.");
    }

    // Rule [Variable -> term]
    for (size_t i = 0; i < term_rules_count; i++) {
        LAGraph_rule_WCNF term_rule = rules[term_rules[i]];

        GrB_Index nnz;
        GrB_Matrix_nvals(&nnz, adj_matrices[term_rule.prod_A]);

        GrB_Index row[nnz], col[nnz];
        bool values[nnz];

        GrB_Matrix_extractTuples(row, col, values, &nnz, adj_matrices[term_rule.prod_A]);
        for (size_t j = 0; j < nnz; j++) {
#ifdef DEBUG
            printf("[TERM] SET ELEMENT [TRUE], NONTERM: %d, ROW: %ld, COL: %ld\n",
                   term_rule.nonterm, row[j], col[j]);
#endif
            GrB_Matrix_setElement(T[term_rule.nonterm], true, row[j], col[j]);
        }
    }

    // Rule [Variable -> eps]
    for (size_t i = 0; i < eps_rules_count; i++) {
        LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];

        for (size_t j = 0; j < n; j++) {
#ifdef DEBUG
            printf("[EPS] SET ELEMENT [TRUE], NONTERM: %d, INDEX: %ld\n",
                   eps_rule.nonterm, j);
#endif
            GrB_Matrix_setElement(T[eps_rule.nonterm], true, j, j);
        }
    }

    // Rule [Variable -> Variable1 Variable2]
    nnz = calloc(nonterms_count, sizeof(uint64_t));
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            GRB_TRY(GrB_mxm(T[bin_rule.nonterm], GrB_NULL, GxB_LOR_BOOL,
                            GxB_ANY_PAIR_BOOL, T[bin_rule.prod_A], T[bin_rule.prod_B],
                            GrB_NULL));

            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, T[bin_rule.nonterm]));

            changed = changed | (nnz[bin_rule.nonterm] != new_nnz);
#ifdef DEBUG
            printf("[TERM1 TERM2] MULTIPLY, S: %d, A: %d, B: %d, "
                   "I: %ld\n",
                   bin_rule.nonterm, bin_rule.prod_A, bin_rule.prod_B, i);
#endif
            nnz[bin_rule.nonterm] = new_nnz;
        }
    }

#ifdef DEBUG
    for (size_t i = 0; i < nonterms_count; i++) {
        printf("MATRIX WITH INDEX %ld:\n", i);
        GxB_print(T[i], 5);
    }
#endif

    GrB_Matrix_dup(output, T[0]);
    LG_FREE_WORK;
    return GrB_SUCCESS;
}