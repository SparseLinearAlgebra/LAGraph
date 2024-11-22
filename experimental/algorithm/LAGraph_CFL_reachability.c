// Code is based on the "A matrix-based CFPQ algorithm" described in the following paper:
//  * Rustam Azimov, Semyon Grigorev, "Context-Free Path Querying Using Linear Algebra"
//  * URL: https://disser.spbu.ru/files/2022/disser_azimov.pdf

#define ERROR_RULE(msg)                                                                  \
    {                                                                                    \
        LG_ASSERT_MSGF(false, GrB_INVALID_VALUE, "Rule with index %ld is invalid. " msg, \
                       i);                                                               \
    }

#define LG_FREE_WORK                                                                     \
    do {                                                                                 \
        free(nnz);                                                                       \
    } while (0)

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        LG_FREE_WORK;                                                                    \
                                                                                         \
        for (size_t i = 0; i < T_size; i++) {                                            \
            GrB_free(&T[i]);                                                             \
        }                                                                                \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

// LAGraph_CFL_reachability: Context-Free Language Reachability Matrix-Based
// Algorithm
//
// Determines the set of vertex pairs (u, v) in graph (which represented by adjacency
// matrices) such that there is a path from u to v whose edge labels form a word from the
// language generated by the context-free grammar (which represented by `rules`).
//
// Terminals and non-terminals enumerated by integers starting with zero.
// Start nonterminal is nonterminal with index 0
//
// Example:
//
// Graph:
// ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
// │ 0 ├───► 1 ├───► 2 ├───► 3 ├───► 4 │
// └───┘ a └─┬─┘ a └───┘ b └─▲─┘ b └───┘
//           │               │
//           │     ┌───┐     │
//           └─────► 5 ├─────┘
//             a   └───┘   b
//
// Grammar: S -> aSb | ab
//
// There are paths from node [1] to node [3] that form the word "ab" ([1]-a->[2]-b->[3]
// and [1]-a->[5]-b->[3]). "ab" word is in language generated by our context-free grammar,
// so pair (1, 3) of nodes will be in result.
//
// Notice: it doesn't matter how many paths exist from node [A] to node [B] that form word
// of language. If at least on path exists - pair ([A], [B]) of nodes will be in result.
//
// In contrast, path from [1] node to [4] node form word "abb" ([1]-a->[2]-b->[3]-b->[4]).
// "aab" word is not in language, so pair (1, 4) of nodes will not be in result.
//
// With these graph and grammar we get result:
// (0, 4) - because exists path that form word of our language (0-1-2-3-4 or 0-1-5-3-4)
// (1, 3) - because exists path (1-2-3 or 1-5-3)
GrB_Info LAGraph_CFL_reachability(
    // output
    GrB_Matrix *outputs, // Array of matrices with result.
                         // Size of array must be equal to nonterms_count
                         //
                         // outputs[k]: (i, j) = true <=> there is a path from i to j
                         // nodes whose edge labels form a word derivable
                         // from the start nonterminal 'k' of the specified CFG.
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
    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "Outputs array is null.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "Rules array is null.");
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "Adjacency matrix array is null.");

    for (size_t i = 0; i < terms_count; i++) {
        LG_ASSERT_MSGF(adj_matrices[i] != NULL, GrB_NULL_POINTER,
                       "Adjacency matrix with index %ld is null.", i);
    }

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

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // Check range on rules
        if (rule.nonterm < 0 || rule.nonterm >= nonterms_count) {
            ERROR_RULE("Nonterm must be in range [0, nonterms_count).");
        }

        // [Variable -> eps]
        if (is_rule_eps) {
            eps_rules[eps_rules_count++] = i;

            continue;
        }

        // [Variable -> term]
        if (is_rule_term) {
            term_rules[term_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= terms_count) {
                ERROR_RULE("Term must be in range [-1, nonterms_count)");
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin) {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= nonterms_count || rule.prod_B < -1 ||
                rule.prod_B >= nonterms_count) {
                ERROR_RULE("Nonterm must be in range [0, nonterms_count).");
            }

            continue;
        }

        // [Variable -> _ B]
        ERROR_RULE("[Variable -> _ B] type of rule is not acceptable.");
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

    for (size_t i = 0; i < nonterms_count; i++) {
        outputs[i] = T[i];
    }
    LG_FREE_WORK;
    return GrB_SUCCESS;
}