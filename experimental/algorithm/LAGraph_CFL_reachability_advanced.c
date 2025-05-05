//------------------------------------------------------------------------------
// LAGraph_CFL_reachability.c: Context-Free Language Reachability Matrix-Based
// Algorithm
// ------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State University.

//------------------------------------------------------------------------------

// Code is based on the "A matrix-based CFPQ algorithm" described in the
// following paper: * Rustam Azimov, Semyon Grigorev, "Context-Free Path
// Querying Using Linear Algebra", URL:
// https://disser.spbu.ru/files/2022/disser_azimov.pdf

#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        LAGraph_Free((void **)&nnzs, msg);                                               \
        GrB_free(&true_scalar);                                                          \
        GrB_free(&identity_matrix);                                                      \
        LAGraph_Free((void **)&T, msg);                                                  \
        LAGraph_Free((void **)&indexes, msg);                                            \
    }

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        for (size_t i = 0; i < nonterms_count; i++) {                                    \
            GrB_free(&T[i]);                                                             \
        }                                                                                \
                                                                                         \
        LG_FREE_WORK;                                                                    \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

#define ERROR_RULE(msg)                                                                  \
    {                                                                                    \
        LG_ASSERT_MSGF(false, GrB_INVALID_VALUE, "Rule with index %ld is invalid. " msg, \
                       i);                                                               \
    }

#define ADD_TO_MSG(...)                                                                  \
    {                                                                                    \
        if (msg_len == 0) {                                                              \
            msg_len +=                                                                   \
                snprintf(msg, LAGRAPH_MSG_LEN,                                           \
                         "LAGraph failure (file %s, line %d): ", __FILE__, __LINE__);    \
        }                                                                                \
        if (msg_len < LAGRAPH_MSG_LEN) {                                                 \
            msg_len += snprintf(msg + msg_len, LAGRAPH_MSG_LEN - msg_len, __VA_ARGS__);  \
        }                                                                                \
    }

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                                                 \
    {                                                                                    \
        rule.len_indexes_str += snprintf(rule.indexes_str + rule.len_indexes_str,        \
                                         LAGRAPH_MSG_LEN - rule.len_indexes_str,         \
                                         rule.count == 0 ? "%ld" : ", %ld", i);          \
        rule.count++;                                                                    \
    }

#define IS_ISO(matrix, msg)                                                              \
    {                                                                                    \
        GxB_Matrix_iso(&iso_flag, matrix);                                               \
        GRB_TRY(GrB_Matrix_nvals(&new_nnz, matrix));                                     \
        if (!iso_flag && new_nnz) {                                                      \
            GxB_print(matrix, 1);                                                        \
            printf(msg);                                                                 \
        }                                                                                \
    }

#define SKIP_IF_NULL(matrix)                                                             \
    GrB_Matrix_nvals(&new_nnz, matrix);                                                  \
    if (new_nnz == 0) {                                                                  \
        continue;                                                                        \
    }

#define TIMER_START()                                                                    \
    {                                                                                    \
        start_time = LAGraph_WallClockTime();                                            \
    }

#define TIMER_STOP(label, accumulator)                                                   \
    {                                                                                    \
        end_time = LAGraph_WallClockTime();                                              \
        printf("%s %.3fs\n", label, end_time - start_time);                              \
        if (accumulator != NULL) {                                                       \
            *(accumulator) += (end_time - start_time);                                   \
        }                                                                                \
    }

// LAGraph_CFL_reachability: Context-Free Language Reachability Matrix-Based Algorithm
//
// This function determines the set of vertex pairs (u, v) in a graph (represented by
// adjacency matrices) such that there is a path from u to v, where the edge labels form a
// word from the language generated by the context-free grammar (represented by `rules`).
//
// Terminals and non-terminals are enumerated by integers starting from zero.
// The start non-terminal is the non-terminal with index 0.
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
// generated by our context-free grammar, so the pairs (1, 3) and (1, 2) will be included
// in the result.
//
// Note: It doesn't matter how many paths exist from node [A] to node [B] that form a word
// in the language. If at least one path exists, the pair ([A], [B]) will be included in
// the result.
//
// In contrast, the path from node [1] to node [4] forms the word "abb"
// ([1]-a->[2]-b->[3]-b->[4]) and the word "abbb" ([1]-a->[5]-b->[2]-b->[3]-b->[4]).
// The words "aab" and "abbb" are not in the language, so the pair (1, 4) will not be
// included in the result.
//
// With this graph and grammar, we obtain the following results:
// (0, 4) - because there exists a path (0-1-2-3-4) that forms the word "aabb"
// (1, 3) - because there exists a path (1-2-3) that forms "ab"
// (1, 2) - because there exists a path (1-5-2) that forms the word "ab"
// (0, 3) - because there exists a path (0-1-5-2-3) that forms the word "aabb"
GrB_Info LAGraph_CFL_reachability_adv(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) = true if and only if there is a path
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
    int32_t terms_count,            // The total number of terminal symbols in the CFG.
    int32_t nonterms_count, // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    size_t rules_count,             // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
) {
    // Declare workspace and clear the msg string, if not NULL
    GrB_Matrix *T;
    GrB_Matrix *delta_matrices;
    GrB_Matrix *matrices;
    GrB_Matrix *temp_matrices;
    bool t_empty_flags[nonterms_count]; // t_empty_flags[i] == true <=> T[i] is empty
    GrB_Matrix identity_matrix = NULL;
    uint64_t *nnzs = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    bool iso_flag = false;
    GrB_Index *indexes = NULL;

    GrB_Scalar true_scalar;
    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);

    LG_TRY(LAGraph_Calloc((void **)&T, nonterms_count, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&delta_matrices, nonterms_count, sizeof(GrB_Matrix),
                          msg));
    LG_TRY(LAGraph_Calloc((void **)&matrices, nonterms_count, sizeof(GrB_Matrix), msg));
    LG_TRY(
        LAGraph_Calloc((void **)&temp_matrices, nonterms_count, sizeof(GrB_Matrix), msg));

    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "The number of terminals must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "The number of non-terminals must be greater than zero.");
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "The number of rules must be greater than zero.");
    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "The outputs array cannot be null.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "The rules array cannot be null.");
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "The adjacency matrices array cannot be null.");

    // Find null adjacency matrices
    bool found_null = false;
    for (int32_t i = 0; i < terms_count; i++) {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null) {
            ADD_TO_MSG("Adjacency matrices with these indexes are null: ");
            ADD_TO_MSG("%d", i);
        } else {
            ADD_TO_MSG(", %d", i);
        }

        found_null = true;
    }

    if (found_null) {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    GrB_Index n;
    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    // Create nonterms matrices
    for (int32_t i = 0; i < nonterms_count; i++) {
        GRB_TRY(GrB_Matrix_new(&T[i], GrB_BOOL, n, n));
        GRB_TRY(GrB_Matrix_new(&delta_matrices[i], GrB_BOOL, n, n));
        GRB_TRY(GrB_Matrix_new(&matrices[i], GrB_BOOL, n, n));
        GRB_TRY(GrB_Matrix_new(&temp_matrices[i], GrB_BOOL, n, n));
        t_empty_flags[i] = true;
    }

    // Arrays for processing rules
    size_t eps_rules[rules_count], eps_rules_count = 0;   // [Variable -> eps]
    size_t term_rules[rules_count], term_rules_count = 0; // [Variable -> term]
    size_t bin_rules[rules_count], bin_rules_count = 0;   // [Variable -> AB]

    // Process rules
    typedef struct {
        size_t count;
        size_t len_indexes_str;
        char indexes_str[LAGRAPH_MSG_LEN];
    } rule_error_s;
    rule_error_s term_err = {0};
    rule_error_s nonterm_err = {0};
    rule_error_s invalid_err = {0};
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // Check that all rules are well-formed
        if (rule.nonterm < 0 || rule.nonterm >= nonterms_count) {
            ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
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
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin) {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= nonterms_count || rule.prod_B < -1 ||
                rule.prod_B >= nonterms_count) {
                ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
            }

            continue;
        }

        // [Variable -> _ B]
        ADD_INDEX_TO_ERROR_RULE(invalid_err, i);
    }

    if (term_err.count + nonterm_err.count + invalid_err.count > 0) {
        ADD_TO_MSG("Count of invalid rules: %ld.\n",
                   term_err.count + nonterm_err.count + invalid_err.count);

        if (nonterm_err.count > 0) {
            ADD_TO_MSG("Non-terminals must be in range [0, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indexes_str)
        }
        if (term_err.count > 0) {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indexes_str)
        }
        if (invalid_err.count > 0) {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indexes_str)
        }

        LG_FREE_ALL;
        return GrB_INVALID_VALUE;
    }

    // Rule [Variable -> term]
    for (size_t i = 0; i < term_rules_count; i++) {
        LAGraph_rule_WCNF term_rule = rules[term_rules[i]];
        GrB_Index adj_matrix_nnz = 0;
        GRB_TRY(GrB_Matrix_nvals(&adj_matrix_nnz, adj_matrices[term_rule.prod_A]));

        if (adj_matrix_nnz == 0) {
            continue;
        }

        GxB_eWiseUnion(delta_matrices[term_rule.nonterm], GrB_NULL, GrB_NULL,
                       GxB_PAIR_BOOL, delta_matrices[term_rule.nonterm], true_scalar,
                       adj_matrices[term_rule.prod_A], true_scalar, GrB_NULL);

        t_empty_flags[term_rule.nonterm] = false;

#ifdef DEBUG_CFL_REACHBILITY
        GxB_Matrix_iso(&iso_flag, T[term_rule.nonterm]);
        printf("[TERM] eWiseUnion: NONTERM: %d (ISO: %d)\n", term_rule.nonterm, iso_flag);
#endif
    }

    GrB_Vector v_diag;
    GRB_TRY(GrB_Vector_new(&v_diag, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v_diag, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&identity_matrix, v_diag, 0));
    GRB_TRY(GrB_free(&v_diag));

    // Rule [Variable -> eps]
    for (size_t i = 0; i < eps_rules_count; i++) {
        LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];

        GxB_eWiseUnion(delta_matrices[eps_rule.nonterm], GrB_NULL, GxB_PAIR_BOOL,
                       GxB_PAIR_BOOL, delta_matrices[eps_rule.nonterm], true_scalar,
                       identity_matrix, true_scalar, GrB_NULL);

        t_empty_flags[eps_rule.nonterm] = false;

#ifdef DEBUG_CFL_REACHBILITY
        GxB_Matrix_iso(&iso_flag, T[eps_rule.nonterm]);
        printf("[EPS] eWiseUnion: NONTERM: %d (ISO: %d)\n", eps_rule.nonterm, iso_flag);
#endif
    }

    // Rule [Variable -> Variable1 Variable2]
    LG_TRY(LAGraph_Calloc((void **)&nnzs, nonterms_count, sizeof(uint64_t), msg));

    double start_time, end_time;
    bool changed = true;
    size_t iteration = 0;
    double mxm1 = 0.0;
    double wise1 = 0.0;
    double mxm2 = 0.0;
    double wise2 = 0.0;
    double rsub = 0.0;
    GrB_Index new_nnz;
    while (changed) {
        iteration++;
        changed = false;

        printf("\n--- ITERATARION %d ---\n", iteration);

        for (size_t i = 0; i < nonterms_count; i++) {
            GRB_TRY(GrB_Matrix_new(&temp_matrices[i], GrB_BOOL, n, n));
        }

        TIMER_START();
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            SKIP_IF_NULL(matrices[bin_rule.prod_A]);
            SKIP_IF_NULL(delta_matrices[bin_rule.prod_B]);

            GrB_Index left_nnz;
            GrB_Index right_nnz;

            GRB_TRY(GrB_mxm(temp_matrices[bin_rule.nonterm], GrB_NULL, GrB_LOR,
                            GrB_LOR_LAND_SEMIRING_BOOL, matrices[bin_rule.prod_A],
                            delta_matrices[bin_rule.prod_B], GrB_NULL));

            IS_ISO(temp_matrices[bin_rule.nonterm], "ALERT");
            IS_ISO(matrices[bin_rule.prod_A], "ALERT");
            IS_ISO(delta_matrices[bin_rule.prod_B], "ALERT");
        }
        TIMER_STOP("MXM 1", &mxm1);

        TIMER_START()
        for (size_t i = 0; i < nonterms_count; i++) {
            SKIP_IF_NULL(delta_matrices[i]);

            GrB_Matrix_nvals(&new_nnz, matrices[i]);
            if (new_nnz == 0) {
                GrB_Matrix_apply(matrices[i], GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL,
                                 delta_matrices[i], GrB_NULL);
                IS_ISO(matrices[i], "ALERT WISE 1");
                continue;
            }

            GrB_eWiseAdd(matrices[i], GrB_NULL, GrB_NULL, GxB_ANY_BOOL, matrices[i],
                         delta_matrices[i], GrB_NULL);

            IS_ISO(matrices[i], "ALERT WISE 1");
        }
        TIMER_STOP("WISE 1", &wise1);

        TIMER_START()
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            SKIP_IF_NULL(delta_matrices[bin_rule.prod_A]);
            SKIP_IF_NULL(matrices[bin_rule.prod_B]);

            GRB_TRY(GrB_mxm(temp_matrices[bin_rule.nonterm], GrB_NULL, GxB_ANY_BOOL,
                            GxB_ANY_PAIR_BOOL, delta_matrices[bin_rule.prod_A],
                            matrices[bin_rule.prod_B], GrB_NULL))
        }
        TIMER_STOP("MXM 2", &mxm2);

        TIMER_START();
        for (size_t i = 0; i < nonterms_count; i++) {
            GrB_Matrix_new(&delta_matrices[i], GrB_BOOL, n, n);
            SKIP_IF_NULL(temp_matrices[i]);

            GrB_Matrix_apply(delta_matrices[i], GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL,
                             temp_matrices[i], GrB_NULL);

            IS_ISO(delta_matrices[i], "WISE 2");
        }
        TIMER_STOP("WISE 2 (copy)", &wise2);

        TIMER_START();
        for (size_t i = 0; i < nonterms_count; i++) {
            SKIP_IF_NULL(delta_matrices[i]);

            GrB_Matrix_apply(delta_matrices[i], matrices[i], GrB_NULL, GrB_IDENTITY_BOOL,
                             delta_matrices[i], GrB_DESC_C);
            IS_ISO(delta_matrices[i], "WISE 3");
        }
        TIMER_STOP("WISE 3 (MASK)", &rsub);

        for (size_t i = 0; i < nonterms_count; i++) {
            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, matrices[i]));
            if (new_nnz != 0)
                t_empty_flags[i] = false;

            changed = changed || (nnzs[i] != new_nnz);
            nnzs[i] = new_nnz;
        }

#ifdef DEBUG_CFL_REACHBILITY
        GxB_Matrix_iso(&iso_flag, T[bin_rule.nonterm]);
        printf("[TERM1 TERM2] MULTIPLY, S: %d, A: %d, B: %d, "
               "I: %ld (ISO: %d)\n",
               bin_rule.nonterm, bin_rule.prod_A, bin_rule.prod_B, i, iso_flag);
#endif
    }

    printf("MXM1: %.3f, wise1: %.3f, MXM2: %.3f, wise2: %.3f, rsub: %.3f", mxm1, wise1,
           mxm2, wise2, rsub);

#ifdef DEBUG_CFL_REACHBILITY
    for (int32_t i = 0; i < nonterms_count; i++) {
        printf("MATRIX WITH INDEX %d:\n", i);
        GxB_print(T[i], GxB_SUMMARY);
    }
#endif

    for (int32_t i = 0; i < nonterms_count; i++) {
        outputs[i] = matrices[i];
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
