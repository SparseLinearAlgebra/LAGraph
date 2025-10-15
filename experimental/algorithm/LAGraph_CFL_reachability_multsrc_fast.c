//------------------------------------------------------------------------------
// LAGraph_CFL_reachability_multsrc.c: Optimized Multiple-Source Context-Free
// Language Reachability Matrix-Based Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Ilhom Kombaev, Vlasenco Daniel, Semyon Grigoriev, St. Petersburg State University.

//------------------------------------------------------------------------------

// Code is based on the "multiple-source CFPQ", described in the following paper:
//  * Arseny Terekhov et al., "Multiple-Source Context-Free Path Querying in Terms of Linear Algebra"
//  * URL: https://openproceedings.org/2021/conf/edbt/p48.pdf

// #define DEBUG_CFL_REACHABILITY

#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        LAGraph_Free((void **) &nnzs_T, msg);                                            \
        LAGraph_Free((void **) &nnzs_TSrc_B, msg);                                       \
        LAGraph_Free((void **) &nnzs_TSrc_C, msg);                                       \
        LAGraph_Free((void **) &ones_vec, msg);                                          \
        LAGraph_Free((void **) &T, msg);                                                 \
        LAGraph_Free((void **) &TSrc, msg);                                              \
        LAGraph_Free((void **) &MSrc, msg);                                              \
        LAGraph_Free((void **) &identity_matrix, msg);                                   \
        LAGraph_Free ((void **) &A, msg);                                                \
        LAGraph_Free ((void **) &a, msg);                                                \
        GrB_free(&true_scalar);                                                          \
    }

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        for (size_t i = 0; i < nonterms_count; i++) {                                    \
            GrB_free(&T[i].base);                                                             \
            GrB_free(&dT[i].base);                                                          \
            GrB_free(&TSrc[i].base);                                                          \
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

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                                                \
    {                                                                                   \
        rule.len_indices_str += snprintf(rule.indices_str + rule.len_indices_str,       \
        LAGRAPH_MSG_LEN - rule.len_indices_str,                                         \
                                         rule.count == 0 ? "%ld" : ", %ld", i);         \
        rule.count++;                                                                   \
    }

#define PRINT_MATRIX(_m) {                                                              \
    for (size_t _i = 0; _i < n; _i++) {                                                 \
        for (size_t _j = 0; _j < n; _j++) {                                             \
            if (GxB_Matrix_isStoredElement(_m, _i, _j) == GrB_SUCCESS) {                \
                printf("1 ");                                                           \
            }                                                                           \
            else {                                                                      \
                printf("0 ");                                                           \
            }                                                                           \
        }                                                                               \
        printf("\n");                                                                   \
    }                                                                                   \
}

#define PRINT_VECTOR(_v) {                                                              \
    for (size_t _i = 0; _i < n; _i++) {                                                 \
        if (GxB_Vector_isStoredElement(_v, _i) == GrB_SUCCESS) {                        \
            printf("1 ");                                                               \
        }                                                                               \
        else {                                                                          \
            printf("0 ");                                                               \
        }                                                                               \
    }                                                                                   \
    printf("\n");                                                                       \
}

#define PRINT_RULE(_r) {                                                                \
    printf("%c -> ", _r.nonterm + 'A' - 1);                                             \
    if (_r.prod_A == -1 && _r.prod_B == -1) {                                           \
        printf("eps\n");                                                                \
    } else if (_r.prod_A != -1 && _r.prod_B == -1) {                                    \
        printf("%c\n", _r.prod_A + 'a');                                                \
    } else {                                                                            \
        printf("%c %c\n", _r.prod_A + 'A'- 1 ,_r.prod_B + 'A' - 1);                     \
    }                                                                                   \
}

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

// LAGraph_CFL_reachability_multsrc: Optimized Multiple-Source Context-Free
// Language Reachability Matrix-Based Algorithm
//
// This function determines the set of vertex pairs (u, v) in a graph (represented by
// adjacency matrices) such that u is in the given set of source vertices `src`, and
// there is a path from u to v, where the edge labels form a word from the language
// generated by the context-free grammar (represented by `rules`).

GrB_Info LAGraph_CFL_reachability_multsrc_fast
(
    // Output
    GrB_Matrix *output, // A handle for the matrix containing results.
                        //
                        // output: (i, j) = true if and only if there is a path
                        // from node i to node j whose edge labels form a word
                        // derivable from the specified CFG.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, terms_count - 1]).
    GrB_Index *src,                 // Array of source vertices
    // These parameters will become unsigned
    int32_t src_count,              // The total number of source vertices
    int32_t terms_count,            // The total number of terminal symbols in the CFG.
    int32_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    size_t rules_count,             // The total number of rules in the CFG.
    char *msg,                      // Message string for error reporting.
    int8_t opt_mask                 // Optimizations mask
)
{
    // Declare workspace and clear the msg string, if not NULL
    CFL_Matrix *T;
    CFL_Matrix *dT;
    CFL_Matrix *TSrc;
    CFL_Matrix *Adj;
    CFL_Matrix MSrc;
    CFL_Matrix M1;
    CFL_Matrix M2;
    CFL_Matrix Temp1;
    CFL_Matrix Temp2;
    CFL_Matrix A;
    GrB_Vector a;
    GrB_Index n; // number of vertices in the graph
    CFL_Matrix iden;
    GrB_Matrix identity_matrix = NULL;
    GrB_Index *nnzs_T = NULL;
    GrB_Index *nnzs_TSrc_B = NULL;
    GrB_Index *nnzs_TSrc_C = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    bool iso_flag = false;
    GrB_Scalar true_scalar;
    GrB_Vector ones_vec;

    // Will change the interface and omit this check in the future
    if (nonterms_count < 0 || terms_count < 0)
        return GrB_INVALID_VALUE;

    if (!nonterms_count || !rules_count)
        return GrB_INVALID_VALUE;

    bool t_is_empty[nonterms_count]; // t_is_empty[i] == true <=> T[i] is empty
    bool t_src_is_empty[nonterms_count]; // t_src_is_empty[i] == true <=> TSrc[i] is empty

    if (!output || !rules || !adj_matrices || !src)
        return GrB_NULL_POINTER;
    if (src_count <= 0)
        return GrB_INVALID_VALUE;

    // Find null adjacency matrices
    bool found_null = false;
    for (int32_t i = 0; i < terms_count; i++) {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null) {
            ADD_TO_MSG("Adjacency matrices with these indices are null: ");
            ADD_TO_MSG("%d", i);
        } else {
            ADD_TO_MSG(", %d", i);
        }

        found_null = true;
    }
    if (found_null) {
        return GrB_NULL_POINTER;
    }

    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));
    if (n < src_count) return GrB_INVALID_VALUE;

    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);

    LG_TRY(LAGraph_Calloc((void **) &T, nonterms_count, sizeof(CFL_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **) &dT, nonterms_count, sizeof(CFL_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **) &TSrc, nonterms_count, sizeof(CFL_Matrix), msg));

    LG_TRY(LAGraph_Calloc((void **) &Adj, terms_count, sizeof(CFL_Matrix), msg));

    GRB_TRY(GrB_Vector_new(&ones_vec, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(ones_vec, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&identity_matrix, ones_vec, 0));

    LG_TRY(LAGraph_Calloc((void **) &nnzs_T, nonterms_count, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Calloc((void **) &nnzs_TSrc_B, nonterms_count, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Calloc((void **) &nnzs_TSrc_C, nonterms_count, sizeof(GrB_Index), msg));

    for (int32_t i = 0; i < terms_count; i++) {
        Adj[i] = CFL_matrix_from_base(adj_matrices[i]);
    }

    // Create nonterms matrices
    for (int32_t i = 0; i < nonterms_count; i++) {
        GrB_Matrix matrix;

        // dT[i] = CFL_matrix_from_base(adj_matrices[i]);
        dT[i] = CFL_matrix_create(n, n);

        GRB_TRY(GrB_Matrix_new(&matrix, GrB_BOOL, n, n));
        T[i] = ((opt_mask & OPT_LAZY) || (opt_mask & OPT_BLOCK))
                          ? CFL_matrix_from_base_lazy(matrix)
                          : CFL_matrix_from_base(matrix);

        TSrc[i] = CFL_matrix_create(n, n);
    }

    for (int32_t i = 0; i < src_count; i++) {
        GrB_Matrix_setElement(TSrc[0].base, true, src[i], src[i]);
    }

    t_src_is_empty[0] = false;

    GRB_TRY(GrB_Matrix_dup(&MSrc.base, TSrc[0].base));
    MSrc = CFL_matrix_from_base(MSrc.base);
    GRB_TRY(GrB_Vector_new(&a, GrB_BOOL, n));

    M1 = CFL_matrix_create(n, n);
    M2 = CFL_matrix_create(n, n);
    Temp1 = CFL_matrix_create(n, n);
    Temp2 = CFL_matrix_create(n, n);
    A = CFL_matrix_create(n, n);

    // #ifdef DEBUG_CFL_REACHABILITY
    // printf("MSrc:\n");
    // PRINT_MATRIX(MSrc);
    // #endif

    // Arrays for processing rules
    size_t eps_rules[rules_count], eps_rules_count = 0;   // [Variable -> eps]
    size_t term_rules[rules_count], term_rules_count = 0; // [Variable -> term]
    size_t bin_rules[rules_count], bin_rules_count = 0;   // [Variable -> AB]

    // Process rules
    typedef struct {
        size_t count;
        size_t len_indices_str;
        char indices_str[LAGRAPH_MSG_LEN];
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
            ADD_TO_MSG("indices of invalid rules: %s\n", nonterm_err.indices_str)
        }
        if (term_err.count > 0) {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("indices of invalid rules: %s\n", term_err.indices_str)
        }
        if (invalid_err.count > 0) {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("indices of invalid rules: %.120s\n", invalid_err.indices_str)
        }

        LG_FREE_ALL;
        return GrB_INVALID_VALUE;
    }

    // Rule [Variable -> term]
    for (size_t i = 0; i < term_rules_count; i++) {
        LAGraph_rule_WCNF term_rule = rules[term_rules[i]];

        // GxB_eWiseUnion(
        //     T[term_rule.nonterm].base, GrB_NULL, GrB_NULL, GxB_PAIR_BOOL,
        //     T[term_rule.nonterm].base, true_scalar, adj_matrices[term_rule.prod_A], true_scalar, GrB_NULL
        // );

        // t_is_empty[term_rule.nonterm] = false;

        GRB_TRY(CFL_wise(&T[term_rule.nonterm], &T[term_rule.nonterm], &Adj[term_rule.prod_A], true, opt_mask));

        // #ifdef DEBUG_CFL_REACHABILITY
        // GxB_Matrix_iso(&iso_flag, T[term_rule.nonterm]);
        // printf("[TERM] eWiseUnion: NONTERM: %d (ISO: %d)\n", term_rule.nonterm, iso_flag);
        // #endif
    }

    iden = CFL_matrix_from_base(identity_matrix);

    // Rule [Variable -> eps]
    for (size_t i = 0; i < eps_rules_count; i++) {
        LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];

        // GxB_eWiseUnion (
        //     T[eps_rule.nonterm].base,GrB_NULL,GxB_PAIR_BOOL,GxB_PAIR_BOOL,
        //     T[eps_rule.nonterm].base,true_scalar,identity_matrix,true_scalar,GrB_NULL
        // );

        // dt?
        GRB_TRY(CFL_wise(&T[eps_rule.nonterm], &T[eps_rule.nonterm], &iden, true, opt_mask));
        
        // t_is_empty[eps_rule.nonterm] = false;

        // #ifdef DEBUG_CFL_REACHABILITY
        // GxB_Matrix_iso(&iso_flag, T[eps_rule.nonterm]);
        // printf("[EPS] eWiseUnion: NONTERM: %d (ISO: %d)\n",
        //         eps_rule.nonterm, iso_flag);
        // #endif
    }

    // needed?
    for (int32_t i = 0; i < nonterms_count; i++) {
        GRB_TRY(GrB_Matrix_dup(&dT[i].base, T[i].base));
    }

    // Rule [Variable -> Variable1 Variable2]
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];
            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("Rule: ");
            // PRINT_RULE(bin_rule)
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("Before M = TSrc^A * T^B:\n");
            // printf("TSrc^A:\n");
            // PRINT_MATRIX(TSrc[bin_rule.nonterm]);
            // printf("T^B:\n");
            // PRINT_MATRIX(T[bin_rule.prod_A]);
            // #endif


            GRB_TRY(CFL_mxm(&M1, &TSrc[bin_rule.nonterm], &T[bin_rule.prod_A], false, false, opt_mask));


            GRB_TRY(CFL_wise(&T[bin_rule.nonterm], &T[bin_rule.nonterm], &dT[bin_rule.nonterm], false, opt_mask));


            // ? <- this means that without this step, unit tests pass,
            // although algorithm requires it.
            GRB_TRY(CFL_mxm(&M2, &TSrc[bin_rule.nonterm], &dT[bin_rule.prod_A], false, false, opt_mask));

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("After M = TSrc^A * T^B:\n");
            // printf("M:\n");
            // PRINT_MATRIX(M);
            // printf("Before T^A = T^A + M * T^C:\n");
            // printf("T^A:\n");
            // PRINT_MATRIX(T[bin_rule.nonterm]);
            // printf("T^C:\n");
            // PRINT_MATRIX(T[bin_rule.prod_B]);
            // #endif

            // GrB_BinaryOp acc_op = t_is_empty[bin_rule.nonterm] ? GrB_NULL : GxB_ANY_BOOL;


            GRB_TRY(CFL_mxm(&Temp1, &M1, &dT[bin_rule.prod_B], false, false, opt_mask));

            // GRB_TRY(GrB_mxm(T[bin_rule.nonterm], GrB_NULL, acc_op, GxB_ANY_PAIR_BOOL,
            //             M, T[bin_rule.prod_B], GrB_NULL));


            // ?
            GRB_TRY(CFL_mxm(&Temp2, &M2, &T[bin_rule.prod_B], false, false, opt_mask));


            GRB_TRY(CFL_wise(&dT[bin_rule.nonterm], &Temp1, &Temp2, false, opt_mask));

            // GRB_TRY(GrB_eWiseAdd(dT[bin_rule.nonterm].base, GrB_NULL, GrB_NULL, GxB_ANY_BOOL,
            //             Temp1.base, Temp2.base, GrB_NULL));


            // ?
            GRB_TRY(CFL_rsub(&dT[bin_rule.nonterm], &T[bin_rule.nonterm], opt_mask));


            // GrB_eWiseAdd(dT[bin_rule.nonterm].base, GrB_NULL, GrB_NULL, GrB_PLUS_BOOL,
            //             dT[bin_rule.nonterm].base, T[bin_rule.nonterm].base, GrB_NULL);

            // CFL_matrix_update(&dT[bin_rule.nonterm]);

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("After T^A = T^A + M * T^C:\n");
            // printf("T^A:\n");
            // PRINT_MATRIX(T[bin_rule.nonterm]);
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("Before TSrc^B = TSrc^B + TSrc^A:\n");
            // printf("TSrc^A:\n");
            // PRINT_MATRIX(TSrc[bin_rule.nonterm]);
            // printf("TSrc^B:\n");
            // PRINT_MATRIX(TSrc[bin_rule.prod_A]);
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("After TSrc^B = TSrc^B + TSrc^A:\n");
            // printf("TSrc^B:\n");
            // PRINT_MATRIX(TSrc[bin_rule.prod_A]);
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("Before A = dest(M)\n");
            // printf("M:\n");
            // PRINT_MATRIX(M);
            // #endif


            // Update source vertices matrix to find appropriate paths only
            // M1[i, j] == 1 => A[j, j] == 1
            GRB_TRY(GrB_vxm(a, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, ones_vec, M1.base, GrB_NULL));
            GRB_TRY(GrB_Matrix_free(&A.base));
            GRB_TRY(GrB_Matrix_diag(&A.base, a, 0));
            A = CFL_matrix_from_base(A.base);


            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("After A = dest(M)\n");
            // printf("a:\n");
            // PRINT_VECTOR(a);
            // printf("A:\n");
            // PRINT_MATRIX(A);
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("Before TSrc^C = TSrc^c + A\n");
            // printf("TSrc^C:\n");
            // PRINT_MATRIX(TSrc[bin_rule.prod_B])
            // printf("A:\n");
            // PRINT_MATRIX(A)
            // #endif


            GRB_TRY(CFL_wise(&TSrc[bin_rule.prod_A], &TSrc[bin_rule.prod_A], &TSrc[bin_rule.nonterm], false, opt_mask));

            // GRB_TRY(GrB_eWiseAdd(TSrc[bin_rule.prod_A].base, GrB_NULL, GrB_NULL, GxB_ANY_BOOL,
            //             TSrc[bin_rule.prod_A].base, TSrc[bin_rule.nonterm].base, GrB_NULL));


            GRB_TRY(CFL_wise(&TSrc[bin_rule.prod_B], &TSrc[bin_rule.prod_B], &A, false, opt_mask));

            // GRB_TRY(GrB_eWiseAdd(TSrc[bin_rule.prod_B].base, GrB_NULL, GrB_NULL, GxB_ANY_BOOL,
            //     TSrc[bin_rule.prod_B].base, A.base, GrB_NULL));



            // #ifdef DEBUG_CFL_REACHABILITY
            // printf("After TSrc^C = TSrc^c + A\n");
            // printf("TSrc^C:\n");
            // PRINT_MATRIX(TSrc[bin_rule.prod_B])
            // #endif

            // #ifdef DEBUG_CFL_REACHABILITY
            // GxB_Matrix_iso(&iso_flag, T[bin_rule.nonterm]);
            // printf("[TERM1 TERM2] MULTIPLY, S: %d, A: %d, B: %d, "
            //        "I: %ld (ISO: %d)\n",
            //        bin_rule.nonterm, bin_rule.prod_A, bin_rule.prod_B, i, iso_flag);
            // #endif

            // Check if any of the matrices changed. If not, job is done.
            GrB_Index nnz_T, nnz_TSrc_B, nnz_TSrc_C;

            GRB_TRY(GrB_Matrix_nvals(&nnz_T, T[bin_rule.nonterm].base));
            GRB_TRY(GrB_Matrix_nvals(&nnz_TSrc_B, TSrc[bin_rule.prod_A].base));
            GRB_TRY(GrB_Matrix_nvals(&nnz_TSrc_C, TSrc[bin_rule.prod_B].base));
            
            if (nnz_T != 0) t_is_empty[bin_rule.nonterm] = false;
            if (nnz_TSrc_B != 0) t_src_is_empty[bin_rule.prod_A] = false;
            if (nnz_TSrc_C != 0) t_src_is_empty[bin_rule.prod_B] = false;

            changed = changed || (nnzs_T[bin_rule.nonterm] != nnz_T);
            changed = changed || (nnzs_TSrc_B[bin_rule.prod_A] != nnz_TSrc_B);
            changed = changed || (nnzs_TSrc_C[bin_rule.prod_B] != nnz_TSrc_C);

            nnzs_T[bin_rule.nonterm] = nnz_T;
            nnzs_TSrc_B[bin_rule.prod_A] = nnz_TSrc_B;
            nnzs_TSrc_C[bin_rule.prod_B] = nnz_TSrc_C;
        }
    }

    // #ifdef DEBUG_CFL_REACHABILITY
    // printf("Before MSrc = MSrc * T^S\n");
    // printf("MSrc:\n");
    // PRINT_MATRIX(MSrc)
    // printf("T^S:\n");
    // PRINT_MATRIX(T[0]);
    // #endif

    GRB_TRY(CFL_mxm(&MSrc, &MSrc, &T[0], false, false, opt_mask));
    // GRB_TRY(GrB_mxm(MSrc.base, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL,
    //         MSrc.base, T[0].base, GrB_NULL));

    // #ifdef DEBUG_CFL_REACHABILITY
    // printf("After MSrc = MSrc * T^S\n");
    // printf("MSrc, output:\n");
    // PRINT_MATRIX(MSrc)
    // #endif

    GRB_TRY(GrB_Matrix_dup(output, MSrc.base));

    LG_FREE_ALL;
    return GrB_SUCCESS;
}
