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

// clang-format off
#if BENCH_CFL_REACHBILITY
    #define IS_ISO(matrix, str)                                                              \
    {                                                                                    \
        bool iso_flag;                                                                   \
        GrB_Index nnz;                                                                   \
        GxB_Matrix_iso(&iso_flag, matrix);                                               \
        GrB_Matrix_nvals(&nnz, matrix);                                                  \
        if (!iso_flag && nnz) {                                                          \
            printf("-----ISO ALERT----- (%s)\n", str);                                   \
            GxB_print(matrix, 1);                                                        \
            printf("-------------------\n");                                             \
        }                                                                                \
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

    #define IS_ROW(matrix, str)                                                              \
    {                                                                                    \
        int32_t orientation;                                                             \
        GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT);                     \
        if (orientation != GrB_ROWMAJOR) {                                               \
            printf("-----NOT A ROW----- (%s)\n", str);                                   \
            GxB_print(matrix, 1);                                                        \
            printf("-------------------\n");                                             \
        }                                                                                \
    }

    #define IS_COL(matrix, str)                                                              \
    {                                                                                    \
        int32_t orientation;                                                             \
        GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT);                     \
        if (orientation != GrB_COLMAJOR) {                                               \
            printf("-----NOT A COL----- (%s)\n", str);                                   \
            GxB_print(matrix, 1);                                                        \
            printf("-------------------\n");                                             \
        }                                                                                \
    }
#else
    #define IS_ISO(matrix, str)
    #define TIMER_START()
    #define TIMER_STOP(label, accumulator)
    #define IS_ROW(matrix, str)
    #define IS_COL(matrix, str)
#endif
// clang-format on

#define SKIP_IF_NULL(matrix)                                                             \
    GrB_Matrix_nvals(&new_nnz, matrix);                                                  \
    if (new_nnz == 0) {                                                                  \
        continue;                                                                        \
    }

#define TO_COL(matrix) GrB_set(matrix, GrB_COLMAJOR, GrB_STORAGE_ORIENTATION_HINT)
#define TO_ROW(matrix) GrB_set(matrix, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)

#define TRY(GrB_method)                                                                  \
    {                                                                                    \
        GrB_Info LG_GrB_Info = GrB_method;                                               \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                 \
            return LG_GrB_Info;                                                          \
        }                                                                                \
    }

typedef struct {
    GrB_Matrix base;
    GrB_Matrix base_row;
    GrB_Matrix base_col;
    GrB_Index nvals;
    GrB_Index size;
    int32_t format;
    bool is_both;
} Matrix;

void matrix_update(Matrix *matrix) {
    GrB_Matrix_nvals(&matrix->nvals, matrix->base);
    GrB_Matrix_nrows(&matrix->size, matrix->base);
    GrB_get(matrix->base, &matrix->format, GrB_STORAGE_ORIENTATION_HINT);
}

Matrix matrix_from_base(GrB_Matrix matrix) {
    Matrix result;
    result.base = matrix;
    result.base_row = matrix;
    result.base_col = NULL;
    result.nvals = 0;
    result.size = 0;
    result.format = GrB_ROWMAJOR;
    result.is_both = false;
    matrix_update(&result);
    return result;
}

void matrix_format_sync(Matrix *matrix) {
    if (!matrix->is_both) {
        // if (matrix->format == GrB_ROWMAJOR) {
        //     GrB_Matrix_assign(matrix->base_row, GrB_NULL, GrB_NULL, matrix->base,
        //     GrB_ALL,
        //                       matrix->size, GrB_ALL, matrix->size, GrB_NULL);
        // } else {
        //     GrB_Matrix_assign(matrix->base_col, GrB_NULL, GrB_NULL, matrix->base,
        //     GrB_ALL,
        //                       matrix->size, GrB_ALL, matrix->size, GrB_NULL);
        // }

        return;
    }

    // GrB_Matrix_assign(matrix->base_row, GrB_NULL, GrB_NULL, matrix->base, GrB_ALL,
    //                   matrix->size, GrB_ALL, matrix->size, GrB_NULL);
    // GrB_Matrix_assign(matrix->base_col, GrB_NULL, GrB_NULL, matrix->base, GrB_ALL,
    //                   matrix->size, GrB_ALL, matrix->size, GrB_NULL);

    GrB_Matrix *new_matrix =
        matrix->format == GrB_ROWMAJOR ? &matrix->base_col : &matrix->base_row;
    GrB_Matrix *old_matrix =
        matrix->format == GrB_ROWMAJOR ? &matrix->base_row : &matrix->base_col;

    GrB_Matrix_assign(*new_matrix, GrB_NULL, GrB_NULL, *old_matrix, GrB_ALL, matrix->size,
                      GrB_ALL, matrix->size, GrB_NULL);

    // is_matrix_equal(matrix->base_row, matrix->base_col);
    // is_matrix_equal(matrix->base_row, matrix->base);
}

void matrix_to_format(Matrix *matrix, int32_t format, bool is_both) {
    // Matrix contain both formats so just switch base matrix
    if (matrix->is_both) {
        matrix->base = format == GrB_ROWMAJOR ? matrix->base_row : matrix->base_col;
        matrix->format = format;
        return;
    }

    // No changes required
    if (matrix->format == format) {
        return;
    }

    // Matrix contain just one matrix and format is not same
    GrB_Matrix *new_matrix =
        matrix->format == GrB_ROWMAJOR ? &matrix->base_col : &matrix->base_row;
    GrB_Matrix *old_matrix =
        matrix->format == GrB_ROWMAJOR ? &matrix->base_row : &matrix->base_col;

    if (is_both) {
        GrB_Matrix_new(new_matrix, GrB_BOOL, matrix->size, matrix->size);
        GrB_Matrix_assign(*new_matrix, GrB_NULL, GrB_NULL, *old_matrix, GrB_ALL,
                          matrix->size, GrB_ALL, matrix->size, GrB_NULL);
        matrix->is_both = true;
    } else {
        *new_matrix = *old_matrix;
        *old_matrix = NULL;
    }

    matrix->base = format == GrB_ROWMAJOR ? matrix->base_row : matrix->base_col;
    format == GrB_ROWMAJOR ? TO_ROW(matrix->base) : TO_COL(matrix->base);
    matrix->format = format;
    return;
}

GrB_Info matrix_clear(Matrix *A) { return GrB_Matrix_clear(A->base); }

GrB_Info matrix_clear_format(Matrix *A) {
    if (!A->is_both) {
        matrix_clear(A);
    }

    matrix_to_format(A, GrB_ROWMAJOR, false);
    GrB_Info result = matrix_clear(A);

    if (result < GrB_SUCCESS) {
        return result;
    }

    matrix_to_format(A, GrB_COLMAJOR, false);
    return matrix_clear(A);
}

GrB_Info matrix_clear_empty(Matrix *A) {
    if (A->nvals == 0) {
        return GrB_SUCCESS;
    }

    return matrix_clear_format(A);
}

GrB_Info matrix_dup(Matrix *output, Matrix *input) {
    return GrB_Matrix_assign(output->base, GrB_NULL, GrB_NULL, input->base, GrB_ALL,
                             input->size, GrB_ALL, input->size, GrB_NULL);
}

GrB_Info matrix_dup_format(Matrix *output, Matrix *input) {
    if (!output->is_both) {
        return matrix_dup(output, input);
    }

    matrix_to_format(output, GrB_ROWMAJOR, false);
    GrB_Info result = matrix_dup(output, input);

    if (result < GrB_SUCCESS) {
        return result;
    }

    matrix_to_format(output, GrB_COLMAJOR, false);
    return matrix_dup(output, input);
}

GrB_Info matrix_dup_empty(Matrix *output, Matrix *input) {
    if (input->nvals == 0) {
        return matrix_clear_empty(output);
    }

    return matrix_dup_format(output, input);
}

GrB_Info matrix_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    GrB_Info result = GrB_mxm(output->base, GrB_NULL, accum ? GxB_ANY_BOOL : GrB_NULL,
                              GxB_ANY_PAIR_BOOL, first->base, second->base, GrB_NULL);
    IS_ISO(output->base, "MXM output");
    return result;
}

GrB_Info matrix_mxm_format(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    int32_t desired_orientation =
        first->nvals > second->nvals ? GrB_COLMAJOR : GrB_ROWMAJOR;

    if (!first->is_both && first->format != desired_orientation &&
        !(first->nvals > second->nvals / 3.0)) {
        GrB_Info result = matrix_mxm(output, first, second, accum);
        return result;
    }

    matrix_to_format(first, desired_orientation, true);
    matrix_to_format(second, desired_orientation, false);
    matrix_to_format(output, desired_orientation, false);
    GrB_Info result = matrix_mxm(output, first, second, accum);
    return result;
}

GrB_Info matrix_mxm_empty(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (first->nvals == 0 || second->nvals == 0)
        return GrB_SUCCESS;

    return matrix_mxm_format(output, first, second, accum);
}

GrB_Info matrix_rmxm_format(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    Matrix *temp = first;
    first = second;
    second = temp;

    int32_t desired_orientation =
        first->nvals > second->nvals ? GrB_ROWMAJOR : GrB_COLMAJOR;

    if (!first->is_both && first->format != desired_orientation &&
        !(first->nvals > second->nvals / 3.0)) {
        GrB_Info result = matrix_mxm(output, second, first, accum);
        return result;
    }

    matrix_to_format(first, desired_orientation, true);
    matrix_to_format(second, desired_orientation, false);
    matrix_to_format(output, desired_orientation, false);
    GrB_Info result = matrix_mxm(output, second, first, accum);
    return result;
}

GrB_Info matrix_rmxm_empty(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (first->nvals == 0 || second->nvals == 0)
        return GrB_SUCCESS;

    return matrix_rmxm_format(output, first, second, accum);
}

GrB_Info matrix_wise(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    GrB_BinaryOp accum_op = accum ? GxB_ANY_BOOL : GrB_NULL;

    return GrB_eWiseAdd(output->base, GrB_NULL, accum_op, GxB_ANY_BOOL, first->base,
                        second->base, GrB_NULL);
}

GrB_Info matrix_wise_format(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (!output->is_both) {
        return matrix_wise(output, first, second, accum);
    }

    matrix_to_format(output, GrB_ROWMAJOR, false);
    GrB_Info result = matrix_wise(output, first, second, accum);

    if (result < GrB_SUCCESS) {
        return result;
    }

    matrix_to_format(output, GrB_COLMAJOR, false);
    return matrix_wise(output, first, second, accum);
}

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (first->nvals == 0 && second->nvals == 0) {
        if (accum) {
            return GrB_SUCCESS;
        }

        return matrix_clear_format(output);
    }

    if (first->nvals == 0) {
        if (accum) {
            return matrix_wise_format(output, output, second, false);
        }

        return matrix_dup_format(output, second);
    }

    if (second->nvals == 0) {
        if (accum) {
            return matrix_wise_format(output, output, first, false);
        }

        return matrix_dup_format(output, first);
    }

    return matrix_wise_format(output, first, second, accum);
}

GrB_Info matrix_rsub(Matrix *output, Matrix *mask) {
    return GrB_eWiseAdd(output->base, mask->base, GrB_NULL, GxB_ANY_BOOL, output->base,
                        output->base, GrB_DESC_RSC);
}

GrB_Info matrix_rsub_format(Matrix *output, Matrix *mask) {
    Matrix *larger_matrix = output->nvals > mask->nvals ? output : mask;
    matrix_to_format(output, larger_matrix->format, false);
    matrix_to_format(mask, larger_matrix->format, false);

    if (!output->is_both) {
        return matrix_rsub(output, mask);
    }

    printf("LOOOOOO\n\n");

    matrix_to_format(output, GrB_ROWMAJOR, false);
    GrB_Info result = matrix_rsub(output, mask);

    if (result < GrB_SUCCESS) {
        return result;
    }

    matrix_to_format(output, GrB_COLMAJOR, false);
    return matrix_rsub(output, mask);
}

GrB_Info matrix_rsub_empty(Matrix *output, Matrix *mask) {
    if (mask->nvals == 0 || output->nvals == 0) {
        return GrB_SUCCESS;
    }

    return matrix_rsub_format(output, mask);
}

// LAGraph_CFL_reachability: Context-Free Language Reachability Matrix-Based Algorithm
//
// This function determines the set of vertex pairs (u, v) in a graph (represented by
// adjacency matrices) such that there is a path from u to v, where the edge labels
// form a word from the language generated by the context-free grammar (represented by
// `rules`).
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
// There are paths from node [1] to node [3] and from node [1] to node [2] that form
// the word "ab" ([1]-a->[2]-b->[3] and [1]-a->[5]-b->[2]). The word "ab" is in the
// language generated by our context-free grammar, so the pairs (1, 3) and (1, 2) will
// be included in the result.
//
// Note: It doesn't matter how many paths exist from node [A] to node [B] that form a
// word in the language. If at least one path exists, the pair ([A], [B]) will be
// included in the result.
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
    Matrix *delta_matrices;
    Matrix *matrices;
    Matrix *temp_matrices;
    GrB_Matrix identity_matrix = NULL;
    uint64_t *nnzs = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    GrB_Index *indexes = NULL;

    GrB_Scalar true_scalar;
    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);

    LG_TRY(LAGraph_Calloc((void **)&T, nonterms_count, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&delta_matrices, nonterms_count, sizeof(Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&matrices, nonterms_count, sizeof(Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&temp_matrices, nonterms_count, sizeof(Matrix), msg));

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
        GrB_Matrix matrix;

        GRB_TRY(GrB_Matrix_new(&T[i], GrB_BOOL, n, n));

        GRB_TRY(GrB_Matrix_new(&matrix, GrB_BOOL, n, n));
        delta_matrices[i] = matrix_from_base(matrix);

        GRB_TRY(GrB_Matrix_new(&matrix, GrB_BOOL, n, n));
        matrices[i] = matrix_from_base(matrix);

        GRB_TRY(GrB_Matrix_new(&matrix, GrB_BOOL, n, n));
        temp_matrices[i] = matrix_from_base(matrix);
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

        GxB_eWiseUnion(delta_matrices[term_rule.nonterm].base, GrB_NULL, GrB_NULL,
                       GxB_PAIR_BOOL, delta_matrices[term_rule.nonterm].base, true_scalar,
                       adj_matrices[term_rule.prod_A], true_scalar, GrB_NULL);
        matrix_update(&delta_matrices[term_rule.nonterm]);

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

        GxB_eWiseUnion(delta_matrices[eps_rule.nonterm].base, GrB_NULL, GxB_PAIR_BOOL,
                       GxB_PAIR_BOOL, delta_matrices[eps_rule.nonterm].base, true_scalar,
                       identity_matrix, true_scalar, GrB_NULL);
        matrix_update(&delta_matrices[eps_rule.nonterm]);

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
    while (changed) {
        iteration++;
        changed = false;

#if BENCH_CFL_REACHBILITY
        printf("\n--- ITERATARION %ld ---\n", iteration);
#endif

        for (int32_t i = 0; i < nonterms_count; i++) {
            GRB_TRY(matrix_clear_empty(&temp_matrices[i]));
            matrix_update(&temp_matrices[i]);
        }

        TIMER_START();
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            matrix_mxm_empty(&temp_matrices[bin_rule.nonterm], &matrices[bin_rule.prod_A],
                             &delta_matrices[bin_rule.prod_B], false);
            matrix_update(&temp_matrices[bin_rule.nonterm]);
        }
        TIMER_STOP("MXM 1", &mxm1);

        TIMER_START()
        for (int32_t i = 0; i < nonterms_count; i++) {
            matrix_wise_empty(&matrices[i], &matrices[i], &delta_matrices[i], false);
            matrix_update(&matrices[i]);
        }
        TIMER_STOP("WISE 1", &wise1);

        TIMER_START()
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            matrix_rmxm_empty(&temp_matrices[bin_rule.nonterm],
                              &delta_matrices[bin_rule.prod_A],
                              &matrices[bin_rule.prod_B], true);
            matrix_update(&temp_matrices[bin_rule.nonterm]);
        }
        TIMER_STOP("MXM 2", &mxm2);

        TIMER_START();
        for (int32_t i = 0; i < nonterms_count; i++) {
            matrix_dup_empty(&delta_matrices[i], &temp_matrices[i]);
        }
        TIMER_STOP("WISE 2 (copy)", &wise2);

        TIMER_START();
        for (int32_t i = 0; i < nonterms_count; i++) {
            matrix_rsub_empty(&delta_matrices[i], &matrices[i]);
            matrix_update(&delta_matrices[i]);
        }
        TIMER_STOP("WISE 3 (MASK)", &rsub);

        for (int32_t i = 0; i < nonterms_count; i++) {
            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, matrices[i].base));

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

#if BENCH_CFL_REACHBILITY
    printf("MXM1: %.3f, wise1: %.3f, MXM2: %.3f, wise2: %.3f, rsub: %.3f", mxm1, wise1,
           mxm2, wise2, rsub);
#endif

#ifdef DEBUG_CFL_REACHBILITY
    for (int32_t i = 0; i < nonterms_count; i++) {
        printf("MATRIX WITH INDEX %d:\n", i);
        GxB_print(T[i], GxB_SUMMARY);
    }
#endif

    for (int32_t i = 0; i < nonterms_count; i++) {
        outputs[i] = matrices[i].base;
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
