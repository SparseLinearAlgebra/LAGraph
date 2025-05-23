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
        for (int32_t i = 0; i < nonterms_count; i++) {                                   \
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

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)

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

enum Matrix_block { CELL, VEC_HORIZ, VEC_VERT };

typedef struct Matrix {
    GrB_Matrix base;
    GrB_Matrix base_row;
    GrB_Matrix base_col;
    struct Matrix *base_matrices;
    size_t base_matrices_count;
    GrB_Index nvals;
    GrB_Index nrows;
    GrB_Index ncols;
    int32_t format;
    enum Matrix_block block_type;
    bool is_both;
    bool is_lazy;
} Matrix;

void matrix_update(Matrix *matrix) {
    if (matrix->base_matrices_count == 0) {
        GrB_Matrix_nvals(&matrix->nvals, matrix->base);
    } else {
        size_t new_nnz = 0;
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            new_nnz += matrix->base_matrices[i].nvals;
        }

        matrix->nvals = new_nnz;
    }
    GrB_Matrix_nrows(&matrix->nrows, matrix->base);
    GrB_Matrix_ncols(&matrix->ncols, matrix->base);

    if (matrix->nrows > matrix->ncols) {
        matrix->block_type = VEC_VERT;
    }

    if (matrix->ncols > matrix->nrows) {
        matrix->block_type = VEC_HORIZ;
    }
    // GrB_get(matrix->base, &matrix->format, GrB_STORAGE_ORIENTATION_HINT);
}

Matrix matrix_from_base(GrB_Matrix matrix) {
    Matrix result;
    result.base = matrix;
    result.base_row = matrix;
    result.base_col = NULL;
    result.base_matrices = malloc(sizeof(Matrix) * 40);
    result.base_matrices_count = 0;
    result.nvals = 0;
    result.nrows = 0;
    result.ncols = 0;
    result.block_type = CELL;
    result.format = GrB_ROWMAJOR;
    result.is_both = false;
    matrix_update(&result);
    return result;
}

Matrix matrix_create(GrB_Index nrows, GrB_Index ncols) {
    GrB_Matrix _result;
    GrB_Matrix_new(&_result, GrB_BOOL, nrows, ncols);
    Matrix result = matrix_from_base(_result);

    return result;
};

void matrix_free(Matrix *matrix) {
    free(matrix->base_matrices);
    GrB_free(&matrix->base);
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
        GrB_Matrix_new(new_matrix, GrB_BOOL, matrix->nrows, matrix->ncols);
        GrB_Matrix_assign(*new_matrix, GrB_NULL, GrB_NULL, *old_matrix, GrB_ALL,
                          matrix->nrows, GrB_ALL, matrix->ncols, GrB_NULL);
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

GrB_Info matrix_clear(Matrix *A) {
    GrB_Info result = GrB_Matrix_clear(A->base);
    matrix_update(A);
    return result;
}

GrB_Info matrix_clear_format(Matrix *A) {
    if (!A->is_both) {
        return matrix_clear(A);
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

void block_matrix_hyper_rotate_i(Matrix *matrix, enum Matrix_block format) {
    if (matrix->block_type == CELL) {
        return;
    }

    if (matrix->block_type == format) {
        return;
    }

    GrB_Scalar scalar_true;
    GrB_Scalar_new(&scalar_true, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(scalar_true, true);

    if (matrix->format == VEC_VERT) {
        GrB_Index *nrows = malloc(matrix->nvals * sizeof(GrB_Index));
        GrB_Index *ncols = malloc(matrix->nvals * sizeof(GrB_Index));

        GrB_Matrix_extractTuples_BOOL(nrows, ncols, NULL, &matrix->nvals, matrix->base);

        for (size_t i = 0; i < matrix->nvals; i++) {
            ncols[i] = ncols[i] + ncols[i] / matrix->ncols * matrix->ncols;
            nrows[i] = nrows[i] % matrix->ncols;
        }

        GxB_Matrix_build_Scalar(matrix->base, nrows, ncols, scalar_true, matrix->nvals);
        free(nrows);
        free(ncols);
        GrB_free(&scalar_true);
        return;
    }

    if (matrix->format == VEC_HORIZ) {
        GrB_Index *nrows = malloc(matrix->nvals * sizeof(GrB_Index));
        GrB_Index *ncols = malloc(matrix->nvals * sizeof(GrB_Index));

        GrB_Matrix_extractTuples_BOOL(nrows, ncols, NULL, &matrix->nvals, matrix->base);

        for (size_t i = 0; i < matrix->nvals; i++) {
            nrows[i] = nrows[i] + nrows[i] / matrix->nrows * matrix->nrows;
            ncols[i] = ncols[i] % matrix->ncols;
        }

        GxB_Matrix_build_Scalar(matrix->base, nrows, ncols, scalar_true, matrix->nvals);
        free(nrows);
        free(ncols);
        GrB_free(&scalar_true);
        return;
    }
}

void block_matrix_to_diag(Matrix *diag, Matrix *input) {
    if (input->block_type == CELL) {
        exit(-1);
    }

    GrB_Scalar scalar_true;
    GrB_Scalar_new(&scalar_true, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(scalar_true, true);

    GrB_Index *rows = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Index *cols = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &input->nvals, input->base);

    if (input->block_type == VEC_HORIZ) {
        for (size_t i = 0; i < input->nvals; i++) {
            rows[i] = rows[i] + cols[i] / input->nrows * input->nrows;
        }
    }

    if (input->block_type == VEC_VERT) {
        for (size_t i = 0; i < input->nvals; i++) {
            cols[i] = cols[i] + rows[i] / input->ncols * input->ncols;
        }
    }

    GxB_Matrix_build_Scalar(diag->base, rows, cols, scalar_true, input->nvals);

    free(rows);
    free(cols);
    GrB_free(&scalar_true);
}

GrB_Info matrix_dup_empty(Matrix *output, Matrix *input);

void block_matrix_reduce(Matrix *matrix, Matrix *input) {
    if (input->block_type == CELL) {
        matrix_dup_empty(matrix, input);
    }

    GrB_Scalar scalar_true;
    GrB_Scalar_new(&scalar_true, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(scalar_true, true);

    GrB_Index *rows = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Index *cols = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &input->nvals, input->base);

    if (input->format == VEC_VERT) {
        for (size_t i = 0; i < input->nvals; i++) {
            rows[i] = rows[i] % input->ncols;
        }
    }

    if (input->format == VEC_HORIZ) {
        for (size_t i = 0; i < input->nvals; i++) {
            cols[i] = cols[i] % input->nrows;
        }
    }

    GxB_Matrix_build_Scalar(matrix->base, rows, cols, scalar_true, input->nvals);

    free(rows);
    free(cols);
    GrB_free(&scalar_true);
}

void block_matrix_repeat_into_vector(Matrix *matrix, Matrix *input,
                                     GrB_Index block_count) {
    GrB_Matrix *tiles = malloc(block_count * sizeof(GrB_Matrix));
    for (size_t i = 0; i < block_count; i++) {
        tiles[i] = input->base;
    }

    GxB_Matrix_concat(matrix->base, tiles, block_count, 1, GrB_NULL);
}

GrB_Info matrix_dup(Matrix *output, Matrix *input) {
    GrB_Info result =
        GrB_Matrix_assign(output->base, GrB_NULL, GrB_NULL, input->base, GrB_ALL,
                          input->nrows, GrB_ALL, input->ncols, GrB_NULL);

    matrix_update(output);
    return result;
}

GrB_Info matrix_dup_format(Matrix *output, Matrix *input) {
    if (!output->is_both) {
        Matrix *larger = output->nvals > input->nvals ? output : input;

        matrix_to_format(output, larger->format, false);
        matrix_to_format(input, larger->format, false);

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

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum);

GrB_Info matrix_sort_lazy(Matrix *A, bool reverse) {
    for (size_t i = 0; i < A->base_matrices_count; i++) {
        for (size_t j = i + 1; j < A->base_matrices_count; j++) {
            Matrix first = reverse ? A->base_matrices[i] : A->base_matrices[j];
            Matrix second = reverse ? A->base_matrices[j] : A->base_matrices[i];
            if (first.nvals < second.nvals) {
                Matrix temp = A->base_matrices[i];
                A->base_matrices[i] = A->base_matrices[j];
                A->base_matrices[j] = temp;
            }
        }
    }
}

GrB_Info matrix_combine_lazy(Matrix *A, size_t threshold) {
    Matrix *new_matrices = malloc(sizeof(Matrix) * 50);
    size_t new_size = 0;

    matrix_sort_lazy(A, false);

    for (size_t i = 0; i < A->base_matrices_count; i++) {
        if (A->base_matrices[i].nvals <= threshold && new_size > 0) {
            matrix_wise_empty(&new_matrices[new_size - 1], &new_matrices[new_size - 1],
                              &A->base_matrices[i], false);
            GrB_free(&A->base_matrices[i].base);
        } else {
            new_matrices[new_size++] = A->base_matrices[i];
        }
    }

    A->base_matrices = new_matrices;
    A->base_matrices_count = new_size;

    return GrB_SUCCESS;
}

GrB_Info matrix_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum,
                    bool swap) {
    Matrix *left = swap ? second : first;
    Matrix *right = swap ? first : second;

    GrB_Info result = GrB_mxm(output->base, GrB_NULL, accum ? GxB_ANY_BOOL : GrB_NULL,
                              GxB_ANY_PAIR_BOOL, left->base, right->base, GrB_NULL);
    IS_ISO(output->base, "MXM output");
    matrix_update(output);
    return result;
}

GrB_Info matrix_mxm_format(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           bool swap) {
    GrB_Index left_nvals = swap ? second->nvals : first->nvals;
    GrB_Index right_nvals = swap ? first->nvals : second->nvals;

    int32_t desired_orientation = left_nvals <= right_nvals ? GrB_ROWMAJOR : GrB_COLMAJOR;

    if (!first->is_both && first->format != desired_orientation &&
        !(first->nvals > second->nvals / 3.0)) {
        GrB_Info result = matrix_mxm(output, first, second, accum, swap);
        return result;
    }

    matrix_to_format(first, desired_orientation, true);
    matrix_to_format(second, desired_orientation, false);
    matrix_to_format(output, desired_orientation, false);
    GrB_Info result = matrix_mxm(output, first, second, accum, swap);
    return result;
}

GrB_Info matrix_mxm_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          bool swap) {
    if (first->nvals == 0 || second->nvals == 0) {
        if (accum) {
            return GrB_SUCCESS;
        }

        matrix_clear_empty(output);
    }

    return matrix_mxm_format(output, first, second, accum, swap);
}

GrB_Info matrix_mxm_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum,
                         bool swap) {
    if (first->base_matrices_count == 0) {
        return matrix_mxm_empty(output, first, second, accum, swap);
    }

    matrix_sort_lazy(first, true);
    matrix_combine_lazy(first, second->nvals);

    GrB_Matrix *accs = malloc(sizeof(GrB_Matrix) * first->base_matrices_count);
    Matrix *acc_matrices = malloc(sizeof(Matrix) * first->base_matrices_count);
    for (size_t i = 0; i < first->base_matrices_count; i++) {
        GrB_Matrix_new(&accs[i], GrB_BOOL, output->nrows, output->ncols);
        acc_matrices[i] = matrix_from_base(accs[i]);
    }

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        matrix_mxm_empty(&acc_matrices[i], &first->base_matrices[i], second, accum, swap);
    }

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        for (size_t j = i + 1; j < first->base_matrices_count; j++) {
            if (acc_matrices[i].nvals < acc_matrices[j].nvals) {
                Matrix temp = acc_matrices[i];
                acc_matrices[i] = acc_matrices[j];
                acc_matrices[j] = temp;
            }
        }
    }

    GrB_Matrix acc;
    GrB_Matrix_new(&acc, GrB_BOOL, first->nrows, first->ncols);
    Matrix acc_matrix = matrix_from_base(acc);

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        matrix_wise_empty(&acc_matrix, &acc_matrix, &acc_matrices[i], false);
        GrB_free(&acc_matrices[i].base);
    }

    if (accum) {
        return matrix_wise_empty(output, output, &acc_matrix, false);
    }

    GrB_Info result = matrix_dup_empty(output, &acc_matrix);
    GrB_free(&acc_matrix.base);

    return result;
}

GrB_Info matrix_mxm_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          bool swap) {
    if (first->block_type == CELL && second->block_type == CELL) {
        matrix_mxm_lazy(output, first, second, accum, swap);
    }

    if (first->block_type == CELL) {
        block_matrix_hyper_rotate_i(second, swap ? VEC_VERT : VEC_HORIZ);
        block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ);

        matrix_mxm_lazy(output, first, second, accum, swap);

        return GrB_SUCCESS;
    }

    if (second->block_type == CELL) {
        block_matrix_hyper_rotate_i(first, swap ? VEC_HORIZ : VEC_VERT);
        block_matrix_hyper_rotate_i(output, swap ? VEC_HORIZ : VEC_VERT);

        matrix_mxm_lazy(output, first, second, accum, swap);

        return GrB_SUCCESS;
    }

    GrB_Index size = first->nrows > first->ncols ? first->nrows : first->ncols;
    GrB_Matrix _diag;
    GrB_Matrix_new(&_diag, GrB_BOOL, size, size);
    Matrix diag = matrix_from_base(_diag);
    block_matrix_to_diag(&diag, second);

    block_matrix_hyper_rotate_i(first, swap ? VEC_VERT : VEC_HORIZ);
    block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ);

    matrix_mxm_lazy(output, first, &diag, accum, swap);
}

GrB_Info matrix_wise(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    GrB_BinaryOp accum_op = accum ? GxB_ANY_BOOL : GrB_NULL;

    GrB_Info result = GrB_eWiseAdd(output->base, GrB_NULL, accum_op, GxB_ANY_BOOL,
                                   first->base, second->base, GrB_NULL);

    matrix_update(output);
    return result;
}

GrB_Info matrix_wise_format(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (!output->is_both) {
        Matrix *larger = output->nvals > first->nvals ? output : first;
        larger = larger->nvals > second->nvals ? larger : second;

        matrix_to_format(output, larger->format, false);
        matrix_to_format(first, larger->format, false);
        matrix_to_format(second, larger->format, false);

        return matrix_wise(output, first, second, accum);
    }

    matrix_to_format(output, GrB_ROWMAJOR, false);
    matrix_to_format(first, output->format, false);
    matrix_to_format(second, output->format, false);
    GrB_Info result = matrix_wise(output, first, second, accum);

    if (result < GrB_SUCCESS) {
        return result;
    }

    matrix_to_format(output, GrB_COLMAJOR, false);
    matrix_to_format(first, output->format, false);
    matrix_to_format(second, output->format, false);
    return matrix_wise(output, first, second, accum);
}

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (first->nvals == 0 && second->nvals == 0) {
        if (accum) {
            return GrB_SUCCESS;
        }

        return matrix_clear_empty(output);
    }

    if (first->nvals == 0) {
        if (accum) {
            return matrix_wise_empty(output, output, second, false);
        }

        return matrix_dup_empty(output, second);
    }

    if (second->nvals == 0) {
        if (accum) {
            return matrix_wise_empty(output, output, first, false);
        }

        return matrix_dup_empty(output, first);
    }

    return matrix_wise_format(output, first, second, accum);
}

GrB_Info matrix_wise_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (first->base_matrices_count == 0) {
        first->base_matrices_count = 1;
        first->base_matrices[0] = matrix_from_base(first->base);
    }

    GrB_Matrix _other;
    GrB_Matrix_new(&_other, GrB_BOOL, output->nrows, output->ncols);
    Matrix other = matrix_from_base(_other);
    matrix_dup_empty(&other, second);

    size_t other_nvals = other.nvals >= 10 ? other.nvals : 10;

    while (true) {
        bool found = false;

        for (size_t i = 0; i < first->base_matrices_count; i++) {
            size_t self_nvals = first->base_matrices[i].nvals >= 10 ? first->nvals : 10;

            if (other.nvals / 10 <= self_nvals && self_nvals <= other.nvals * 10) {
                matrix_wise_empty(&other, &other, &first->base_matrices[i], false);
                GrB_free(&first->base_matrices[i].base);
                for (size_t j = i + 1; j < first->base_matrices_count; j++) {
                    first->base_matrices[j - 1] = first->base_matrices[j];
                }
                first->base_matrices_count--;
                found = true;
                break;
            }
        }

        if (found) {
            continue;
        }

        first->base_matrices[first->base_matrices_count++] = other;
        break;
    }

    return GrB_SUCCESS;
}

// - Any operation on two hyper vectors is performed block-wise
// - When hyper vector is added in-place to a cell, then sum of hyper vector's blocks is
// added to a cell
// - When cell is added in-place to a hyper vector, then cell is added to each of
// hyper vector's blocks
GrB_Info matrix_wise_block(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    if (output != first) {
        fprintf(stderr, "Matrix wise currently support only iadd operation");
        exit(-122);
    }

    if (first->block_type == CELL && second->block_type == CELL) {
        matrix_wise_lazy(output, first, second, accum);
    }

    // second is vector
    if (first->block_type == CELL) {
        Matrix temp_reduced = matrix_create(first->nrows, first->ncols);
        block_matrix_reduce(&temp_reduced, second);

        GrB_Info info = matrix_wise_lazy(output, first, &temp_reduced, accum);
        matrix_free(&temp_reduced);
        return info;
    }

    // first is vector
    if (second->block_type == CELL) {
        Matrix temp_vector = matrix_create(first->nrows, first->ncols);
        GrB_Index block_count = first->nrows > first->ncols ? first->nrows : first->ncols;
        block_matrix_repeat_into_vector(&temp_vector, second, block_count);

        GrB_Info info = matrix_wise_lazy(output, first, &temp_vector, accum);
        matrix_free(&temp_vector);
        return info;
    }

    // both are vector
    block_matrix_hyper_rotate_i(second, first->block_type);
    return matrix_wise_block(output, first, second, accum);
}

GrB_Info matrix_rsub(Matrix *output, Matrix *mask) {
    GrB_Info result = GrB_eWiseAdd(output->base, mask->base, GrB_NULL, GxB_ANY_BOOL,
                                   output->base, output->base, GrB_DESC_RSC);

    matrix_update(output);
    return result;
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

GrB_Info matrix_rsub_lazy(Matrix *output, Matrix *mask) {
    if (mask->base_matrices_count == 0) {
        return matrix_rsub_empty(output, mask);
    }

    matrix_sort_lazy(mask, true);
    matrix_combine_lazy(mask, output->nvals);

    for (size_t i = 0; i < mask->base_matrices_count; i++) {
        matrix_rsub_empty(output, &mask->base_matrices[i]);
    }

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_block(Matrix *output, Matrix *mask) {
    if (output->block_type == CELL && output->block_type != CELL ||
        output->block_type != CELL && mask->block_type == CELL) {
        fprintf(stderr, "Don't support rsub operation between cell and vector");
        exit(-1);
    }

    if (output->block_type == CELL) {
        return matrix_rsub_lazy(output, mask);
    }

    block_matrix_hyper_rotate_i(output, mask->block_type);
    return matrix_rsub_lazy(output, mask);
}

void matrix_print_lazy(Matrix *A) {
    if (A->base_matrices_count == 0) {
        GxB_print(A->base, 1);
        return;
    }

    if (A->base_matrices_count == 1) {
        GxB_print(A->base_matrices[0].base, 1);
        return;
    }

    GrB_Matrix _temp;
    GrB_Matrix_new(&_temp, GrB_BOOL, A->nrows, A->ncols);
    Matrix temp = matrix_from_base(_temp);
    for (size_t i = 0; i < A->base_matrices_count; i++) {
        matrix_wise_empty(&temp, &temp, &A->base_matrices[i], false);
    }

    A = &temp;
    GxB_print(A->base, 1);
    GrB_free(&_temp);
}

#define IS_NONTERM(index)                                                                \
    {                                                                                    \
        for (size_t m = 0; m < rules_count; m++) {                                       \
            if (rules[m].nonterm == index)                                               \
                return true;                                                             \
        }                                                                                \
                                                                                         \
        return false;                                                                    \
    }

bool is_nonterm(int index, const LAGraph_rule_WCNF *rules, size_t rules_count) {
    for (size_t i = 0; i < rules_count; i++) {
        if (rules[i].nonterm == index) {
            return true;
        }
    }

    return false;
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
    int32_t nonterms_count,
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    size_t rules_count,             // The total number of rules in the CFG.
    char *msg,                      // Message string for error reporting.
    int8_t optimizations            // Optimizations flags
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

    int32_t symbols_amount = 0;
    for (size_t i = 0; i < rules_count; i++) {
        symbols_amount =
            rules[i].nonterm + 1 > symbols_amount ? rules[i].nonterm + 1 : symbols_amount;
        symbols_amount =
            rules[i].prod_A + 1 > symbols_amount ? rules[i].prod_A + 1 : symbols_amount;
        symbols_amount =
            rules[i].prod_B + 1 > symbols_amount ? rules[i].prod_B + 1 : symbols_amount;
    }

    GrB_Scalar true_scalar;
    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);

    LG_TRY(LAGraph_Calloc((void **)&T, symbols_amount, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&delta_matrices, symbols_amount, sizeof(Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&matrices, symbols_amount, sizeof(Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&temp_matrices, symbols_amount, sizeof(Matrix), msg));

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
    for (int32_t i = 0; i < symbols_amount; i++) {
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
    for (int32_t i = 0; i < symbols_amount; i++) {
        GrB_Matrix matrix;

        GRB_TRY(GrB_Matrix_new(&T[i], GrB_BOOL, n, n));

        GRB_TRY(GrB_Matrix_new(&matrix, GrB_BOOL, n, n));
        delta_matrices[i] = matrix_from_base(matrix);

        GrB_Matrix_dup(&matrices[i].base, adj_matrices[i]);

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
        if (rule.nonterm < 0 || rule.nonterm >= symbols_amount) {
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

            if (rule.prod_A < -1 || rule.prod_A >= symbols_amount) {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin) {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= symbols_amount || rule.prod_B < -1 ||
                rule.prod_B >= symbols_amount) {
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

    typedef GrB_Info (*matrix_mxm_fn)(Matrix *output, Matrix *first, Matrix *second,
                                      bool accum, bool swap);
    typedef GrB_Info (*matrix_wise_fn)(Matrix *output, Matrix *first, Matrix *second,
                                       bool accum);
    typedef GrB_Info (*matrix_rsub_fn)(Matrix *input, Matrix *mask);

    matrix_mxm_fn mxm;
    matrix_wise_fn wise;
    matrix_rsub_fn rsub;

    if (optimizations & OPT_LAZY) {
        mxm = matrix_mxm_lazy;
        wise = matrix_wise_lazy;
        rsub = matrix_rsub_lazy;
    } else if (optimizations & OPT_EMPTY) {
        mxm = matrix_mxm_empty;
        wise = matrix_wise_empty;
        rsub = matrix_rsub_empty;
    } else if (optimizations & OPT_FORMAT) {
        mxm = matrix_mxm_format;
        wise = matrix_wise_format;
        rsub = matrix_rsub_format;
    } else {
        mxm = matrix_mxm;
        wise = matrix_wise;
        rsub = matrix_rsub;
    }

    double start_time, end_time;
    bool changed = true;
    size_t iteration = 0;
    double mxm1 = 0.0;
    double wise1 = 0.0;
    double mxm2 = 0.0;
    double wise2 = 0.0;
    double rsubt = 0.0;
    while (changed) {
        iteration++;
        changed = false;

#if BENCH_CFL_REACHBILITY
        printf("\n--- ITERATARION %ld ---\n", iteration);
#endif

        for (int32_t i = 0; i < nonterms_count; i++) {
            GRB_TRY(matrix_clear_empty(&temp_matrices[i]));
        }

        TIMER_START();
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];
            Matrix *A = &matrices[bin_rule.prod_A];
            Matrix *B = &delta_matrices[bin_rule.prod_B];
            Matrix *C = &temp_matrices[bin_rule.nonterm];

            mxm(C, A, B, true, false);
            // matrix_print_lazy(A);
            // matrix_print_lazy(B);
            // matrix_print_lazy(C);
        }
        TIMER_STOP("MXM 1", &mxm1);

        TIMER_START()
        for (int32_t i = 0; i < nonterms_count; i++) {
            Matrix *A = &delta_matrices[i];
            Matrix *C = &matrices[i];

            wise(C, C, A, false);
            // matrix_print_lazy(A);
            // matrix_print_lazy(C);
        }
        TIMER_STOP("WISE 1", &wise1);

        TIMER_START()
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];
            Matrix *A = &matrices[bin_rule.prod_B];
            Matrix *B = &delta_matrices[bin_rule.prod_A];
            Matrix *C = &temp_matrices[bin_rule.nonterm];
            // printf("ITER: %ld\n", i);
            mxm(C, A, B, true, true);
            // matrix_print_lazy(A);
            // matrix_print_lazy(B);
            // matrix_print_lazy(C);
        }
        TIMER_STOP("MXM 2", &mxm2);

        TIMER_START();
        for (int32_t i = 0; i < nonterms_count; i++) {
            matrix_dup_empty(&delta_matrices[i], &temp_matrices[i]);
        }
        TIMER_STOP("WISE 2 (copy)", &wise2);

        TIMER_START();
        for (int32_t i = 0; i < nonterms_count; i++) {
            Matrix *A = &matrices[i];
            Matrix *C = &delta_matrices[i];

            rsub(C, A);
            // matrix_print_lazy(A);
            // matrix_print_lazy(C);
        }
        TIMER_STOP("WISE 3 (MASK)", &rsubt);

        for (int32_t i = 0; i < nonterms_count; i++) {
            size_t new_nnz = 0;
            if (matrices[i].base_matrices_count == 0) {
                new_nnz = matrices[i].nvals;
            } else {
                for (size_t j = 0; j < matrices[i].base_matrices_count; j++) {
                    new_nnz += matrices[i].base_matrices[j].nvals;
                }
            }

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
           mxm2, wise2, rsubt);
#endif

#ifdef DEBUG_CFL_REACHBILITY
    for (int32_t i = 0; i < nonterms_count; i++) {
        printf("MATRIX WITH INDEX %d:\n", i);
        GxB_print(T[i], GxB_SUMMARY);
    }
#endif

    for (int32_t i = 0; i < nonterms_count; i++) {
        if (matrices[i].base_matrices_count == 0) {
            outputs[i] = matrices[i].base;
        } else {
            GrB_Matrix _acc;
            GrB_Matrix_new(&_acc, GrB_BOOL, matrices[i].nrows, matrices[i].ncols);
            Matrix acc = matrix_from_base(_acc);

            for (size_t j = 0; j < matrices[i].base_matrices_count; j++) {
                matrix_wise_empty(&acc, &acc, &matrices[i].base_matrices[j], false);
                GrB_free(&matrices[i].base_matrices[j].base);
            }

            outputs[i] = acc.base;
        }
        // outputs[i] = matrices[i].base;
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
