//------------------------------------------------------------------------------
// LAGraph_CFL_optimized_matrix_opt.c: Implementation of operations for
// Optimized Context-Free Language Reachability Matrix-Based Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vlasenco Daniel, Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

// Code is an implementation of optimized matrix operations for CFL_reachability
// algorithms, described in the following paper:
// * Ilia Muravev, "Optimization of the Context-Free Language Reachability Matrix-Based
// Algorithm" and based on the python implementation from:
// https://github.com/FormalLanguageConstrainedPathQuerying/CFPQ_PyAlgo/tree/murav/optimize-matrix

#include "LG_internal.h"
#include <LAGraphX.h>

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

typedef CFL_Matrix Matrix;
typedef enum CFL_Matrix_block Matrix_block;

#define TO_COL(matrix) GrB_set(matrix, GrB_COLMAJOR, GrB_STORAGE_ORIENTATION_HINT)
#define TO_ROW(matrix) GrB_set(matrix, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)

GrB_Info matrix_to_format(Matrix *matrix, int32_t format, bool is_both) {
    // Matrix contain both formats so just switch base matrix
    if (matrix->is_both) {
        matrix->base = format == GrB_ROWMAJOR ? matrix->base_row : matrix->base_col;
        matrix->format = format;
        return GrB_SUCCESS;
    }

    // No changes required
    if (matrix->format == format) {
        return GrB_SUCCESS;
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
    return GrB_SUCCESS;
}

// clear methods

GrB_Info matrix_clear(Matrix *A) {
    GrB_Info result = GrB_Matrix_clear(A->base);
    result;
    CFL_matrix_update(A);
    return result;
}

GrB_Info matrix_clear_format(Matrix *A, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        return matrix_clear(A);
    }

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

GrB_Info matrix_clear_empty(Matrix *A, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        return matrix_clear_format(A, optimizations);
    }

    if (A->nvals == 0) {
        return GrB_SUCCESS;
    }

    return matrix_clear_format(A, optimizations);
}

// duplicate methods

GrB_Info matrix_dup(Matrix *output, Matrix *input) {
    if (output == input) {
        return GrB_SUCCESS;
    }

    GrB_Info result = GrB_Matrix_apply(output->base, GrB_NULL, GrB_NULL,
                                       GrB_IDENTITY_BOOL, input->base, GrB_NULL);
    CFL_matrix_update(output);

    return result;
}

GrB_Info matrix_dup_format(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        return matrix_dup(output, input);
    }

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

GrB_Info matrix_dup_empty(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        return matrix_dup_format(output, input, optimizations);
    }

    if (input->nvals == 0) {
        return matrix_clear_empty(output, optimizations);
    }

    return matrix_dup_format(output, input, optimizations);
}

GrB_Info block_matrix_hyper_rotate_i(Matrix *matrix, enum CFL_Matrix_block format);

GrB_Info matrix_dup_block(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        return matrix_dup_empty(output, input, optimizations);
    }

    if (output->block_type == CELL && input->block_type == CELL) {
        return matrix_dup_empty(output, input, optimizations);
    }

    block_matrix_hyper_rotate_i(input, output->block_type);
    return matrix_dup_empty(output, input, optimizations);
}

// block optimization specific methods

GrB_Info block_matrix_hyper_rotate_i(Matrix *matrix, enum CFL_Matrix_block format) {
    if (matrix->is_lazy) {
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            block_matrix_hyper_rotate_i(&matrix->base_matrices[i], format);
        }

        CFL_matrix_update(matrix);
        return GrB_SUCCESS;
    }

    if (matrix->block_type == CELL) {
        return GrB_SUCCESS;
    }

    if (matrix->block_type == format) {
        return GrB_SUCCESS;
    }

    GrB_Scalar scalar_true;
    GrB_Scalar_new(&scalar_true, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(scalar_true, true);

    if (matrix->block_type == VEC_VERT) {
        // fix: change to lagraph malloc
        GrB_Index *nrows = malloc(matrix->nvals * sizeof(GrB_Index));
        GrB_Index *ncols = malloc(matrix->nvals * sizeof(GrB_Index));

        GrB_Matrix_extractTuples_BOOL(nrows, ncols, NULL, &matrix->nvals, matrix->base);

        for (size_t i = 0; i < matrix->nvals; i++) {
            ncols[i] = ncols[i] + nrows[i] / matrix->ncols * matrix->ncols;
            nrows[i] = nrows[i] % matrix->ncols;
        }

        GrB_Matrix new;
        GrB_Matrix_new(&new, GrB_BOOL, matrix->ncols, matrix->nrows);
        GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals);
        CFL_matrix_free(matrix);
        *matrix =
            matrix->is_lazy ? CFL_matrix_from_base_lazy(new) : CFL_matrix_from_base(new);
        free(nrows);
        free(ncols);
        GrB_free(&scalar_true);
        return GrB_SUCCESS;
    }

    if (matrix->block_type == VEC_HORIZ) {
        GrB_Index *nrows = malloc(matrix->nvals * sizeof(GrB_Index));
        GrB_Index *ncols = malloc(matrix->nvals * sizeof(GrB_Index));

        GrB_Matrix_extractTuples_BOOL(nrows, ncols, NULL, &matrix->nvals, matrix->base);

        for (size_t i = 0; i < matrix->nvals; i++) {
            nrows[i] = nrows[i] + ncols[i] / matrix->nrows * matrix->nrows;
            ncols[i] = ncols[i] % matrix->nrows;
        }

        GrB_Matrix new;
        GrB_Matrix_new(&new, GrB_BOOL, matrix->ncols, matrix->nrows);
        GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals);
        CFL_matrix_free(matrix);
        *matrix =
            matrix->is_lazy ? CFL_matrix_from_base_lazy(new) : CFL_matrix_from_base(new);
        free(nrows);
        free(ncols);
        GrB_free(&scalar_true);
        return GrB_SUCCESS;
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

void block_matrix_reduce(Matrix *matrix, Matrix *input, int8_t optimizations) {
    if (input->block_type == CELL) {
        matrix_dup_block(matrix, input, optimizations);
    }

    GrB_Scalar scalar_true;
    GrB_Scalar_new(&scalar_true, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(scalar_true, true);

    GrB_Index *rows = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Index *cols = malloc(input->nvals * sizeof(GrB_Index));
    GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &input->nvals, input->base);

    if (input->block_type == VEC_VERT) {
        for (size_t i = 0; i < input->nvals; i++) {
            rows[i] = rows[i] % input->ncols;
        }
    }

    if (input->block_type == VEC_HORIZ) {
        for (size_t i = 0; i < input->nvals; i++) {
            cols[i] = cols[i] % input->nrows;
        }
    }

    GxB_Matrix_build_Scalar(matrix->base, rows, cols, scalar_true, input->nvals);
    CFL_matrix_update(matrix);

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
    CFL_matrix_update(matrix);
}

// lazy optimization specific methods

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations);

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

    return GrB_SUCCESS;
}

GrB_Matrix CFL_matrix_lazy_to_base(Matrix *matrix, int8_t optimizations) {
    GrB_Matrix _acc;
    GrB_Matrix_new(&_acc, GrB_BOOL, matrix->nrows, matrix->ncols);
    Matrix acc = CFL_matrix_from_base(_acc);

    matrix_sort_lazy(matrix, false);
    for (size_t j = 0; j < matrix->base_matrices_count; j++) {
        CFL_wise(&acc, &acc, &matrix->base_matrices[j], false, optimizations);
        GrB_free(&matrix->base_matrices[j].base);
    }

    return acc.base;
}

GrB_Info matrix_combine_lazy(Matrix *A, size_t threshold, int8_t optimizations) {
    Matrix *new_matrices = malloc(sizeof(Matrix) * 50);
    size_t new_size = 0;

    matrix_sort_lazy(A, false);

    for (size_t i = 0; i < A->base_matrices_count; i++) {
        if (new_size == 0 || A->base_matrices[i].nvals > threshold) {
            new_matrices[new_size++] = A->base_matrices[i];
            continue;
        }

        matrix_wise_empty(&new_matrices[new_size - 1], &new_matrices[new_size - 1],
                          &A->base_matrices[i], false, optimizations);
        GrB_free(&A->base_matrices[i].base);
    }

    A->base_matrices = new_matrices;
    A->base_matrices_count = new_size;
    CFL_matrix_update(A);

    return GrB_SUCCESS;
}

// create and update methods

void CFL_matrix_update(Matrix *matrix) {
    if (!matrix->is_lazy) {
        GrB_Matrix_nvals(&matrix->nvals, matrix->base);
    } else {
        size_t new_nnz = 0;
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            GrB_Matrix_nvals(&matrix->base_matrices[i].nvals,
                             matrix->base_matrices[i].base);
            new_nnz += matrix->base_matrices[i].nvals;
        }

        matrix->nvals = new_nnz;
    }

    if (!matrix->is_lazy) {
        GrB_Matrix_nrows(&matrix->nrows, matrix->base);
        GrB_Matrix_ncols(&matrix->ncols, matrix->base);
    } else {
        GrB_Matrix_nrows(&matrix->nrows, matrix->base_matrices[0].base);
        GrB_Matrix_ncols(&matrix->ncols, matrix->base_matrices[0].base);
    }

    if (matrix->nrows == matrix->ncols)
        matrix->block_type = CELL;
    else
        matrix->block_type = matrix->nrows > matrix->ncols ? VEC_VERT : VEC_HORIZ;
}

// Create optimizied matrix from base GrB_Matrix
Matrix CFL_matrix_from_base(GrB_Matrix matrix) {
    Matrix result;

    result.base = matrix;
    result.nvals = 0; // We will get actual info in update functoin
    result.nrows = 0;
    result.ncols = 0;

    // Format optimization fields
    result.base_row = matrix;
    result.base_col = NULL;
    result.format = GrB_ROWMAJOR;
    result.is_both = false;

    // Lazy addition optimization fields
    result.is_lazy = false;
    result.base_matrices = malloc(sizeof(CFL_Matrix) * 40); // TODO: dynamic size
    result.base_matrices_count = 0;

    // Block optimization fields
    result.block_type = CELL;

    CFL_matrix_update(&result);
    return result;
}

Matrix CFL_matrix_from_base_lazy(GrB_Matrix matrix) {
    Matrix result = CFL_matrix_from_base(matrix);

    Matrix lazy_result = CFL_matrix_from_base(matrix);
    lazy_result.is_lazy = true;
    lazy_result.base_matrices[0] = result;
    lazy_result.base_matrices_count = 1;
    CFL_matrix_update(&lazy_result);

    return lazy_result;
}

Matrix CFL_matrix_create(GrB_Index nrows, GrB_Index ncols) {
    GrB_Matrix _result;
    GrB_Matrix_new(&_result, GrB_BOOL, nrows, ncols);
    Matrix result = CFL_matrix_from_base(_result);

    return result;
}

// TODO: free all base_matrices, free format matrices
void CFL_matrix_free(Matrix *matrix) {
    free(matrix->base_matrices);
    GrB_free(&matrix->base);
}

// mxm operations

GrB_Info matrix_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum,
                    bool swap) {
    Matrix *left = swap ? second : first;
    Matrix *right = swap ? first : second;

    GrB_Info result = GrB_mxm(output->base, GrB_NULL, accum ? GxB_ANY_BOOL : GrB_NULL,
                              GxB_ANY_PAIR_BOOL, left->base, right->base, GrB_NULL);
    // IS_ISO(output->base, "MXM output");
    CFL_matrix_update(output);
    return result;
}

GrB_Info matrix_mxm_format(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        return matrix_mxm(output, first, second, accum, swap);
    }

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
                          bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        return matrix_mxm_format(output, first, second, accum, swap, optimizations);
    }

    if (first->nvals == 0 || second->nvals == 0) {
        if (accum) {
            return GrB_SUCCESS;
        }

        if (output->nvals == 0) {
            return GrB_SUCCESS;
        }

        matrix_clear_empty(output, optimizations);
    }

    return matrix_mxm_format(output, first, second, accum, swap, optimizations);
}

GrB_Info matrix_mxm_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum,
                         bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        return matrix_mxm_empty(output, first, second, accum, swap, optimizations);
    }

    if (!first->is_lazy) {
        return matrix_mxm_empty(output, first, second, accum, swap, optimizations);
    }

    matrix_combine_lazy(first, second->nvals, optimizations);
    matrix_sort_lazy(first, false);

    GrB_Matrix *accs = malloc(sizeof(GrB_Matrix) * first->base_matrices_count);
    Matrix *acc_matrices = malloc(sizeof(Matrix) * first->base_matrices_count);
    for (size_t i = 0; i < first->base_matrices_count; i++) {
        GrB_Matrix_new(&accs[i], GrB_BOOL, swap ? second->nrows : first->nrows,
                       swap ? first->ncols : second->ncols);
        acc_matrices[i] = CFL_matrix_from_base(accs[i]);
    }

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        matrix_mxm_empty(&acc_matrices[i], &first->base_matrices[i], second, false, swap,
                         optimizations);
    }

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        for (size_t j = i + 1; j < first->base_matrices_count; j++) {
            if (acc_matrices[i].nvals > acc_matrices[j].nvals) {
                Matrix temp = acc_matrices[i];
                acc_matrices[i] = acc_matrices[j];
                acc_matrices[j] = temp;
            }
        }
    }

    GrB_Matrix acc;
    GrB_Matrix_new(&acc, GrB_BOOL, swap ? second->nrows : first->nrows,
                   swap ? first->ncols : second->ncols);
    Matrix acc_matrix = CFL_matrix_from_base(acc);

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        matrix_wise_empty(&acc_matrix, &acc_matrix, &acc_matrices[i], false,
                          optimizations);
        GrB_free(&acc_matrices[i].base);
    }

    if (accum) {
        return matrix_wise_empty(output, output, &acc_matrix, false, optimizations);
    }

    GrB_Info result = matrix_dup_block(output, &acc_matrix, optimizations);
    GrB_free(&acc_matrix.base);

    return result;
}

GrB_Info matrix_wise_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations);

GrB_Info matrix_mxm_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        return matrix_mxm_lazy(output, first, second, accum, swap, optimizations);
    }

    if (first->block_type == CELL && second->block_type == CELL) {
        return matrix_mxm_lazy(output, first, second, accum, swap, optimizations);
    }

    if (first->block_type == CELL) {
        block_matrix_hyper_rotate_i(second, swap ? VEC_VERT : VEC_HORIZ);
        block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ);

        Matrix temp = CFL_matrix_create(swap ? second->nrows : first->nrows,
                                        swap ? first->ncols : second->ncols);
        matrix_mxm_lazy(&temp, first, second, accum, swap, optimizations);
        matrix_wise_block(output, output, &temp, false, optimizations);
        CFL_matrix_free(&temp);

        return GrB_SUCCESS;
    }

    if (second->block_type == CELL) {
        block_matrix_hyper_rotate_i(first, swap ? VEC_HORIZ : VEC_VERT);
        block_matrix_hyper_rotate_i(output, swap ? VEC_HORIZ : VEC_VERT);

        Matrix temp = CFL_matrix_create(swap ? second->nrows : first->nrows,
                                        swap ? first->ncols : first->ncols);
        matrix_mxm_lazy(&temp, first, second, accum, swap, optimizations);
        matrix_wise_block(output, output, &temp, false, optimizations);
        CFL_matrix_free(&temp);

        return GrB_SUCCESS;
    }

    GrB_Index size = first->nrows > first->ncols ? first->nrows : first->ncols;
    GrB_Matrix _diag;
    GrB_Matrix_new(&_diag, GrB_BOOL, size, size);
    Matrix diag = CFL_matrix_from_base(_diag);
    block_matrix_to_diag(&diag, second);
    CFL_matrix_update(&diag);

    block_matrix_hyper_rotate_i(first, swap ? VEC_VERT : VEC_HORIZ);
    block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ);

    Matrix temp = CFL_matrix_create(first->nrows, diag.ncols);
    matrix_mxm_lazy(&temp, first, &diag, false, swap, optimizations);
    return matrix_wise_block(output, output, &temp, false, optimizations);
}

// wise operations

GrB_Info matrix_wise(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    GrB_BinaryOp accum_op = accum ? GxB_ANY_BOOL : GrB_NULL;
    if (output == first)
        accum_op = GrB_NULL;

    GrB_Info result = GrB_eWiseAdd(output->base, GrB_NULL, accum_op, GxB_ANY_BOOL,
                                   first->base, second->base, GrB_NULL);

    CFL_matrix_update(output);
    return result;
}

GrB_Info matrix_wise_format(Matrix *output, Matrix *first, Matrix *second, bool accum,
                            int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        return matrix_wise(output, first, second, accum);
    }

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

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        return matrix_wise_format(output, first, second, accum, optimizations);
    }

    if (output == first) {
        if (first->nvals == 0) {
            return matrix_dup_empty(first, second, optimizations);
        }

        if (second->nvals == 0) {
            return GrB_SUCCESS;
        }

        return matrix_wise_format(output, first, second, accum, optimizations);
    }

    if (first->nvals == 0 && second->nvals == 0) {
        if (accum || output->nvals == 0) {
            return GrB_SUCCESS;
        }

        return matrix_clear_empty(output, optimizations);
    }

    if (first->nvals == 0) {
        if (accum) {
            return matrix_wise_empty(output, output, second, false, optimizations);
        }

        return matrix_dup_empty(output, second, optimizations);
    }

    if (second->nvals == 0) {
        if (accum) {
            return matrix_wise_empty(output, output, first, false, optimizations);
        }

        return matrix_dup_empty(output, first, optimizations);
    }

    return matrix_wise_format(output, first, second, accum, optimizations);
}

GrB_Info matrix_wise_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        return matrix_wise_empty(output, first, second, accum, optimizations);
    }

    if (!first->is_lazy && !second->is_lazy) {
        return matrix_wise_empty(output, first, second, accum, optimizations);
    }

    if (!first->is_lazy && second->is_lazy) {
        for (size_t i = 0; i < second->base_matrices_count; i++) {
            matrix_wise_empty(output, first, &second->base_matrices[i], false,
                              optimizations);
        }

        return GrB_SUCCESS;
    }

    GrB_Matrix _other;
    GrB_Matrix_new(&_other, GrB_BOOL, output->nrows, output->ncols);
    Matrix other = CFL_matrix_from_base(_other);
    matrix_dup_empty(&other, second, optimizations);

    size_t other_nvals = other.nvals >= 10 ? other.nvals : 10;

    while (true) {
        bool found = false;

        for (size_t i = 0; i < first->base_matrices_count; i++) {
            size_t self_nvals =
                first->base_matrices[i].nvals >= 10 ? first->base_matrices[i].nvals : 10;

            if (other_nvals / 10 <= self_nvals && self_nvals <= other_nvals * 10) {
                matrix_wise_empty(&other, &other, &first->base_matrices[i], accum,
                                  optimizations);
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
GrB_Info matrix_wise_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        return matrix_wise_lazy(output, first, second, accum, optimizations);
    }

    if (output != first) {
        fprintf(stderr, "Matrix wise currently support only iadd operation");
        exit(-122);
    }

    if (first->block_type == CELL && second->block_type == CELL) {
        return matrix_wise_lazy(output, first, second, accum, optimizations);
    }

    // second is vector
    if (first->block_type == CELL) {
        Matrix temp_reduced = CFL_matrix_create(first->nrows, first->ncols);
        block_matrix_reduce(&temp_reduced, second, optimizations);

        GrB_Info info =
            matrix_wise_lazy(output, first, &temp_reduced, accum, optimizations);
        CFL_matrix_free(&temp_reduced);
        return info;
    }

    // first is vector
    if (second->block_type == CELL) {
        // LG_SET_BURBLE(true);
        Matrix temp_vector = CFL_matrix_create(first->nrows, first->ncols);
        GrB_Index block_count = first->nrows > first->ncols ? first->nrows : first->ncols;
        block_matrix_repeat_into_vector(&temp_vector, second, block_count);

        block_matrix_hyper_rotate_i(&temp_vector, first->block_type);
        GrB_Info info =
            matrix_wise_lazy(output, first, &temp_vector, accum, optimizations);
        CFL_matrix_free(&temp_vector);
        return info;
    }

    // both are vector
    block_matrix_hyper_rotate_i(second, first->block_type);
    return matrix_wise_lazy(output, first, second, accum, optimizations);
}

// rsub methods

GrB_Info matrix_rsub(Matrix *output, Matrix *mask) {
    GrB_Info result = GrB_eWiseAdd(output->base, mask->base, GrB_NULL, GxB_ANY_BOOL,
                                   output->base, output->base, GrB_DESC_RSC);

    CFL_matrix_update(output);
    return result;
}

GrB_Info matrix_rsub_format(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        return matrix_rsub(output, mask);
    }

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

GrB_Info matrix_rsub_empty(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        return matrix_rsub_format(output, mask, optimizations);
    }

    if (mask->nvals == 0 || output->nvals == 0) {
        return GrB_SUCCESS;
    }

    return matrix_rsub_format(output, mask, optimizations);
}

GrB_Info matrix_rsub_lazy(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        return matrix_rsub_empty(output, mask, optimizations);
    }

    if (!mask->is_lazy) {
        return matrix_rsub_empty(output, mask, optimizations);
    }

    matrix_combine_lazy(mask, output->nvals, optimizations);
    matrix_sort_lazy(mask, true);

    for (size_t i = 0; i < mask->base_matrices_count; i++) {
        matrix_rsub_empty(output, &mask->base_matrices[i], optimizations);
    }

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_block(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        return matrix_rsub_lazy(output, mask, optimizations);
    }

    if ((output->block_type == CELL && mask->block_type != CELL) ||
        (output->block_type != CELL && mask->block_type == CELL)) {
        fprintf(stderr, "Don't support rsub operation between cell and vector");
        exit(-1);
    }

    if (output->block_type == CELL) {
        return matrix_rsub_lazy(output, mask, optimizations);
    }

    block_matrix_hyper_rotate_i(output, mask->block_type);
    return matrix_rsub_lazy(output, mask, optimizations);
}

// utility methods

// void matrix_print_lazy(Matrix *A) {
//     return;
//     GxB_Print_Level pr = 1;

//     if (!A->is_lazy) {
//         CFL_matrix_update(A);
//         GxB_print(A->base, pr);
//         // printf("nnz: %ld\n", A->nvals);
//         return;
//     }

//     if (A->base_matrices_count == 1) {
//         CFL_matrix_update(A);
//         // printf("nnz: %ld\n", A->nvals);
//         GxB_print(A->base_matrices[0].base, pr);
//         return;
//     }

//     GrB_Matrix _temp;
//     GrB_Matrix_new(&_temp, GrB_BOOL, A->nrows, A->ncols);
//     Matrix temp = CFL_matrix_from_base(_temp);
//     for (size_t i = 0; i < A->base_matrices_count; i++) {
//         matrix_wise_empty(&temp, &temp, &A->base_matrices[i], false, optimizations);
//     }

//     A = &temp;
//     GxB_print(A->base, pr);
//     // printf("nnz: %ld\n", A->nvals);
//     GrB_free(&_temp);
// }

void print_graph_info(Matrix *matrices, size_t count) {
    GrB_Index nnz = 0;

    for (size_t i = 0; i < count; i++) {
        Matrix *A = &matrices[i];
        CFL_matrix_update(A);
        nnz += A->nvals;
    }

    printf("NNZ: %ld\n", nnz);
}

// order of optimizations: block -> lazy -> empty -> format

GrB_Info CFL_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum, bool swap,
                 int8_t optimizations) {
    return matrix_mxm_block(output, first, second, accum, swap, optimizations);
}

GrB_Info CFL_wise(Matrix *output, Matrix *first, Matrix *second, bool accum,
                  int8_t optimizations) {
    return matrix_wise_block(output, first, second, accum, optimizations);
}

GrB_Info CFL_rsub(Matrix *output, Matrix *mask, int8_t optimizations) {
    return matrix_rsub_block(output, mask, optimizations);
}

GrB_Info CFL_dup(Matrix *output, Matrix *input, int8_t optimizations) {
    return matrix_dup_block(output, input, optimizations);
}

GrB_Info CFL_clear(Matrix *A, int8_t optimizations) {
    return matrix_clear_empty(A, optimizations);
}