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