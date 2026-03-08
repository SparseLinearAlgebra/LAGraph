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

#include "LAGraph_CFL_optimized_matrix_opt.h"
#include "LG_internal.h"
#include <LAGraphX.h>

#define BENCH_CFL_REACHBILITY false

#if BENCH_CFL_REACHBILITY
    #define IS_ISO(matrix, str)                                                          \
        {                                                                                \
            bool iso_flag;                                                               \
            GrB_Index nnz;                                                               \
            TRY(GxB_Matrix_iso(&iso_flag, matrix));                                      \
            TRY(GrB_Matrix_nvals(&nnz, matrix));                                         \
            if (!iso_flag && nnz) {                                                      \
                printf("-----ISO ALERT----- (%s)\n", str);                               \
                TRY(GxB_print(matrix, 1));                                               \
                printf("-------------------\n");                                         \
                TRY(-122);                                                               \
            }                                                                            \
        }

    #define IS_ROW(matrix, str)                                                          \
        {                                                                                \
            if (matrix != NULL) {                                                        \
                int32_t orientation;                                                     \
                TRY(GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT));        \
                if (orientation != GrB_ROWMAJOR) {                                       \
                    printf("-----NOT A ROW----- (%s)\n", str);                           \
                    TRY(GxB_print(matrix, 1));                                           \
                    printf("-------------------\n");                                     \
                    TRY(-122);                                                           \
                }                                                                        \
            }                                                                            \
        }

    #define IS_COL(matrix, str)                                                          \
        {                                                                                \
            if (matrix != NULL) {                                                        \
                int32_t orientation;                                                     \
                TRY(GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT));        \
                if (orientation != GrB_COLMAJOR) {                                       \
                    printf("-----NOT A COL----- (%s)\n", str);                           \
                    TRY(GxB_print(matrix, 1));                                           \
                    printf("-------------------\n");                                     \
                    TRY(-122);                                                           \
                }                                                                        \
            }                                                                            \
        }
#else
    #define IS_ISO(matrix, str)
    #define IS_ROW(matrix, str)
    #define IS_COL(matrix, str)
#endif

#define TRY(GrB_method)                                                                  \
    {                                                                                    \
        GrB_Info LG_GrB_Info = GrB_method;                                               \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                 \
            fprintf(stderr, "LAGraph failure (file %s, line %d): (%d) \n", __FILE__,     \
                    __LINE__, LG_GrB_Info);                                              \
            return (LG_GrB_Info);                                                        \
        }                                                                                \
    }

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

typedef CFL_Matrix Matrix;
typedef enum CFL_Matrix_block Matrix_block;

#define TO_COL(matrix) GrB_set(matrix, GrB_COLMAJOR, GrB_STORAGE_ORIENTATION_HINT)
#define TO_ROW(matrix) GrB_set(matrix, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)

GrB_Info matrix_print_lazy(Matrix *A, int8_t optimizations);

GrB_Info matrix_to_format(Matrix *matrix, int32_t format, bool is_both) {
    TRY(CFL_matrix_update(matrix));

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
    GrB_Matrix new_matrix, old_matrix;
    if (matrix->format == GrB_ROWMAJOR) {
        new_matrix = matrix->base_col;
        old_matrix = matrix->base_row;
    } else {
        new_matrix = matrix->base_row;
        old_matrix = matrix->base_col;
    }

    if (is_both) {
        TRY(GrB_Matrix_new(&new_matrix, GrB_BOOL, matrix->nrows, matrix->ncols));
        TRY(GrB_set(new_matrix, format, GrB_STORAGE_ORIENTATION_HINT))
        TRY(GrB_Matrix_assign(new_matrix, GrB_NULL, GrB_NULL, old_matrix, GrB_ALL,
                              matrix->nrows, GrB_ALL, matrix->ncols, GrB_NULL));
        if (format == GrB_ROWMAJOR) {
            matrix->base_row = new_matrix;
            matrix->base_col = matrix->base;
        } else {
            matrix->base_row = matrix->base;
            matrix->base_col = new_matrix;
        }
        matrix->is_both = true;
        matrix->base = new_matrix;
    } else {
        if (format == GrB_ROWMAJOR) {
            matrix->base_row = matrix->base;
            matrix->base_col = NULL;
            TRY(TO_ROW(matrix->base));
        } else {
            matrix->base_row = NULL;
            matrix->base_col = matrix->base;
            TRY(TO_COL(matrix->base));
        }
    }

    matrix->format = format;

    TRY(CFL_matrix_update(matrix));

    return GrB_SUCCESS;
}

// clear methods

GrB_Info matrix_clear(Matrix *A) {
    TRY(GrB_Matrix_clear(A->base));
    TRY(CFL_matrix_update(A));
    return GrB_SUCCESS;
}

GrB_Info matrix_clear_format(Matrix *A, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        TRY(matrix_clear(A));
        return GrB_SUCCESS;
    }

    if (!A->is_both) {
        TRY(matrix_clear(A));
        return GrB_SUCCESS;
    }

    TRY(matrix_to_format(A, GrB_ROWMAJOR, false));
    TRY(matrix_clear(A));

    TRY(matrix_to_format(A, GrB_COLMAJOR, false));
    TRY(matrix_clear(A));

    return GrB_SUCCESS;
}

GrB_Info matrix_clear_empty(Matrix *A, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        TRY(matrix_clear_format(A, optimizations));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(A));
    if (A->nvals == 0) {
        return GrB_SUCCESS;
    }

    TRY(matrix_clear_format(A, optimizations));
    return GrB_SUCCESS;
}

// duplicate methods

GrB_Info matrix_dup(Matrix *output, Matrix *input) {
    if (output == input) {
        return GrB_SUCCESS;
    }

    TRY(GrB_Matrix_apply(output->base, GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL, input->base,
                         GrB_NULL));
    TRY(CFL_matrix_update(output));

    return GrB_SUCCESS;
}

GrB_Info matrix_dup_format(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        TRY(matrix_dup(output, input));
        return GrB_SUCCESS;
    }

    if (!output->is_both) {
        TRY(CFL_matrix_update(output));
        TRY(CFL_matrix_update(input));
        Matrix *larger = output->nvals > input->nvals ? output : input;

        TRY(matrix_to_format(output, larger->format, false));
        TRY(matrix_to_format(input, larger->format, false));

        TRY(matrix_dup(output, input));

        return GrB_SUCCESS;
    }

    TRY(matrix_to_format(output, GrB_ROWMAJOR, false));
    TRY(matrix_dup(output, input));
    TRY(matrix_to_format(output, GrB_COLMAJOR, false));
    TRY(matrix_dup(output, input));

    return GrB_SUCCESS;
}

GrB_Info matrix_dup_empty(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        TRY(matrix_dup_format(output, input, optimizations));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(input));
    if (input->nvals == 0) {
        TRY(matrix_clear_empty(output, optimizations));
        return GrB_SUCCESS;
    }

    TRY(matrix_dup_format(output, input, optimizations));
    return GrB_SUCCESS;
}

GrB_Info matrix_dup_lazy(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        TRY(matrix_dup_empty(output, input, optimizations));
        return GrB_SUCCESS;
    }

    if (!input->is_lazy) {
        TRY(matrix_dup_empty(output, input, optimizations));
        return GrB_SUCCESS;
    }

    for (size_t i = 0; i < input->base_matrices_count; i++) {
        TRY(matrix_dup_empty(output->base_matrices[i], input->base_matrices[i],
                             optimizations));
    }
    output->base_matrices_count = input->base_matrices_count;

    return GrB_SUCCESS;
}

GrB_Info block_matrix_hyper_rotate_i(Matrix *matrix, enum CFL_Matrix_block format);

GrB_Info matrix_dup_block(Matrix *output, Matrix *input, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        TRY(matrix_dup_lazy(output, input, optimizations));
        return GrB_SUCCESS;
    }

    if (output->block_type == CELL && input->block_type == CELL) {
        TRY(matrix_dup_lazy(output, input, optimizations));
        return GrB_SUCCESS;
    }

    TRY(block_matrix_hyper_rotate_i(input, output->block_type));
    TRY(matrix_dup_lazy(output, input, optimizations));
    return GrB_SUCCESS;
}

// block optimization specific methods

GrB_Info block_matrix_hyper_rotate_i(Matrix *matrix, enum CFL_Matrix_block format) {
    if (matrix->is_lazy) {
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            TRY(block_matrix_hyper_rotate_i(matrix->base_matrices[i], format));
        }

        TRY(CFL_matrix_update(matrix));
        return GrB_SUCCESS;
    }

    if (matrix->block_type == CELL) {
        return GrB_SUCCESS;
    }

    if (matrix->block_type == format) {
        return GrB_SUCCESS;
    }

    GrB_Scalar scalar_true;
    TRY(GrB_Scalar_new(&scalar_true, GrB_BOOL));
    TRY(GrB_Scalar_setElement_BOOL(scalar_true, true));

    TRY(CFL_matrix_update(matrix));

    GrB_Index *nrows;
    GrB_Index *ncols;
    TRY(LAGraph_Calloc((void **)&nrows, matrix->nvals, sizeof(GrB_Index), NULL));
    TRY(LAGraph_Calloc((void **)&ncols, matrix->nvals, sizeof(GrB_Index), NULL));

    TRY(GrB_Matrix_extractTuples_BOOL(nrows, ncols, NULL, &matrix->nvals, matrix->base));

    if (matrix->block_type == VEC_VERT) {
        for (size_t i = 0; i < matrix->nvals; i++) {
            ncols[i] = ncols[i] + nrows[i] / matrix->ncols * matrix->ncols;
            nrows[i] = nrows[i] % matrix->ncols;
        }
    } else {
        for (size_t i = 0; i < matrix->nvals; i++) {
            nrows[i] = nrows[i] + ncols[i] / matrix->nrows * matrix->nrows;
            ncols[i] = ncols[i] % matrix->nrows;
        }
    }

    GrB_Matrix new;
    TRY(GrB_Matrix_new(&new, GrB_BOOL, matrix->ncols, matrix->nrows));
    if (matrix->is_both) {
        TRY(GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals));
        TRY(GrB_Matrix_free(&matrix->base_row));
        TRY(GrB_Matrix_free(&matrix->base_col));
        matrix->base = NULL;
        matrix->base = new;
        matrix->base_row = new;
        TRY(TO_ROW(matrix->base));
        TRY(GrB_Matrix_new(&new, GrB_BOOL, matrix->ncols, matrix->nrows));
        TRY(TO_COL(new));
        TRY(GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals));
        matrix->base_col = new;

        int format;
        TRY(GrB_get(matrix->base_row, &format, GrB_STORAGE_ORIENTATION_HINT));
        if (format != GrB_ROWMAJOR) {
            fprintf(stderr, "WRONG FORMAT\n");
            exit(-1);
        }
        TRY(GrB_get(matrix->base_col, &format, GrB_STORAGE_ORIENTATION_HINT));
        if (format != GrB_COLMAJOR) {
            fprintf(stderr, "WRONG FORMAT\n");
            exit(-1);
        }
        matrix->format = GrB_ROWMAJOR;
    } else {
        TRY(GrB_Matrix_free(&matrix->base));
        matrix->base = new;
        if (matrix->base_row != NULL) {
            TRY(GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals));
            matrix->base_row = new;
        } else {
            TRY(TO_COL(new));
            TRY(GxB_Matrix_build_Scalar(new, nrows, ncols, scalar_true, matrix->nvals));
            matrix->base_col = new;
        }
        matrix->format = GrB_ROWMAJOR;
    }
    TRY(CFL_matrix_update(matrix));
    // TRY(CFL_matrix_free(matrix_p));
    // TRY(CFL_matrix_from_base(matrix_p, new));

    TRY(LAGraph_Free((void **)&nrows, NULL));
    TRY(LAGraph_Free((void **)&ncols, NULL));
    TRY(GrB_free(&scalar_true));

    return GrB_SUCCESS;
}

GrB_Info block_matrix_to_diag(Matrix *diag, Matrix *input) {
    if (input->block_type == CELL) {
        return GrB_SUCCESS;
    }

    GrB_Scalar scalar_true;
    TRY(GrB_Scalar_new(&scalar_true, GrB_BOOL));
    TRY(GrB_Scalar_setElement_BOOL(scalar_true, true));

    TRY(CFL_matrix_update(input));

    GrB_Index *rows;
    GrB_Index *cols;
    TRY(LAGraph_Calloc((void **)&rows, input->nvals, sizeof(GrB_Index), NULL));
    TRY(LAGraph_Calloc((void **)&cols, input->nvals, sizeof(GrB_Index), NULL));
    TRY(GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &input->nvals, input->base));

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

    TRY(GxB_Matrix_build_Scalar(diag->base, rows, cols, scalar_true, input->nvals));
    TRY(CFL_matrix_update(diag));

    TRY(LAGraph_Free((void **)&rows, NULL));
    TRY(LAGraph_Free((void **)&cols, NULL));
    TRY(GrB_free(&scalar_true));

    return GrB_SUCCESS;
}

GrB_Info block_matrix_reduce(Matrix *matrix, Matrix *input, int8_t optimizations) {
    if (input->block_type == CELL) {
        TRY(matrix_dup_block(matrix, input, optimizations));
        return GrB_SUCCESS;
    }

    GrB_Scalar scalar_true;
    TRY(GrB_Scalar_new(&scalar_true, GrB_BOOL));
    TRY(GrB_Scalar_setElement_BOOL(scalar_true, true));

    TRY(CFL_matrix_update(input));

    GrB_Index *rows;
    GrB_Index *cols;
    TRY(LAGraph_Calloc((void **)&rows, input->nvals, sizeof(GrB_Index), NULL));
    TRY(LAGraph_Calloc((void **)&cols, input->nvals, sizeof(GrB_Index), NULL));
    TRY(GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &input->nvals, input->base));

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

    TRY(GxB_Matrix_build_Scalar(matrix->base, rows, cols, scalar_true, input->nvals));
    TRY(CFL_matrix_update(matrix));

    TRY(LAGraph_Free((void **)&rows, NULL));
    TRY(LAGraph_Free((void **)&cols, NULL));
    TRY(GrB_free(&scalar_true));

    return GrB_SUCCESS;
}

GrB_Info block_matrix_repeat_into_vector(Matrix *matrix, Matrix *input,
                                         GrB_Index block_count) {
    GrB_Matrix *tiles;
    TRY(LAGraph_Calloc((void **)&tiles, block_count, sizeof(GrB_Matrix), NULL));
    for (size_t i = 0; i < block_count; i++) {
        tiles[i] = input->base;
    }

    if (matrix->block_type == VEC_VERT) {
        TRY(GxB_Matrix_concat(matrix->base, tiles, block_count, 1, GrB_NULL));
    } else {
        TRY(GxB_Matrix_concat(matrix->base, tiles, 1, block_count, GrB_NULL));
    }
    TRY(CFL_matrix_update(matrix));
    TRY(LAGraph_Free((void **)&tiles, NULL));

    return GrB_SUCCESS;
}

// lazy optimization specific methods

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations);

GrB_Info matrix_sort_lazy(Matrix *A, bool reverse) {
    for (size_t i = 0; i < A->base_matrices_count; i++) {
        for (size_t j = i + 1; j < A->base_matrices_count; j++) {
            Matrix *first = reverse ? A->base_matrices[i] : A->base_matrices[j];
            Matrix *second = reverse ? A->base_matrices[j] : A->base_matrices[i];
            TRY(CFL_matrix_update(first));
            TRY(CFL_matrix_update(second));

            if (first->nvals < second->nvals) {
                Matrix *temp = A->base_matrices[i];
                A->base_matrices[i] = A->base_matrices[j];
                A->base_matrices[j] = temp;
            }
        }
    }

    return GrB_SUCCESS;
}

GrB_Info CFL_matrix_to_base(Matrix **matrix_p, Matrix *input, int8_t optimizations) {
    if (!input->is_lazy) {
        Matrix *result;
        TRY(CFL_matrix_create(&result, input->nrows, input->ncols));
        TRY(CFL_dup(result, input, optimizations));
        *matrix_p = result;

        return GrB_SUCCESS;
    }

    Matrix *matrix;
    TRY(CFL_matrix_create_lazy(&matrix, input->nrows, input->ncols));
    TRY(CFL_matrix_free(&matrix->base_matrices[0]));
    matrix->base_matrices_count = 0;

    for (size_t i = 0; i < input->base_matrices_count; i++) {
        Matrix *base;
        TRY(CFL_matrix_create(&base, input->nrows, input->ncols));
        matrix->base_matrices[i] = base;
        matrix->base_matrices_count++;
    }

    TRY(CFL_matrix_update(matrix));
    TRY(CFL_dup(matrix, input, optimizations));
    TRY(matrix_sort_lazy(matrix, false));

    Matrix *acc;
    TRY(CFL_matrix_create(&acc, input->nrows, input->ncols));
    for (size_t j = 0; j < matrix->base_matrices_count; j++) {
        TRY(CFL_wise(acc, acc, matrix->base_matrices[j], false, optimizations));
    }

    TRY(CFL_matrix_free(&matrix));

    *matrix_p = acc;

    return GrB_SUCCESS;
}

GrB_Info matrix_combine_lazy(Matrix *A, size_t threshold, int8_t optimizations) {
    Matrix **new_matrices;
    TRY(LAGraph_Calloc((void **)&new_matrices, 50, sizeof(Matrix *), NULL));

    size_t new_size = 0;

    TRY(matrix_sort_lazy(A, false));

    for (size_t i = 0; i < A->base_matrices_count; i++) {
        TRY(CFL_matrix_update(A->base_matrices[i]));
        if (new_size == 0 || A->base_matrices[i]->nvals > threshold) {
            new_matrices[new_size++] = A->base_matrices[i];
            continue;
        }

        TRY(matrix_wise_empty(new_matrices[new_size - 1], new_matrices[new_size - 1],
                              A->base_matrices[i], false, optimizations));
        TRY(CFL_matrix_free(&A->base_matrices[i]));
    }

    TRY(LAGraph_Free((void **)&A->base_matrices, NULL));
    A->base_matrices = new_matrices;
    A->base_matrices_count = new_size;
    TRY(CFL_matrix_update(A));

    return GrB_SUCCESS;
}

// create and update methods
GrB_Info CFL_matrix_update(Matrix *matrix) {
    if (!matrix->is_lazy) {
        TRY(GrB_Matrix_nvals(&matrix->nvals, matrix->base));
    } else {
        size_t new_nnz = 0;
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            TRY(GrB_Matrix_nvals(&(matrix->base_matrices[i]->nvals),
                                 matrix->base_matrices[i]->base));
            new_nnz += matrix->base_matrices[i]->nvals;
        }

        matrix->nvals = new_nnz;
    }

    if (!matrix->is_lazy) {
        TRY(GrB_Matrix_nrows(&matrix->nrows, matrix->base));
        TRY(GrB_Matrix_ncols(&matrix->ncols, matrix->base));
    } else {
        TRY(GrB_Matrix_nrows(&matrix->nrows, matrix->base_matrices[0]->base));
        TRY(GrB_Matrix_ncols(&matrix->ncols, matrix->base_matrices[0]->base));
    }

    if (matrix->nrows == matrix->ncols)
        matrix->block_type = CELL;
    else
        matrix->block_type = matrix->nrows > matrix->ncols ? VEC_VERT : VEC_HORIZ;

#ifdef BENCH_CFL_REACHBILITY
    if (matrix->is_lazy) {
        for (size_t i = 0; i < matrix->base_matrices_count; i++) {
            IS_ISO(matrix->base_matrices[i]->base, "lazy matrix base");
            IS_ROW(matrix->base_matrices[i]->base_row, "");
            IS_COL(matrix->base_matrices[i]->base_col, "");
        }
    } else {
        IS_ISO(matrix->base, "");
        IS_ROW(matrix->base_row, "");
        IS_COL(matrix->base_col, "");
    }
#endif

    return GrB_SUCCESS;
}

// Create optimizied matrix from base GrB_Matrix
// Input: base
// Output: matrix
GrB_Info CFL_matrix_from_base(Matrix **matrix, GrB_Matrix base) {
    Matrix *result;
    TRY(LAGraph_Calloc((void **)&result, 1, sizeof(Matrix), NULL));

    result->base = base;
    result->nvals = 0; // We will get actual info in update function
    result->nrows = 0;
    result->ncols = 0;

    // Format optimization fields
    result->base_row = base;
    result->base_col = NULL;
    result->format = GrB_ROWMAJOR;
    result->is_both = false;

    // Lazy addition optimization fields
    result->is_lazy = false;
    result->base_matrices = NULL;
    result->base_matrices_count = 0;

    // Block optimization fields
    result->block_type = CELL;

    TRY(CFL_matrix_update(result));
    *matrix = result;

    return GrB_SUCCESS;
}

// Create optimizied lazy matrix from base GrB_Matrix
// Input: base
// Output: matrix_p that holds new allocated matrix.
//
// Base matrix will be in array of lazy matrix
GrB_Info CFL_matrix_from_base_lazy(Matrix **matrix_p, GrB_Matrix base) {
    Matrix *result;
    TRY(CFL_matrix_from_base(&result, base));

    TRY(CFL_matrix_from_base(matrix_p, base));

    Matrix *matrix = *matrix_p;
    matrix->is_lazy = true;
    TRY(LAGraph_Calloc(
        (void **)&matrix->base_matrices, 40, sizeof(CFL_Matrix *),
        NULL)); // this is enough for this centry i guess, because 40th matrices must
                // have 10^40 nvals for being putten in base_matrices, this is 2^132
    matrix->base_matrices[0] = result;
    matrix->base_matrices_count = 1;
    matrix->base = NULL;
    TRY(CFL_matrix_update(matrix));

    return GrB_SUCCESS;
}

GrB_Info CFL_matrix_create(Matrix **matrix, GrB_Index nrows, GrB_Index ncols) {
    GrB_Matrix _result;
    TRY(GrB_Matrix_new(&_result, GrB_BOOL, nrows, ncols));
    TRY(CFL_matrix_from_base(matrix, _result));

    return GrB_SUCCESS;
}

GrB_Info CFL_matrix_create_lazy(Matrix **matrix, GrB_Index nrows, GrB_Index ncols) {
    GrB_Matrix _result;
    TRY(GrB_Matrix_new(&_result, GrB_BOOL, nrows, ncols));
    TRY(CFL_matrix_from_base_lazy(matrix, _result))

    return GrB_SUCCESS;
}

GrB_Info CFL_matrix_free(Matrix **matrix_p) {
    if (*matrix_p == NULL) {
        return GrB_SUCCESS;
    }

    Matrix *matrix = *matrix_p;
    if (matrix->is_both) {
        TRY(GrB_Matrix_free(&(matrix->base_col)));
        TRY(GrB_Matrix_free(&(matrix->base_row)));
    } else {
        TRY(GrB_Matrix_free(&(matrix->base)));
    }

    for (size_t i = 0; i < matrix->base_matrices_count; i++) {
        TRY(CFL_matrix_free(&(matrix->base_matrices[i])));
    }
    TRY(LAGraph_Free((void **)&matrix->base_matrices, NULL));
    TRY(LAGraph_Free((void **)matrix_p, NULL));

    return GrB_SUCCESS;

    GrB_Matrix_free(&matrix->base);
}

// mxm operations

GrB_Info matrix_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum,
                    bool swap) {
    Matrix *left = swap ? second : first;
    Matrix *right = swap ? first : second;

    // CFL_matrix_update(first);
    // CFL_matrix_update(second);
    // CFL_matrix_update(output);
    TRY(GrB_mxm(output->base, GrB_NULL, accum ? GxB_ANY_BOOL : GrB_NULL,
                GxB_ANY_PAIR_BOOL, left->base, right->base, GrB_NULL));
    // IS_ISO(output->base, "MXM output");
    TRY(CFL_matrix_update(output));

    return GrB_SUCCESS;
}

GrB_Info matrix_mxm_format(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        TRY(matrix_mxm(output, first, second, accum, swap))
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));
    GrB_Index left_nvals = swap ? second->nvals : first->nvals;
    GrB_Index right_nvals = swap ? first->nvals : second->nvals;

    int32_t desired_orientation = left_nvals < right_nvals ? GrB_ROWMAJOR : GrB_COLMAJOR;

    if (!first->is_both && first->format != desired_orientation &&
        !(first->nvals > second->nvals / 3.0)) {
        TRY(matrix_mxm(output, first, second, accum, swap))
        return GrB_SUCCESS;
    }

    TRY(matrix_to_format(first, desired_orientation, true));
    TRY(matrix_to_format(second, desired_orientation, false));
    TRY(matrix_to_format(output, desired_orientation, false));
    TRY(matrix_mxm(output, first, second, accum, swap));

    return GrB_SUCCESS;
}

GrB_Info matrix_mxm_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        TRY(matrix_mxm_format(output, first, second, accum, swap, optimizations))
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));

    if (first->nvals == 0 || second->nvals == 0) {
        if (accum) {
            return GrB_SUCCESS;
        }

        TRY(CFL_matrix_update(output));
        if (output->nvals == 0) {
            return GrB_SUCCESS;
        }

        TRY(matrix_clear_empty(output, optimizations));
        return GrB_SUCCESS;
    }

    TRY(matrix_mxm_format(output, first, second, accum, swap, optimizations));
    return GrB_SUCCESS;
}

GrB_Info matrix_mxm_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum,
                         bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        TRY(matrix_mxm_empty(output, first, second, accum, swap, optimizations));
        return GrB_SUCCESS;
    }

    if (!first->is_lazy) {
        TRY(matrix_mxm_empty(output, first, second, accum, swap, optimizations));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(second));

    TRY(matrix_combine_lazy(first, second->nvals, optimizations));
    TRY(matrix_sort_lazy(first, false));

    GrB_Matrix *accs;
    TRY(LAGraph_Calloc((void **)&accs, first->base_matrices_count, sizeof(GrB_Matrix),
                       NULL));
    Matrix **acc_matrices;
    TRY(LAGraph_Calloc((void **)&acc_matrices, first->base_matrices_count,
                       sizeof(Matrix *), NULL));
    for (size_t i = 0; i < first->base_matrices_count; i++) {
        TRY(GrB_Matrix_new(&accs[i], GrB_BOOL, swap ? second->nrows : first->nrows,
                           swap ? first->ncols : second->ncols));
        TRY(CFL_matrix_from_base(&acc_matrices[i], accs[i]))
    }

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        TRY(matrix_mxm_empty(acc_matrices[i], first->base_matrices[i], second, false,
                             swap, optimizations));
    }

    GrB_Matrix acc;
    GrB_Matrix_new(&acc, GrB_BOOL, swap ? second->nrows : first->nrows,
                   swap ? first->ncols : second->ncols);
    Matrix *acc_matrix;
    CFL_matrix_from_base(&acc_matrix, acc);

    for (size_t i = 0; i < first->base_matrices_count; i++) {
        TRY(matrix_wise_empty(acc_matrix, acc_matrix, acc_matrices[i], false,
                              optimizations));
        TRY(CFL_matrix_free(&acc_matrices[i]));
    }
    LAGraph_Free((void **)&accs, NULL);
    LAGraph_Free((void **)&acc_matrices, NULL);

    if (accum) {
        TRY(matrix_wise_empty(output, output, acc_matrix, false, optimizations));
    } else {
        TRY(matrix_dup_block(output, acc_matrix, optimizations));
    }

    TRY(CFL_matrix_free(&acc_matrix));

    return GrB_SUCCESS;
}

GrB_Info matrix_wise_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations);

GrB_Info matrix_mxm_block(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          bool swap, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        TRY(matrix_mxm_lazy(output, first, second, accum, swap, optimizations));
        return GrB_SUCCESS;
    }

    if (first->block_type == CELL && second->block_type == CELL) {
        TRY(matrix_mxm_lazy(output, first, second, accum, swap, optimizations));
        return GrB_SUCCESS;
    }

    if (first->block_type == CELL) {
        TRY(block_matrix_hyper_rotate_i(second, swap ? VEC_VERT : VEC_HORIZ));
        TRY(block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ));

        Matrix *temp;
        TRY(CFL_matrix_create(&temp, swap ? second->nrows : first->nrows,
                              swap ? first->ncols : second->ncols));
        TRY(matrix_mxm_lazy(temp, first, second, accum, swap, optimizations));
        TRY(matrix_wise_block(output, output, temp, false, optimizations));
        TRY(CFL_matrix_free(&temp));

        return GrB_SUCCESS;
    }

    if (second->block_type == CELL) {
        TRY(block_matrix_hyper_rotate_i(first, swap ? VEC_HORIZ : VEC_VERT));
        TRY(block_matrix_hyper_rotate_i(output, swap ? VEC_HORIZ : VEC_VERT));

        Matrix *temp;
        TRY(CFL_matrix_create(&temp, swap ? second->nrows : first->nrows,
                              swap ? first->ncols : first->ncols));
        TRY(matrix_mxm_lazy(temp, first, second, accum, swap, optimizations));
        TRY(matrix_wise_block(output, output, temp, false, optimizations));
        TRY(CFL_matrix_free(&temp));

        return GrB_SUCCESS;
    }

    GrB_Index size = first->nrows > first->ncols ? first->nrows : first->ncols;

    Matrix *diag;
    TRY(CFL_matrix_create(&diag, size, size));
    TRY(block_matrix_to_diag(diag, second));
    TRY(CFL_matrix_update(diag));

    TRY(block_matrix_hyper_rotate_i(first, swap ? VEC_VERT : VEC_HORIZ));
    TRY(block_matrix_hyper_rotate_i(output, swap ? VEC_VERT : VEC_HORIZ));

    Matrix *temp;
    TRY(CFL_matrix_create(&temp, first->nrows, first->ncols));
    TRY(matrix_mxm_lazy(temp, first, diag, false, swap, optimizations));
    TRY(matrix_wise_block(output, output, temp, false, optimizations));

    TRY(CFL_matrix_free(&temp));
    TRY(CFL_matrix_free(&diag));

    return GrB_SUCCESS;
}

// wise operations

GrB_Info matrix_wise(Matrix *output, Matrix *first, Matrix *second, bool accum) {
    GrB_BinaryOp accum_op = accum ? GxB_ANY_BOOL : GrB_NULL;
    if (output == first)
        accum_op = GrB_NULL;

    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));
    TRY(CFL_matrix_update(output));

    TRY(GrB_Matrix_eWiseAdd_BinaryOp(output->base, GrB_NULL, accum_op, GxB_ANY_BOOL,
                                     first->base, second->base, GrB_NULL));

    TRY(CFL_matrix_update(output));

    return GrB_SUCCESS;
}

GrB_Info matrix_wise_format(Matrix *output, Matrix *first, Matrix *second, bool accum,
                            int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        TRY(matrix_wise(output, first, second, accum))
        return GrB_SUCCESS;
    }

    if (!output->is_both) {
        TRY(CFL_matrix_update(output));
        TRY(CFL_matrix_update(first));
        TRY(CFL_matrix_update(second));

        Matrix *larger = output->nvals > first->nvals ? output : first;
        larger = larger->nvals > second->nvals ? larger : second;

        TRY(matrix_to_format(output, larger->format, false));
        TRY(matrix_to_format(first, larger->format, false));
        TRY(matrix_to_format(second, larger->format, false));

        TRY(matrix_wise(output, first, second, accum));

        return GrB_SUCCESS;
    }

    TRY(matrix_to_format(output, GrB_ROWMAJOR, false));
    TRY(matrix_to_format(first, output->format, false));
    TRY(matrix_to_format(second, output->format, false));

    TRY(matrix_wise(output, first, second, accum));

    TRY(matrix_to_format(output, GrB_COLMAJOR, false));
    TRY(matrix_to_format(first, output->format, false));
    TRY(matrix_to_format(second, output->format, false));

    TRY(matrix_wise(output, first, second, accum));

    return GrB_SUCCESS;
}

GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second, bool accum,
                           int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        TRY(matrix_wise_format(output, first, second, accum, optimizations));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));
    TRY(CFL_matrix_update(output));

    if (output == first) {
        if (first->nvals == 0) {
            TRY(matrix_dup_empty(first, second, optimizations));
            return GrB_SUCCESS;
        }

        if (second->nvals == 0) {
            return GrB_SUCCESS;
        }

        TRY(matrix_wise_format(output, first, second, accum, optimizations))
        return GrB_SUCCESS;
    }

    if (first->nvals == 0 && second->nvals == 0) {
        if (accum || output->nvals == 0) {
            return GrB_SUCCESS;
        }

        TRY(matrix_clear_empty(output, optimizations));
        return GrB_SUCCESS;
    }

    if (first->nvals == 0) {
        if (accum) {
            TRY(matrix_wise_empty(output, output, second, false, optimizations));
            return GrB_SUCCESS;
        }

        TRY(matrix_dup_empty(output, second, optimizations));
        return GrB_SUCCESS;
    }

    if (second->nvals == 0) {
        if (accum) {
            TRY(matrix_wise_empty(output, output, first, false, optimizations));
            return GrB_SUCCESS;
        }

        TRY(matrix_dup_empty(output, first, optimizations));
        return GrB_SUCCESS;
    }

    TRY(matrix_wise_format(output, first, second, accum, optimizations));
    return GrB_SUCCESS;
}

GrB_Info matrix_wise_lazy(Matrix *output, Matrix *first, Matrix *second, bool accum,
                          int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        TRY(matrix_wise_empty(output, first, second, accum, optimizations))
        return GrB_SUCCESS;
    }

    if (!first->is_lazy && !second->is_lazy) {
        TRY(matrix_wise_empty(output, first, second, accum, optimizations));
        return GrB_SUCCESS;
    }

    if (!first->is_lazy && second->is_lazy) {
        for (size_t i = 0; i < second->base_matrices_count; i++) {
            TRY(matrix_wise_empty(output, first, second->base_matrices[i], true,
                                  optimizations));
        }

        return GrB_SUCCESS;
    }

    Matrix *other;
    TRY(CFL_matrix_create(&other, output->nrows, output->ncols));
    TRY(matrix_dup_empty(other, second, optimizations));

    size_t other_nvals = other->nvals >= 10 ? other->nvals : 10;

    while (true) {
        bool found = false;

        for (size_t i = 0; i < first->base_matrices_count; i++) {
            TRY(CFL_matrix_update(first->base_matrices[i]));
            size_t self_nvals = first->base_matrices[i]->nvals >= 10
                                    ? first->base_matrices[i]->nvals
                                    : 10;

            if (other_nvals / 10 <= self_nvals && self_nvals <= other_nvals * 10) {
                TRY(matrix_wise_empty(other, other, first->base_matrices[i], accum,
                                      optimizations));
                TRY(CFL_matrix_free(&first->base_matrices[i]));
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

    TRY(matrix_sort_lazy(first, false));
    TRY(CFL_matrix_update(output));

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
        TRY(matrix_wise_lazy(output, first, second, accum, optimizations));
        return GrB_SUCCESS;
    }

    if (output != first) {
        // fprintf(stderr, "Matrix wise currently support only iadd operation");
        return GrB_INVALID_VALUE;
    }

    if (first->block_type == CELL && second->block_type == CELL) {
        TRY(matrix_wise_lazy(output, first, second, accum, optimizations));
        return GrB_SUCCESS;
    }

    // second is vector
    if (first->block_type == CELL) {
        Matrix *temp_reduced;
        TRY(CFL_matrix_create(&temp_reduced, first->nrows, first->ncols));
        TRY(block_matrix_reduce(temp_reduced, second, optimizations));

        TRY(matrix_wise_lazy(output, first, temp_reduced, accum, optimizations));
        TRY(CFL_matrix_free(&temp_reduced));

        return GrB_SUCCESS;
    }

    // first is vector
    if (second->block_type == CELL) {
        // LG_SET_BURBLE(true);
        Matrix *temp_vector;
        TRY(CFL_matrix_create(&temp_vector, first->nrows, first->ncols));

        GrB_Index block_count = 0;
        if (first->nrows > first->ncols) {
            block_count = first->nrows / first->ncols;
        } else {
            block_count = first->ncols / first->nrows;
        }
        TRY(block_matrix_repeat_into_vector(temp_vector, second, block_count));
        TRY(block_matrix_hyper_rotate_i(temp_vector, first->block_type));
        TRY(matrix_wise_lazy(output, first, temp_vector, accum, optimizations));

        TRY(CFL_matrix_free(&temp_vector));

        return GrB_SUCCESS;
    }

    // both are vector
    TRY(block_matrix_hyper_rotate_i(second, first->block_type));
    TRY(matrix_wise_lazy(output, first, second, accum, optimizations));

    return GrB_SUCCESS;
}

// rsub methods

GrB_Info matrix_rsub(Matrix *output, Matrix *mask) {
    TRY(GrB_eWiseAdd(output->base, mask->base, GrB_NULL, GxB_ANY_BOOL, output->base,
                     output->base, GrB_DESC_RSC));

    TRY(CFL_matrix_update(output));

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_format(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_FORMAT)) {
        TRY(matrix_rsub(output, mask));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(output));
    TRY(CFL_matrix_update(mask));

    Matrix *larger_matrix = output->nvals > mask->nvals ? output : mask;
    TRY(matrix_to_format(output, larger_matrix->format, false));
    TRY(matrix_to_format(mask, larger_matrix->format, false));

    TRY(matrix_rsub(output, mask));

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_empty(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_EMPTY)) {
        TRY(matrix_rsub_format(output, mask, optimizations));
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(mask));
    TRY(CFL_matrix_update(output));
    if (mask->nvals == 0 || output->nvals == 0) {
        return GrB_SUCCESS;
    }

    TRY(matrix_rsub_format(output, mask, optimizations));

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_lazy(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_LAZY)) {
        TRY(matrix_rsub_empty(output, mask, optimizations));
        return GrB_SUCCESS;
    }

    if (!mask->is_lazy) {
        TRY(matrix_rsub_empty(output, mask, optimizations))
        return GrB_SUCCESS;
    }

    TRY(CFL_matrix_update(output));
    TRY(matrix_combine_lazy(mask, output->nvals, optimizations));
    TRY(matrix_sort_lazy(mask, true));

    for (size_t i = 0; i < mask->base_matrices_count; i++) {
        TRY(matrix_rsub_empty(output, mask->base_matrices[i], optimizations));
    }

    return GrB_SUCCESS;
}

GrB_Info matrix_rsub_block(Matrix *output, Matrix *mask, int8_t optimizations) {
    if (!(optimizations & OPT_BLOCK)) {
        TRY(matrix_rsub_lazy(output, mask, optimizations))
        return GrB_SUCCESS;
    }

    if ((output->block_type == CELL && mask->block_type != CELL) ||
        (output->block_type != CELL && mask->block_type == CELL)) {
        // fprintf(stderr, "Don't support rsub operation between cell and vector");
        return GrB_INVALID_VALUE;
    }

    if (output->block_type == CELL) {
        TRY(matrix_rsub_lazy(output, mask, optimizations));
        return GrB_SUCCESS;
    }

    TRY(block_matrix_hyper_rotate_i(output, mask->block_type));
    TRY(matrix_rsub_lazy(output, mask, optimizations));

    return GrB_SUCCESS;
}

// utility methods

GrB_Info matrix_print_lazy(Matrix *A, int8_t optimizations) {
    // return;
    GxB_Print_Level pr = 1;

    if (!A->is_lazy) {
        TRY(CFL_matrix_update(A));
        TRY(GxB_print(A->base, pr));
        // printf("nnz: %ld\n", A->nvals);
        return GrB_SUCCESS;
    }

    if (A->base_matrices_count == 1) {
        TRY(CFL_matrix_update(A));
        // printf("nnz: %ld\n", A->nvals);
        TRY(GxB_print(A->base_matrices[0]->base, pr));
        return GrB_SUCCESS;
    }

    Matrix *temp;
    TRY(CFL_matrix_create(&temp, A->nrows, A->ncols));
    for (size_t i = 0; i < A->base_matrices_count; i++) {
        TRY(matrix_wise_empty(temp, temp, A->base_matrices[i], false, optimizations));
    }

    A = temp;
    TRY(GxB_print(A->base, pr));
    TRY(CFL_matrix_update(A));
    // printf("nnz: %ld\n", A->nvals);
    TRY(CFL_matrix_free(&temp));

    return GrB_SUCCESS;
}

GrB_Info print_graph_info(Matrix **matrices, size_t count) {
    GrB_Index nnz = 0;

    for (size_t i = 0; i < count; i++) {
        Matrix *A = matrices[i];
        TRY(CFL_matrix_update(A));
        nnz += A->nvals;
    }

    printf("NNZ: %ld\n", nnz);
    return GrB_SUCCESS;
}

// order of optimizations: block -> lazy -> empty -> format

GrB_Info CFL_mxm(Matrix *output, Matrix *first, Matrix *second, bool accum, bool swap,
                 int8_t optimizations) {
    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));
    TRY(CFL_matrix_update(output));

    TRY(matrix_mxm_block(output, first, second, accum, swap, optimizations));

    TRY(CFL_matrix_update(output));
    return GrB_SUCCESS;
}

GrB_Info CFL_wise(Matrix *output, Matrix *first, Matrix *second, bool accum,
                  int8_t optimizations) {
    TRY(CFL_matrix_update(first));
    TRY(CFL_matrix_update(second));
    TRY(CFL_matrix_update(output));

    TRY(matrix_wise_block(output, first, second, accum, optimizations));

    TRY(CFL_matrix_update(output));

    return GrB_SUCCESS;
}

GrB_Info CFL_rsub(Matrix *output, Matrix *mask, int8_t optimizations) {
    TRY(CFL_matrix_update(mask));
    TRY(CFL_matrix_update(output));

    TRY(matrix_rsub_block(output, mask, optimizations))

    TRY(CFL_matrix_update(output));
    return GrB_SUCCESS;
}

GrB_Info CFL_dup(Matrix *output, Matrix *input, int8_t optimizations) {
    TRY(CFL_matrix_update(output));
    TRY(CFL_matrix_update(input));

    TRY(matrix_dup_block(output, input, optimizations))

    TRY(CFL_matrix_update(output));
    return GrB_SUCCESS;
}

GrB_Info CFL_clear(Matrix *A, int8_t optimizations) {
    TRY(CFL_matrix_update(A));

    TRY(matrix_clear_empty(A, optimizations))

    TRY(CFL_matrix_update(A));
    return GrB_SUCCESS;
}