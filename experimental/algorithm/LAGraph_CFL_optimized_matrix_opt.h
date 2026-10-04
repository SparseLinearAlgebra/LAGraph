//------------------------------------------------------------------------------
// LAGraph_CFL_optimized_matrix_opt.h: Header with operations for
// Optimized Context-Free Language Reachability Matrix-Based Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vlasenco Daniel, Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

// Code is an implementation of optimized matrix operations for CFL_reachability
// algorithms, described in the following paper:
// * Ilia Muravev, "Optimization of the Context-Free Language Reachability Matrix-Based
// Algorithm" and based on the python implementation from:
// https://github.com/FormalLanguageConstrainedPathQuerying/CFPQ_PyAlgo/tree/murav/optimize-matrix

#include <GraphBLAS.h>

// Several square matrices can be concatenated into a horizontal or vertical vector of
// matrices for fast multiplication
//
// [CELL]      - the matrix is square: n x n
// [VEC_HORIZ] - the matrix is a horizontal block vector: n x (n * k),
//               where k is the number of concatenated matrices
// [VEC_VERT]  - the matrix is a vertical block vector: (n * k) x n,
//               where k is the number of concatenated matrices
enum CFL_Matrix_block { CELL, VEC_HORIZ, VEC_VERT };

// Matrix wrapper for CFL reachability algorithms with optional optimizations
typedef struct CFL_Matrix {
    GrB_Matrix base;      // Underlying GrB_Matrix
    int8_t optimizations; // Optimizations flags
    // Fields of base matrix
    GrB_Index nvals;
    GrB_Index nrows;
    GrB_Index ncols;
    // Fields of format optimization
    GrB_Matrix base_row;
    GrB_Matrix base_col;
    int32_t format;
    bool is_both;
    // Fields of lazy addition optimization
    struct CFL_Matrix **base_matrices;
    size_t base_matrices_count;
    bool is_lazy;
    // Fields of block optimization
    enum CFL_Matrix_block block_type;
} CFL_Matrix;

// Creates a CFL matrix from intialized GrB_Matrix
//
// Creates a CFL_Matrix wrapping an already-initialized GrB_Matrix
// On success, *matrix is set to a newly allocated CFL_Matrix
GrB_Info CFL_matrix_from_base(CFL_Matrix **matrix, GrB_Matrix base);

// Creates a lazy CFL matrix from intialized GrB_Matrix
//
// Creates a CFL_Matrix wrapping an already-initialized GrB_Matrix
// On success, *matrix is set to a newly allocated CFL_Matrix
GrB_Info CFL_matrix_from_base_lazy(CFL_Matrix **matrix, GrB_Matrix base);

// Creates an empty CFL_Matrix of the given dimensions
//
// The underlying GrB_Matrix is allocated internally
// On success, *matrix is set to a newly allocated CFL_Matrix
GrB_Info CFL_matrix_create(CFL_Matrix **matrix, GrB_Index nrows, GrB_Index ncols);

// Creates an empty lazy CFL_Matrix of the given dimensions
//
// The underlying GrB_Matrix is allocated internally
// On success, *matrix is set to a newly allocated CFL_Matrix with new CFL_Matrix
GrB_Info CFL_matrix_create_lazy(CFL_Matrix **matrix, GrB_Index nrows, GrB_Index ncols);

// Free a CFL_Matrix and sets *matrix to NULL
GrB_Info CFL_matrix_free(CFL_Matrix **matrix);

// Recomputes internal fields (nvals, nrows, ncols, format, block_type, etc...)
// from the current state of the underlying base matrix
GrB_Info CFL_matrix_update(CFL_Matrix *matrix);

GrB_Info CFL_mxm(CFL_Matrix *output, CFL_Matrix *first, CFL_Matrix *second, bool accum,
                 bool swap, int8_t optimizations);
GrB_Info CFL_wise(CFL_Matrix *output, CFL_Matrix *first, CFL_Matrix *second, bool accum,
                  int8_t optimizations);

// Computes C = B \ A, i.e. C(i, j) = true iff B(i, j) = true and A(i, j) = false.
GrB_Info CFL_rsub(CFL_Matrix *output, CFL_Matrix *mask, int8_t optimizations);
GrB_Info CFL_dup(CFL_Matrix *output, CFL_Matrix *input, int8_t optimizations);

// Converts a matrix (lazy or regular) into its evaluated base form by merging
// all underlying base matrices. If the input matrix is not lazy, returns its copy
//
// Parameters:
//   matrix_p       - [out] Pointer that will be holds created matrix
//   input          - [in]  Pointer to the source matrix (may be lazy).
//   optimizations  - [in]  Bitmask specifying enabled optimizations.
GrB_Info CFL_matrix_to_base(
    // output
    CFL_Matrix **matrix_p,
    // input
    CFL_Matrix *matrix, int8_t optimizations);

GrB_Info CFL_clear(CFL_Matrix *A, int8_t optimizations);