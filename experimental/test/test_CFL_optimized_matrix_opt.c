//------------------------------------------------------------------------------
// LAGraph/experimental/test/LAGraph_CFL_reachability.c: test cases for
// operations for Optimized Context-Free Language Reachability Matrix-Based
// Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

#define OPT_N_FLAGS 16

char msg[LAGRAPH_MSG_LEN];

typedef CFL_Matrix Matrix;
typedef enum CFL_Matrix_block Matrix_block;

extern GrB_Info matrix_to_format(Matrix *matrix, int32_t format, bool is_bool);
extern GrB_Info matrix_clear_format(Matrix *A, int8_t optimizations);
extern GrB_Info matrix_dup_format(Matrix *output, Matrix *input, int8_t optimizations);
extern GrB_Info block_matrix_hyper_rotate_i(Matrix *matrix, enum CFL_Matrix_block format);
extern void block_matrix_to_diag(Matrix *diag, Matrix *input);
extern void block_matrix_reduce(Matrix *matrix, Matrix *input, int8_t optimizations);

// extern GrB_Info matrix_clear_empty(Matrix *A, int8_t optimizations);
// extern GrB_Info matrix_mxm_empty(Matrix *output, Matrix *first, Matrix *second,
//                                  bool accum, bool swap, int8_t optimizations);
// extern GrB_Info matrix_wise_empty(Matrix *output, Matrix *first, Matrix *second,
//                                   bool accum, int8_t optimizations);

//------------------------------------------------------------------------------
// helpers
//------------------------------------------------------------------------------

static void setup(void) { OK(LAGraph_Init(msg)); }

static void teardown(void) { OK(LAGraph_Finalize(msg)); }

// 0  1  .. 0
// 1  0  .. 0
// .. .. .. ..
// 0  0  .. 0
static Matrix make_simple_matrix(int n) {
    GrB_Matrix A;
    OK(GrB_Matrix_new(&A, GrB_BOOL, n, n));
    OK(GrB_Matrix_setElement_BOOL(A, true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(A, true, 1, 0));
    Matrix M = CFL_matrix_from_base(A);
    return M;
}

// 1  0  .. 0
// 0  1  .. 0
// .. .. .. ..
// 0  0  .. 0
static Matrix make_simple_matrix_inverted(int n) {
    GrB_Matrix A;
    OK(GrB_Matrix_new(&A, GrB_BOOL, n, n));
    OK(GrB_Matrix_setElement_BOOL(A, true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(A, true, 1, 1));
    Matrix M = CFL_matrix_from_base(A);
    return M;
}

// 1  1  .. 1
// 1  1  .. 1
// .. .. .. ..
// 1  1  .. 1
static Matrix make_ones_matrix(int n) {
    GrB_Matrix A;
    OK(GrB_Matrix_new(&A, GrB_BOOL, n, n));

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            OK(GrB_Matrix_setElement_BOOL(A, true, i, j));
        }
    }

    Matrix M = CFL_matrix_from_base(A);
    return M;
}

static void free_matrix(Matrix *M) { CFL_matrix_free(M); }

static bool compare_matrices(Matrix first, Matrix second) {
    Matrix A, B, C;
    A = CFL_matrix_to_base(&first, OPT_LAZY);
    B = CFL_matrix_to_base(&second, OPT_LAZY);
    C = CFL_matrix_create(A.nrows, A.ncols);

    if (A.nvals != B.nvals) {
        CFL_matrix_free(&A);
        CFL_matrix_free(&B);
        CFL_matrix_free(&C);

        return false;
    }

    GrB_eWiseMult(C.base, GrB_NULL, false, GrB_EQ_BOOL, A.base, B.base, GrB_NULL);
    CFL_matrix_update(&C);

    bool result = C.nvals == A.nvals;

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);
    CFL_matrix_free(&C);

    return result;
}

static void test_compare_matrices_function(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = make_simple_matrix(5);
    Matrix B = make_simple_matrix_inverted(5);

    TEST_CHECK(!compare_matrices(A, B));
    TEST_CHECK(compare_matrices(A, A));
    TEST_CHECK(compare_matrices(B, B));

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Сreation and Free
//------------------------------------------------------------------------------

static void test_CFL_create_free(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    CFL_Matrix A = CFL_matrix_create(4, 4);
    TEST_CHECK(A.nrows == 4);
    TEST_CHECK(A.ncols == 4);
    TEST_CHECK(A.base != NULL);

    CFL_matrix_free(&A);
    TEST_CHECK(A.base == NULL);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Multiplication
//------------------------------------------------------------------------------

static void test_CFL_mxm(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (int mask = 0; mask < OPT_N_FLAGS; mask++) {
        Matrix A = make_simple_matrix(2);
        Matrix B = make_simple_matrix(2);
        matrix_to_format(&A, GrB_ROWMAJOR, true);
        Matrix C = CFL_matrix_create(2, 2);

        GrB_Info info = CFL_mxm(&C, &A, &B, false, false, mask);
        TEST_CHECK(info == GrB_SUCCESS);
        TEST_MSG("matrix_mxm_opt failed for mask=%d", mask);

        // Validate result (A*B should have two true element)
        GrB_Index nvals = 0;
        OK(GrB_Matrix_nvals(&nvals, C.base));
        TEST_CHECK(nvals == 2);
        TEST_CHECK(C.nvals == 2);
        TEST_MSG("Unexpected empty result, mask=%d", mask);

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&C);
    }

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Wise
//------------------------------------------------------------------------------

static void test_CFL_wise(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (int mask = 0; mask < OPT_N_FLAGS; mask++) {
        Matrix A = make_simple_matrix(2);
        Matrix B = make_simple_matrix(2);
        matrix_to_format(&A, GrB_ROWMAJOR, true);

        GrB_Info info = CFL_wise(&A, &A, &B, false, mask);
        TEST_CHECK(info == GrB_SUCCESS);
        TEST_MSG("matrix_wise failed for mask=%d", mask);

        // Validate result (A+B should have two true element)
        GrB_Index nvals = 0;
        OK(GrB_Matrix_nvals(&nvals, A.base));
        TEST_CHECK(nvals == 2);
        TEST_CHECK(A.nvals == 2);
        TEST_MSG("Unexpected empty result, mask=%d", mask);

        free_matrix(&A);
        free_matrix(&B);
    }

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Rsub
//------------------------------------------------------------------------------

static void test_CFL_rsub(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (int mask = 0; mask < OPT_N_FLAGS; mask++) {
        Matrix A = make_simple_matrix(2);
        Matrix B = make_simple_matrix(2);
        matrix_to_format(&A, GrB_ROWMAJOR, true);

        GrB_Info info = CFL_rsub(&A, &B, mask);
        TEST_CHECK(info == GrB_SUCCESS);
        TEST_MSG("matrix_rsub failed for mask=%d", mask);

        // Validate result (A-B should have zero true element)
        GrB_Index nvals = 0;
        OK(GrB_Matrix_nvals(&nvals, A.base));
        TEST_CHECK(nvals == 0);
        TEST_CHECK(A.nvals == 0);
        TEST_MSG("Unexpected non-empty result, mask=%d", mask);

        free_matrix(&A);
        free_matrix(&B);
    }

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Clear
//------------------------------------------------------------------------------

static void test_CFL_clear(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (int mask = 0; mask < OPT_N_FLAGS; mask++) {
        Matrix A = make_simple_matrix(3);

        // Clear nonempty matrix
        CFL_clear(&A, mask);
        GrB_Index nvals;
        TEST_CHECK(A.nvals == 0);

        free_matrix(&A);
    }

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Duplicate
//------------------------------------------------------------------------------

static void test_CFL_dup(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (int mask = 0; mask < OPT_N_FLAGS; mask++) {
        CFL_Matrix A = make_simple_matrix(3);
        CFL_Matrix B = CFL_matrix_create(3, 3);

        OK(CFL_dup(&A, &B, mask));

        GrB_Index nvals_A, nvals_B;
        OK(GrB_Matrix_nvals(&nvals_A, A.base));
        OK(GrB_Matrix_nvals(&nvals_B, B.base));

        TEST_CHECK(nvals_A == nvals_B);
        TEST_CHECK(A.nvals == B.nvals);
        free_matrix(&A);
        free_matrix(&B);
    }

    teardown();
#endif
}

static void test_CFL_dup_same(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    CFL_Matrix A = make_simple_matrix(3);
    OK(CFL_dup(&A, &A, 0));

    free_matrix(&A);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Format optimization
//------------------------------------------------------------------------------

static void test_CFL_format_create_rowmajor_matrix(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = make_simple_matrix(5);
    TEST_CHECK(A.is_both == false);
    TEST_CHECK(A.base_col == NULL);
    TEST_CHECK(A.base_row != NULL);
    TEST_CHECK(A.base == A.base_row);
    TEST_CHECK(A.format == GrB_ROWMAJOR);

    CFL_matrix_free(&A);

    teardown();
#endif
}

static void test_CFL_format_to_format(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = make_simple_matrix(5);

    // same
    OK(matrix_to_format(&A, GrB_ROWMAJOR, false));
    TEST_CHECK(A.is_both == false);
    TEST_CHECK(A.base_col == NULL);
    TEST_CHECK(A.base_row != NULL);
    TEST_CHECK(A.base == A.base_row);
    TEST_CHECK(A.format == GrB_ROWMAJOR);

    // switch to another
    OK(matrix_to_format(&A, GrB_COLMAJOR, false));
    TEST_CHECK(A.is_both == false);
    TEST_CHECK(A.base_col != NULL);
    TEST_CHECK(A.base_row == NULL);
    TEST_CHECK(A.base == A.base_col);
    TEST_CHECK(A.format == GrB_COLMAJOR);

    // create both matrices, same format
    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    TEST_CHECK(A.is_both == false);
    TEST_CHECK(A.base_col != NULL);
    TEST_CHECK(A.base_row == NULL);
    TEST_CHECK(A.base == A.base_col);
    TEST_CHECK(A.format == GrB_COLMAJOR);

    // switch to another, both
    OK(matrix_to_format(&A, GrB_ROWMAJOR, true));
    TEST_CHECK(A.is_both == true);
    TEST_CHECK(A.base_col != NULL);
    TEST_CHECK(A.base_row != NULL);
    TEST_CHECK(A.base == A.base_row);
    TEST_CHECK(A.format == GrB_ROWMAJOR);

    // switch to another, when matrix already in both state
    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    TEST_CHECK(A.is_both == true);
    TEST_CHECK(A.base_col != NULL);
    TEST_CHECK(A.base_row != NULL);
    TEST_CHECK(A.base == A.base_col);
    TEST_CHECK(A.format == GrB_COLMAJOR);

    CFL_matrix_free(&A);

    teardown();
#endif
}

static void test_CFL_format_clear_format(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = make_simple_matrix(5);

    // create matrix in both state
    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    OK(matrix_clear_format(&A, OPT_FORMAT));

    TEST_CHECK(A.nvals == 0);

    // NULL Matrix
    GrB_Matrix old_base = A.base;
    GrB_Matrix old_base_row = A.base;
    A.base = NULL;
    A.base_row = NULL;
    GrB_Info result = matrix_clear_format(&A, OPT_FORMAT);
    OK(!result);

    A.base = old_base;
    A.base_row = old_base_row;
    CFL_matrix_free(&A);

    teardown();
#endif
}

static void test_CFL_format_dup_format(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = make_simple_matrix(5);
    Matrix B = CFL_matrix_create(5, 5);

    // create matrix in both state
    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    OK(matrix_dup_format(&A, &B, OPT_FORMAT));

    TEST_CHECK(A.nvals == B.nvals);
    TEST_CHECK(A.base != B.base);

    GrB_Index nvals_A = A.nvals;
    GrB_Index nvals_B = B.nvals;

    GrB_Index *rows_A = malloc(nvals_A * sizeof(GrB_Index));
    GrB_Index *cols_A = malloc(nvals_A * sizeof(GrB_Index));
    bool *vals_A = malloc(nvals_A * sizeof(bool));

    GrB_Index *rows_B = malloc(nvals_B * sizeof(GrB_Index));
    GrB_Index *cols_B = malloc(nvals_B * sizeof(GrB_Index));
    bool *vals_B = malloc(nvals_B * sizeof(bool));

    GrB_Matrix_extractTuples_BOOL(rows_A, cols_A, vals_A, &nvals_A, A.base);
    GrB_Matrix_extractTuples_BOOL(rows_B, cols_B, vals_B, &nvals_B, B.base);

    bool equal = true;
    for (GrB_Index i = 0; i < nvals_A; i++) {
        bool found = false;
        for (GrB_Index j = 0; j < nvals_B; j++) {
            if (rows_A[i] == rows_B[j] && cols_A[i] == cols_B[j] &&
                vals_A[i] == vals_B[j]) {
                found = true;
                break;
            }
        }
        if (!found) {
            equal = false;
            break;
        }
    }
    TEST_CHECK(equal);

    // NULL Matrix
    GrB_Matrix old_base = A.base;
    GrB_Matrix old_base_row = A.base;
    A.base = NULL;
    A.base_row = NULL;
    GrB_Info result = matrix_dup_format(&A, &B, OPT_FORMAT);
    OK(!result);

    A.base = old_base;
    A.base_row = old_base_row;
    CFL_matrix_free(&A);
    CFL_matrix_free(&B);

    // Without optimization flag
    A = make_simple_matrix(5);
    B = CFL_matrix_create(5, 5);
    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    OK(matrix_dup_format(&A, &B, 0));

    teardown();
#endif
}

static void test_CFL_format_mxm_second_greather_then_k(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);
    Matrix B = make_ones_matrix(5);
    Matrix C = CFL_matrix_create(5, 5);

    OK(matrix_to_format(&A, GrB_COLMAJOR, false));
    OK(CFL_mxm(&C, &A, &B, false, false, OPT_FORMAT));
    TEST_CHECK(C.nvals == 0);

    teardown();
#endif
}

static void test_CFL_format_wise_when_both(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);
    Matrix B = make_ones_matrix(5);

    OK(matrix_to_format(&A, GrB_COLMAJOR, true));
    OK(CFL_wise(&A, &A, &B, false, OPT_FORMAT));
    TEST_CHECK(A.nvals == 25);

    // NULL Matrix
    GrB_Matrix old_base = A.base;
    GrB_Matrix old_base_row = A.base;
    A.base = NULL;
    A.base_row = NULL;
    GrB_Info result = CFL_wise(&A, &A, &B, false, OPT_FORMAT);
    A.base = old_base;
    A.base_row = old_base_row;
    OK(!result);

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Empty optimization
//------------------------------------------------------------------------------

static void test_CFL_empty_clear_empty_matrix(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);

    OK(CFL_clear(&A, OPT_EMPTY));
    TEST_CHECK(A.nvals == 0);

    CFL_matrix_free(&A);

    teardown();
#endif
}

static void test_CFL_empty_dup(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);
    Matrix B = make_ones_matrix(5);

    OK(CFL_dup(&A, &B, OPT_EMPTY));
    TEST_CHECK(A.nvals == 25);

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);

    teardown();
#endif
}

static void test_CFL_empty_mxm_one_is_empty(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A, B, C;

    // accum
    A = CFL_matrix_create(5, 5);
    B = make_ones_matrix(5);
    C = make_ones_matrix(5);

    OK(CFL_mxm(&C, &A, &B, true, false, OPT_EMPTY));
    TEST_CHECK(C.nvals == 25);

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);
    CFL_matrix_free(&C);

    // without accum
    A = CFL_matrix_create(5, 5);
    B = make_ones_matrix(5);
    C = make_ones_matrix(5);

    OK(CFL_mxm(&C, &A, &B, false, false, OPT_EMPTY));
    TEST_CHECK(C.nvals == 0);

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);
    CFL_matrix_free(&C);

    // output empty
    A = CFL_matrix_create(5, 5);
    B = make_ones_matrix(5);
    C = CFL_matrix_create(5, 5);

    OK(CFL_mxm(&C, &A, &B, false, false, OPT_EMPTY));
    TEST_CHECK(C.nvals == 0);

    CFL_matrix_free(&A);
    CFL_matrix_free(&B);
    CFL_matrix_free(&C);

    teardown();
#endif
}

// wise(C, A, B, accum)
// matrix may be empty and full
// we have 16 states
static void test_CFL_empty_wise(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A, B, C;

    for (size_t is_first_empty = 0; is_first_empty < 2; is_first_empty++) {
        for (size_t is_second_empty = 0; is_second_empty < 2; is_second_empty++) {
            for (size_t is_output_empty = 0; is_output_empty < 2; is_output_empty++) {
                for (size_t is_accum = 0; is_accum < 2; is_accum++) {
                    A = is_first_empty ? CFL_matrix_create(5, 5) : make_ones_matrix(5);
                    B = is_second_empty ? CFL_matrix_create(5, 5) : make_ones_matrix(5);
                    C = is_output_empty ? CFL_matrix_create(5, 5) : make_ones_matrix(5);

                    CFL_wise(&C, &A, &B, is_accum, OPT_EMPTY);

                    TEST_CHECK(A.nvals == (is_first_empty ? 0 : 25));
                    TEST_CHECK(B.nvals == (is_second_empty ? 0 : 25));

                    if (is_first_empty && is_second_empty) {
                        if (!is_accum) {
                            TEST_CHECK(C.nvals == 0);
                        } else {
                            TEST_CHECK(C.nvals == (is_output_empty ? 0 : 25));
                        }
                    } else {
                        TEST_CHECK(C.nvals == 25);
                    }

                    CFL_matrix_free(&A);
                    CFL_matrix_free(&B);
                    CFL_matrix_free(&C);
                }
            }
        }
    }

    // if output equal first argument (iadd)
    for (size_t is_first_empty = 0; is_first_empty < 2; is_first_empty++) {
        for (size_t is_second_empty = 0; is_second_empty < 2; is_second_empty++) {
            for (size_t is_accum = 0; is_accum < 2; is_accum++) {
                A = is_first_empty ? CFL_matrix_create(5, 5) : make_ones_matrix(5);
                B = is_second_empty ? CFL_matrix_create(5, 5) : make_ones_matrix(5);

                CFL_wise(&A, &A, &B, is_accum, OPT_EMPTY);

                TEST_CHECK(B.nvals == (is_second_empty ? 0 : 25));

                if (is_second_empty) {
                    TEST_CHECK(A.nvals == (is_first_empty ? 0 : 25));
                } else {
                    TEST_CHECK(A.nvals == 25);
                }

                CFL_matrix_free(&A);
                CFL_matrix_free(&B);
            }
        }
    }

    teardown();
#endif
}

static void test_CFL_empty_rsub_both_empty(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A, B;
    A = CFL_matrix_create(5, 5);
    B = CFL_matrix_create(5, 5);

    OK(CFL_rsub(&A, &B, OPT_EMPTY));
    TEST_CHECK(A.nvals == 0);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Lazy optimization
//------------------------------------------------------------------------------

static void test_CFL_lazy_create(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    GrB_Matrix A_base;
    GrB_Matrix_new(&A_base, GrB_BOOL, 5, 5);
    Matrix A = CFL_matrix_from_base_lazy(A_base);
    TEST_CHECK(A.base_matrices != NULL);
    TEST_CHECK(A.base_matrices_count == 1);
    TEST_CHECK(A.is_lazy == true);
    TEST_CHECK(A.nvals == 0);

    teardown();
#endif
}

// Create lazy matrix
// 1 1 1 0 .. 0 (count of 1: 111)
// 1 0 0 0 .. 0 (count of 1: 1)
// 1 1 0 0 .. 0 (count of 1: 11)
// . . . . .. .
// 0 0 0 0 .. 0
static Matrix make_lazy_matrix(size_t base_matrices_count) {
    size_t n = 0;
    Matrix base_matrices[base_matrices_count];
    // 1 -> 1, 2 -> 11, 3 -> 111, 4 -> 1111
    for (size_t i = 0; i < base_matrices_count; i++) {
        n *= 10;
        n++;
    }

    size_t tmp_n = 0;
    for (size_t i = 0; i < base_matrices_count; i++) {
        tmp_n *= 10;
        tmp_n++;
        GrB_Matrix base_matrix;
        GrB_Matrix_new(&base_matrix, GrB_BOOL, n, n);
        for (size_t j = 0; j < tmp_n; j++) {
            // i+1 for future testing of sorting matrices
            GrB_Matrix_setElement_BOOL(base_matrix, true, (i + 1) % base_matrices_count,
                                       j);
        }

        base_matrices[(i + 1) % base_matrices_count] = CFL_matrix_from_base(base_matrix);
    }

    GrB_Matrix _result;
    GrB_Matrix_new(&_result, GrB_BOOL, n, n);
    Matrix result = CFL_matrix_from_base_lazy(_result);

    result.is_lazy = true;
    result.base_matrices_count = base_matrices_count;
    for (size_t i = 0; i < base_matrices_count; i++) {
        result.base_matrices[i] = base_matrices[i];
    }

    CFL_matrix_update(&result);

    return result;
}

static void test_CFL_lazy_mxm(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        for (size_t is_reverse = 0; is_reverse < 2; is_reverse++) {
            Matrix A = make_lazy_matrix(4);
            Matrix B = make_ones_matrix(1111);
            Matrix C = CFL_matrix_create(1111, 1111);

            Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
            Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);
            Matrix C_base = CFL_matrix_to_base(&C, OPT_LAZY);

            CFL_mxm(&C, &A, &B, is_accum, is_reverse, OPT_LAZY);
            CFL_mxm(&C_base, &A_base, &B_base, is_accum, is_reverse, 0);

            TEST_CHECK(compare_matrices(C, C_base));

            free_matrix(&A);
            free_matrix(&B);
            free_matrix(&C);
            free_matrix(&A_base);
            free_matrix(&B_base);
            free_matrix(&C_base);
        }
    }

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        for (size_t is_reverse = 0; is_reverse < 2; is_reverse++) {
            Matrix A = make_lazy_matrix(4);
            Matrix B = CFL_matrix_create(1111, 1111);
            Matrix C = CFL_matrix_create(1111, 1111);

            Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
            Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);
            Matrix C_base = CFL_matrix_to_base(&C, OPT_LAZY);

            CFL_mxm(&C, &A, &B, is_accum, is_reverse, OPT_LAZY);
            CFL_mxm(&C_base, &A_base, &B_base, is_accum, is_reverse, 0);

            TEST_CHECK(compare_matrices(C, C_base));

            free_matrix(&A);
            free_matrix(&B);
            free_matrix(&C);
            free_matrix(&A_base);
            free_matrix(&B_base);
            free_matrix(&C_base);
        }
    }

    teardown();
#endif
}

static void test_CFL_lazy_wise(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        Matrix A = make_lazy_matrix(4);
        Matrix B = CFL_matrix_create(1111, 1111);
        Matrix C = CFL_matrix_create(1111, 1111);

        Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
        Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);
        Matrix C_base = CFL_matrix_to_base(&C, OPT_LAZY);

        CFL_wise(&C, &A, &B, is_accum, OPT_LAZY);
        CFL_wise(&C_base, &A_base, &B_base, is_accum, 0);

        TEST_CHECK(compare_matrices(A, C_base));

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&C);
        free_matrix(&A_base);
        free_matrix(&B_base);
        free_matrix(&C_base);
    }

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        Matrix A = CFL_matrix_create(1111, 1111);
        Matrix B = make_lazy_matrix(4);
        Matrix C = CFL_matrix_create(1111, 1111);

        Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
        Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);
        Matrix C_base = CFL_matrix_to_base(&C, OPT_LAZY);

        CFL_wise(&C, &A, &B, is_accum, OPT_LAZY);
        CFL_wise(&C_base, &A_base, &B_base, is_accum, 0);

        TEST_CHECK(compare_matrices(C, C_base));

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&C);
        free_matrix(&A_base);
        free_matrix(&B_base);
        free_matrix(&C_base);
    }

    teardown();
#endif
}

static void test_CFL_lazy_rsub(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    {
        Matrix A = make_lazy_matrix(4);
        Matrix B = CFL_matrix_create(1111, 1111);

        Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
        Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);

        CFL_rsub(&A, &B, OPT_LAZY);
        CFL_rsub(&A_base, &B_base, 0);

        TEST_CHECK(compare_matrices(A, A_base));

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&A_base);
        free_matrix(&B_base);
    }

    {
        Matrix A = CFL_matrix_create(1111, 1111);
        Matrix B = make_lazy_matrix(4);

        Matrix A_base = CFL_matrix_to_base(&A, OPT_LAZY);
        Matrix B_base = CFL_matrix_to_base(&B, OPT_LAZY);

        CFL_rsub(&A, &B, OPT_LAZY);
        CFL_rsub(&A_base, &B_base, 0);

        TEST_CHECK(compare_matrices(A, A_base));

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&A_base);
        free_matrix(&B_base);
    }

    teardown();
#endif
}

//------------------------------------------------------------------------------
// Block optimization
//------------------------------------------------------------------------------

// Create block vector
static Matrix make_block_vector(size_t n, size_t alpha, Matrix_block format) {
    GrB_Matrix result;
    switch (format) {
    case CELL:
        // 1 1 0 0 0 0 .
        // 0 1 1 0 0 0 .
        // 0 0 1 1 0 0 .
        // 0 0 0 1 1 0 .
        // 0 0 0 0 1 1 .
        // 1 0 0 0 0 1 .
        // 1 1 0 0 0 0 .
        // . . . . . . .
        GrB_Matrix_new(&result, GrB_BOOL, n, n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                bool value = j == (i % n) || j == ((i + 1) % n) ? true : false;
                if (!value)
                    continue;

                GrB_Matrix_setElement_BOOL(result, value, i, j);
            }
        }

        return CFL_matrix_from_base(result);

    case VEC_VERT:
        GrB_Matrix_new(&result, GrB_BOOL, n * alpha, n);
        for (size_t i = 0; i < n * alpha; i++) {
            for (size_t j = 0; j < n; j++) {
                bool value = j == (i + 2 % n) || j == ((i + 3) % n) ? true : false;
                if (!value)
                    continue;

                GrB_Matrix_setElement_BOOL(result, value, i, j);
            }
        }

        return CFL_matrix_from_base(result);

    case VEC_HORIZ:
        GrB_Matrix_new(&result, GrB_BOOL, n, n * alpha);
        for (size_t i = 0; i < n * alpha; i++) {
            for (size_t j = 0; j < n; j++) {
                bool value = j == (i + 5 % n) || j == ((i + 6) % n) ? true : false;
                if (!value)
                    continue;

                GrB_Matrix_setElement_BOOL(result, value, i, j);
            }
        }

        return CFL_matrix_from_base(result);
    }
}

static Matrix make_lazy_block_matrix(size_t base_matrices_count) {
    size_t n = 0;
    Matrix base_matrices[base_matrices_count];
    // 1 -> 1, 2 -> 11, 3 -> 111, 4 -> 1111
    for (size_t i = 0; i < base_matrices_count; i++) {
        n *= 10;
        n++;
    }

    size_t tmp_n = 0;
    for (size_t i = 0; i < base_matrices_count; i++) {
        tmp_n *= 10;
        tmp_n++;
        GrB_Matrix base_matrix;
        GrB_Matrix_new(&base_matrix, GrB_BOOL, n, n);
        for (size_t j = 0; j < tmp_n; j++) {
            // i+1 for future testing of sorting matrices
            GrB_Matrix_setElement_BOOL(base_matrix, true, (i + 1) % base_matrices_count,
                                       j);
        }

        base_matrices[(i + 1) % base_matrices_count] = CFL_matrix_from_base(base_matrix);
    }

    GrB_Matrix _result;
    GrB_Matrix_new(&_result, GrB_BOOL, n, n);
    Matrix result = CFL_matrix_from_base_lazy(_result);

    result.is_lazy = true;
    result.base_matrices_count = base_matrices_count;
    for (size_t i = 0; i < base_matrices_count; i++) {
        result.base_matrices[i] = base_matrices[i];
    }

    CFL_matrix_update(&result);

    return result;
}

// TODO: create equality check with basic mxm
static void test_CFL_block_mxm(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        for (size_t is_reverse = 0; is_reverse < 2; is_reverse++) {
            for (size_t first_format = 0; first_format < 3; first_format++) {
                for (size_t second_format = 0; second_format < 3; second_format++) {
                    Matrix A = make_block_vector(20, 20, first_format);
                    Matrix B = make_block_vector(20, 20, second_format);

                    Matrix C;
                    if (first_format != CELL && second_format != CELL) {
                        size_t nrows = first_format == VEC_HORIZ ? 20 : 20 * 20;
                        size_t ncols = first_format == VEC_HORIZ ? 20 * 20 : 20;
                        C = CFL_matrix_create(nrows, ncols);
                    } else {
                        C = CFL_matrix_create(first_format == CELL ? 20 : 20 * 20,
                                              second_format == CELL ? 20 : 20 * 20);
                    }

                    OK(CFL_mxm(&C, &A, &B, is_accum, is_reverse, OPT_BLOCK));

                    free_matrix(&A);
                    free_matrix(&B);
                    free_matrix(&C);
                }
            }
        }
    }

    teardown();
#endif
}

static void test_CFL_block_dup(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    {
        Matrix A = make_block_vector(20, 20, VEC_HORIZ);
        Matrix B = CFL_matrix_create(20 * 20, 20);

        OK(CFL_dup(&B, &A, OPT_BLOCK));
        TEST_CHECK(B.nvals == A.nvals);
    }

    {
        Matrix A = make_block_vector(20, 20, VEC_VERT);
        Matrix B = CFL_matrix_create(20 * 20, 20);

        OK(CFL_dup(&B, &A, OPT_BLOCK));
        TEST_CHECK(B.nvals == A.nvals);
    }

    {
        Matrix A = make_block_vector(20, 20, VEC_HORIZ);
        Matrix B = CFL_matrix_create(20, 20 * 20);

        OK(CFL_dup(&B, &A, OPT_BLOCK));
        TEST_CHECK(B.nvals == A.nvals);
    }

    {
        Matrix A = make_block_vector(20, 20, VEC_VERT);
        Matrix B = CFL_matrix_create(20, 20 * 20);

        OK(CFL_dup(&B, &A, OPT_BLOCK));
        TEST_CHECK(B.nvals == A.nvals);
    }

    teardown();
#endif
}

static void test_CFL_block_wise(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    {
        Matrix A = make_block_vector(20, 20, VEC_HORIZ);
        Matrix B = make_block_vector(20, 20, VEC_HORIZ);
        Matrix C = make_block_vector(20, 20, VEC_HORIZ);

        TEST_CHECK(CFL_wise(&C, &A, &B, false, OPT_BLOCK) == GrB_INVALID_VALUE);

        free_matrix(&A);
        free_matrix(&B);
        free_matrix(&C);
    }

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        for (size_t is_reverse = 0; is_reverse < 2; is_reverse++) {
            for (size_t first_format = 0; first_format < 3; first_format++) {
                for (size_t second_format = 0; second_format < 3; second_format++) {
                    Matrix A = make_block_vector(20, 20, first_format);
                    Matrix B = make_block_vector(20, 20, second_format);

                    OK(CFL_wise(&A, &A, &B, is_accum, OPT_BLOCK));

                    free_matrix(&A);
                    free_matrix(&B);
                }
            }
        }
    }

    teardown();
#endif
}

static void test_CFL_block_rsub(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    for (size_t is_accum = 0; is_accum < 2; is_accum++) {
        for (size_t is_reverse = 0; is_reverse < 2; is_reverse++) {
            for (size_t first_format = 0; first_format < 3; first_format++) {
                for (size_t second_format = 0; second_format < 3; second_format++) {
                    Matrix A = make_block_vector(20, 20, first_format);
                    Matrix B = make_block_vector(20, 20, second_format);

                    if ((first_format != CELL && second_format == CELL) ||
                        (first_format == CELL && second_format != CELL)) {
                        TEST_CHECK(CFL_rsub(&A, &B, OPT_BLOCK) == GrB_INVALID_VALUE);
                    } else {
                        OK(CFL_rsub(&A, &B, OPT_BLOCK));
                    }

                    free_matrix(&A);
                    free_matrix(&B);
                }
            }
        }
    }

    teardown();
#endif
}

static void test_CFL_block_hyper_rotate(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    // lazy
    {
        GrB_Matrix _A;
        GrB_Matrix_new(&_A, GrB_BOOL, 20 * 20, 20);
        Matrix A = CFL_matrix_from_base_lazy(_A);

        GrB_Matrix _A0;
        GrB_Matrix_new(&_A0, GrB_BOOL, 20 * 20, 20);
        Matrix A0 = CFL_matrix_from_base(_A0);
        GrB_Matrix _A1;
        GrB_Matrix_new(&_A1, GrB_BOOL, 20 * 20, 20);
        Matrix A1 = CFL_matrix_from_base(_A0);
        GrB_Matrix _A2;
        GrB_Matrix_new(&_A2, GrB_BOOL, 20 * 20, 20);
        Matrix A2 = CFL_matrix_from_base(_A0);

        A.base_matrices[0] = A0;
        A.base_matrices[1] = A1;
        A.base_matrices[2] = A2;
        A.base_matrices_count = 3;

        OK(block_matrix_hyper_rotate_i(&A, VEC_HORIZ));
        TEST_CHECK(A.nrows == 20);
        TEST_CHECK(A.ncols == 20 * 20);

        free_matrix(&A);
        free_matrix(&A0);
        free_matrix(&A1);
        free_matrix(&A2);
    }

    // cell
    {
        Matrix A = CFL_matrix_create(5, 5);

        OK(block_matrix_hyper_rotate_i(&A, VEC_HORIZ));

        free_matrix(&A);
    }

    // VERT
    {
        Matrix A = make_block_vector(5, 5, VEC_VERT);
        TEST_CHECK(A.nvals != 0);
        TEST_CHECK(A.is_lazy != true);

        OK(block_matrix_hyper_rotate_i(&A, VEC_HORIZ));

        free_matrix(&A);
    }

    teardown();
#endif
}

static void test_CFL_block_to_diag(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);
    Matrix B = CFL_matrix_create(5, 5);

    GrB_Index nvals_before;
    GrB_Matrix_nvals(&nvals_before, A.base);
    block_matrix_to_diag(&A, &B);
    GrB_Index nvals_after;
    GrB_Matrix_nvals(&nvals_after, A.base);

    // A didn't change
    TEST_CHECK(nvals_after == nvals_before);

    free_matrix(&A);
    free_matrix(&B);

    teardown();
#endif
}

static void test_CFL_block_to_reduce(void) {
#if LAGRAPH_SUITESPARSE
    setup();

    Matrix A = CFL_matrix_create(5, 5);
    Matrix B = make_ones_matrix(5);

    // cell matrix cause just dup
    block_matrix_reduce(&A, &B, OPT_BLOCK);
    TEST_CHECK(A.nvals == B.nvals);

    free_matrix(&A);
    free_matrix(&B);

    teardown();
#endif
}

//------------------------------------------------------------------------------
// TEST LIST
//------------------------------------------------------------------------------

TEST_LIST = {
    {"test_compare_matrices_function", test_compare_matrices_function},
    {"test CFL create and free", test_CFL_create_free},
    {"test CFL mxm", test_CFL_mxm},
    {"test CFL clear", test_CFL_clear},
    {"test CFL dup", test_CFL_dup},
    {"test_CFL_dup_same", test_CFL_dup_same},
    {"test_CFL_wise", test_CFL_wise},
    {"test_CFL_rsub", test_CFL_rsub},
    {"test_CFL_format_create_rowmajor_matrix", test_CFL_format_create_rowmajor_matrix},
    {"test_CFL_format_to_format", test_CFL_format_to_format},
    {"test_CFL_format_clear_format", test_CFL_format_clear_format},
    {"test_CFL_format_dup_format", test_CFL_format_dup_format},
    {"test_CFL_format_wise_when_both", test_CFL_format_wise_when_both},
    {"test_CFL_format_mxm_second_greather_then_k",
     test_CFL_format_mxm_second_greather_then_k},
    {"test_CFL_empty_clear_empty_matrix", test_CFL_empty_clear_empty_matrix},
    {"test_CFL_empty_dup", test_CFL_empty_dup},
    {"test_CFL_empty_mxm_one_is_empty", test_CFL_empty_mxm_one_is_empty},
    {"test_CFL_empty_wise", test_CFL_empty_wise},
    {"test_CFL_empty_rsub_both_empty", test_CFL_empty_rsub_both_empty},
    {"test_CFL_lazy_create", test_CFL_lazy_create},
    {"test_CFL_lazy_mxm", test_CFL_lazy_mxm},
    {"test_CFL_lazy_wise", test_CFL_lazy_wise},
    {"test_CFL_lazy_rsub", test_CFL_lazy_rsub},
    {"test_CFL_block_mxm", test_CFL_block_mxm},
    {"test_CFL_block_wise", test_CFL_block_wise},
    {"test_CFL_block_rsub", test_CFL_block_rsub},
    {"test_CFL_block_dup", test_CFL_block_dup},
    {"test_CFL_block_hyper_rotate", test_CFL_block_hyper_rotate},
    {"test_CFL_block_to_diag", test_CFL_block_to_diag},
    {"test_CFL_block_to_reduce", test_CFL_block_to_reduce},
    {NULL, NULL}};