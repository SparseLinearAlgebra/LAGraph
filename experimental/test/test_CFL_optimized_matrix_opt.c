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

//------------------------------------------------------------------------------
// helpers
//------------------------------------------------------------------------------

static void setup(void) { OK(LAGraph_Init(msg)); }

static void teardown(void) { OK(LAGraph_Finalize(msg)); }

// 0  1  ..
// 1  0  ..
// .. .. ..
static Matrix make_simple_matrix(int n) {
  GrB_Matrix A;
  OK(GrB_Matrix_new(&A, GrB_BOOL, n, n));
  OK(GrB_Matrix_setElement_BOOL(A, true, 0, 1));
  OK(GrB_Matrix_setElement_BOOL(A, true, 1, 0));
  Matrix M = CFL_matrix_from_base(A);
  return M;
}

static void free_matrix(Matrix *M) { CFL_matrix_free(M); }

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

    OK(CFL_dup(&A, &B, OPT_FORMAT));

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

//------------------------------------------------------------------------------
// TEST LIST
//------------------------------------------------------------------------------

TEST_LIST = {{"test CFL create and free", test_CFL_create_free},
             {"test CFL mxm", test_CFL_mxm},
             {"test CFL clear", test_CFL_clear},
             {"test CFL dup", test_CFL_dup},
             {NULL, NULL}};