//------------------------------------------------------------------------------
// LAGraph_RPQMatrix: regular path query algortithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Rodion Suvorov, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

// Code is based on the algorithm described in the following paper:
//  * Diego Arroyuelo, Adrián Gómez-Brandón & Gonzalo Navarro "Evaluating
//    regular path queries on compressed adjacency matrices"
//  * URL: https://link.springer.com/article/10.1007/s00778-024-00885-6

//------------------------------------------------------------------------------
// LAGraph_RPQMatrix: regular path query algortithm
//
// For an edge-labelled directed graph the algorithm computes the nubmer of
// nonzero elements in its reachability matrix.
// The reachability matrix created by following rules:
// * A[i,j] = True if node with index j is reachable from node with index i
//   and concatenation of labels over path between these two labels is a word
//   from specified regular language.
// * A[i,j] = False in other cases.
//
// The algorithm is based on the idea of ​​considering a regular constraint as
// an abstract syntax tree, the leaves of which are matrices of adjacency matrix
// decomposition of the graph, and the internal nodes are the operations of
// conjunction, concatenation, etc.
//
// Example of adjacency matrix decomposition:
//
// Graph:
// (0) --[a]-> (1)
//  |           ^
// [b]    [c]--/
//  |  --/
//  v /
// (2) --[b]-> (3)
//
// Adjacency matrix decomposition of this graph consists of:
// * Adjacency matrix for the label a:
//       0   1   2   3
//   0 |   | T |   |   |
//   1 |   |   |   |   |
//   2 |   |   |   |   |
//   3 |   |   |   |   |
// * Adjacency matrix for the label b:
//       0   1   2   3
//   0 |   |   | T |   |
//   1 |   |   |   |   |
//   2 |   |   |   | T |
//   3 |   |   |   |   |
// * Adjacency matrix for the label c:
//       0   1   2   3
//   0 |   |   |   |   |
//   1 |   |   |   |   |
//   2 |   | T |   |   |
//   3 |   |   |   |   |
//
// The algorithm recursively starts from the root of the given tree and
// performs the operations corresponding to each node on the children of that
// node. As a result of the algorithm's execution, the reachability
// matrix will be stored at the root.
//
// Example of regular expression and its corresponding AST:
//
// Regular expression:
// a/(b|c)*
//
// Abstract syntax tree:
//    ┌─┐
//    │/| (3)
//    └┬┘
// ┌─┬─┴─┬─┐
// │a│   │*│ (2)
// └─┘   └┬┘
//       ┌┴┐
//       │|│ (1)
//       └┬┘
//    ┌─┬─┴─┬─┐
//    │b│   │c│
//    └─┘   └─┘
// The numbers next to the graph nodes show the order in which operations are
// executed. For the decomposition and AST specified above, the resulting
// matrix will have the following structure (Note, that * represents the
// reflexive-transitive closure):
//
//      0   1   2   3
//  0 |   | T |   |   |
//  1 |   |   |   |   |
//  2 |   |   |   |   |
//  3 |   |   |   |   |
//
// So for this example LAGraph_RPQMatrix will return 1.
//
// Full description available at:
//   https://arxiv.org/pdf/2307.14930

#define LG_FREE_WORK \
    {                \
    }

#define LG_FREE_ALL   \
    {                 \
        LG_FREE_WORK ; \
    }

#include "LG_internal.h"
#include "LAGraphX.h"
#include <time.h>
#include <assert.h>

#define OK(s)                                           \
{                                                       \
    GrB_Info info = (s) ;                               \
    if (info != GrB_SUCCESS)                            \
    {                                                   \
        printf("Message: %s\n", msg) ;                  \
        fprintf(stderr, "GraphBLAS error: %d (%s, %d)\n", info, __FILE__, __LINE__) ; \
        return info ;                                   \
    }                                                   \
}

char msg[LAGRAPH_MSG_LEN] ;

#include <stdbool.h>
#include <stdio.h>

GrB_Info LAGraph_RPQMatrix_check(RPQMatrixPlan *plan, GrB_Index *dimension, char *msg)
{
    if (plan == NULL)
    {
        return GrB_SUCCESS ;
    }
    if (plan->op == RPQ_MATRIX_OP_LABEL)
    {
        LG_ASSERT(plan->mat != NULL, GrB_NULL_POINTER) ;
        GrB_Index nrows, ncols ;
        OK(GrB_Matrix_nrows(&nrows, plan->mat)) ;
        OK(GrB_Matrix_ncols(&ncols, plan->mat)) ;
        if (*dimension == -1)
        {
            LG_ASSERT_MSG(nrows == ncols, GrB_INVALID_VALUE,
                          "all the matrices in the graph adjacency matrix decomposition "
                          "should have the same dimensions and be square") ;
            *dimension = ncols ;
        }
        else
        {
            LG_ASSERT_MSG(nrows == *dimension && ncols == *dimension, GrB_INVALID_VALUE,
                          "all the matrices in the graph adjacency matrix decomposition "
                          "should have the same dimensions and be square") ;
        }

        return GrB_SUCCESS ;
    }
    GrB_Info lstatus = LAGraph_RPQMatrix_check(plan->lhs, dimension, msg) ;
    GrB_Info rstatus = LAGraph_RPQMatrix_check(plan->rhs, dimension, msg) ;
    if (rstatus || lstatus)
    {
        return GrB_INVALID_VALUE ;
    }
    return GrB_SUCCESS ;
}

static GrB_Semiring sr = GrB_NULL ;
static GrB_Monoid op = GrB_NULL ;

GrB_Info LAGraph_RPQMatrix_Free(GrB_Matrix *mat) {
    OK(GrB_Matrix_free(mat));
    return GrB_SUCCESS ;
}

static GrB_Info LAGraph_RPQMatrix_storage_to_orientation(RPQMatrixStorage storage, int32_t *orientation)
{
    switch (storage)
    {
    case RPQ_MATRIX_STORAGE_CSC:
        *orientation = GrB_COLMAJOR ;
        return GrB_SUCCESS ;
    case RPQ_MATRIX_STORAGE_CSR:
        *orientation = GrB_ROWMAJOR ;
        return GrB_SUCCESS ;
    default:
        return GrB_INVALID_VALUE ;
    }
}

GrB_Info LAGraph_RPQMatrix_SetGlobalStorageOrientation(RPQMatrixStorage storage)
{
    int32_t orientation = 0 ;
    OK(LAGraph_RPQMatrix_storage_to_orientation(storage, &orientation)) ;
    OK(GrB_Global_set_INT32(GrB_GLOBAL, orientation, GrB_STORAGE_ORIENTATION_HINT)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_SetStorageOrientation(GrB_Matrix mat, RPQMatrixStorage storage)
{
    LG_ASSERT(mat != NULL, GrB_NULL_POINTER) ;
    int32_t orientation = 0 ;
    OK(LAGraph_RPQMatrix_storage_to_orientation(storage, &orientation)) ;
    OK(GrB_Matrix_set_INT32(mat, orientation, GrB_STORAGE_ORIENTATION_HINT)) ;
    OK(GrB_Matrix_wait(mat, GrB_MATERIALIZE)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_DupWithStorageOrientation( GrB_Matrix *dst, GrB_Matrix src, RPQMatrixStorage storage)
{
    LG_ASSERT(dst != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(src != NULL, GrB_NULL_POINTER) ;
    OK(GrB_Matrix_dup(dst, src)) ;
    OK(LAGraph_RPQMatrix_SetStorageOrientation(*dst, storage)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_label(GrB_Matrix *mat, GrB_Index x, GrB_Index i, GrB_Index j)
{
    OK(GrB_Matrix_new(mat, GrB_BOOL, i, j)) ;
    OK(GrB_Matrix_setElement(*mat, true, x, x)) ;
    return (GrB_SUCCESS) ;
}
GrB_Info LAGraph_DestroyRpqMatrixPlan(RPQMatrixPlan *plan)
{
    if (plan == NULL)
    {
        return GrB_SUCCESS ;
    }
    if (plan->res_mat != NULL && plan->mat != plan->res_mat)
    {
        OK(GrB_Matrix_free(&(plan->res_mat))) ;
    }
    GrB_Info lstatus = LAGraph_DestroyRpqMatrixPlan(plan->lhs) ;
    GrB_Info rstatus = LAGraph_DestroyRpqMatrixPlan(plan->rhs) ;
    if (rstatus || lstatus)
    {
        return GrB_INVALID_VALUE ;
    }
    return GrB_SUCCESS ;
} ;

GrB_Info LAGraph_RPQMatrix_solver(RPQMatrixPlan *plan, char *msg) ;

GrB_Info LAGraph_RPQMatrix_reduce(GrB_Index *res, GrB_Matrix mat, uint8_t reduce_type)
{
    GrB_Index nvals ;
    GrB_Vector reduce = GrB_NULL ;

    GrB_Index nrows ;
    OK(GrB_Matrix_nrows(&nrows, mat)) ;
    OK(GrB_Vector_new(&reduce, GrB_BOOL, nrows)) ;

    if (reduce_type == 0)
    {
        OK(GrB_reduce(reduce, GrB_NULL, GrB_NULL, GxB_ANY_BOOL_MONOID, mat, GrB_NULL)) ;
    }
    else if (reduce_type == 1)
    {
        OK(GrB_reduce(reduce, GrB_NULL, GrB_NULL, GxB_ANY_BOOL_MONOID, mat, GrB_DESC_T0)) ;
    }
    else
    {
        OK(GrB_Vector_free(&reduce)) ;
        return GrB_INVALID_VALUE ;
    }

    OK(GrB_Vector_nvals(&nvals, reduce)) ;
    *res = nvals ;

    OK(GrB_Vector_free(&reduce)) ;
    return (GrB_SUCCESS) ;
}

static GrB_Info LAGraph_RPQMatrixLor(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG(plan->op == RPQ_MATRIX_OP_LOR, GrB_INVALID_VALUE, "operator is not lor") ;
    LG_ASSERT_MSG(plan->res_mat == NULL, GrB_INVALID_VALUE, "resulting matrix is already set as lor result") ;

    RPQMatrixPlan *lhs = plan->lhs ;
    RPQMatrixPlan *rhs = plan->rhs ;

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER) ;

    OK(LAGraph_RPQMatrix_solver(lhs, msg)) ;
    OK(LAGraph_RPQMatrix_solver(rhs, msg)) ;

    LG_ASSERT(rhs->res_mat != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(lhs->res_mat != NULL, GrB_NULL_POINTER) ;

    GrB_Matrix lhs_mat = lhs->res_mat ;
    GrB_Matrix rhs_mat = rhs->res_mat ;

    GrB_Index dimension ;
    GrB_Matrix_ncols(&dimension, lhs_mat) ;
    GrB_Matrix res ;
    GrB_Matrix_new(&res, GrB_BOOL, dimension, dimension) ;
    GRB_TRY(GrB_eWiseAdd(res, GrB_NULL, GrB_NULL,
                         GxB_ANY_BOOL, lhs_mat, rhs_mat, GrB_DESC_R)) ;
    plan->res_mat = res ;

    return (GrB_SUCCESS) ;
}

static GrB_Info LAGraph_RPQMatrixConcat(RPQMatrixPlan *plan, char *msg)
{

    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG(plan->op == RPQ_MATRIX_OP_CONCAT, GrB_INVALID_VALUE, "operator is not concat") ;
    LG_ASSERT_MSG(plan->res_mat == NULL, GrB_INVALID_VALUE, "resulting matrix is already set as concat result") ;

    RPQMatrixPlan *lhs = plan->lhs ;
    RPQMatrixPlan *rhs = plan->rhs ;

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER) ;

    OK(LAGraph_RPQMatrix_solver(lhs, msg)) ;
    OK(LAGraph_RPQMatrix_solver(rhs, msg)) ;

    GrB_Matrix lhs_mat = lhs->res_mat ;
    GrB_Matrix rhs_mat = rhs->res_mat ;

    GrB_Index dimension ;
    GrB_Matrix_ncols(&dimension, lhs_mat) ;
    GrB_Matrix res ;
    GrB_Matrix_new(&res, GrB_BOOL, dimension, dimension) ;
    GRB_TRY(GrB_mxm(res, GrB_NULL, GrB_NULL,
                    sr, lhs_mat, rhs_mat, GrB_DESC_R)) ;
    plan->res_mat = res ;
    return (GrB_SUCCESS) ;
}

// A hypersparse seed that is much smaller than an ultra-sparse repeated step
// is normally a bound endpoint.  For these matrices the mask setup costs more
// than the redundant work avoided by the frontier algorithm.
#define RPQ_MATRIX_NAIVE_SEED_RATIO 64
#define RPQ_MATRIX_NAIVE_STEP_DEGREE_DENOMINATOR 64

static GrB_Info LAGraph_RPQMatrixUseNaiveClosure(
    bool *use_naive,
    GrB_Matrix seed,
    GrB_Matrix step,
    char *msg)
{
    LG_ASSERT(use_naive != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(seed != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(step != NULL, GrB_NULL_POINTER) ;

    GrB_Index n = 0, seed_nnz = 0, step_nnz = 0 ;
    int seed_sparsity = LG_SPARSE, step_sparsity = LG_SPARSE ;
    GRB_TRY(GrB_Matrix_nrows(&n, step)) ;
    GRB_TRY(GrB_Matrix_nvals(&seed_nnz, seed)) ;
    GRB_TRY(GrB_Matrix_nvals(&step_nnz, step)) ;
    GRB_TRY(LG_GET_FORMAT_HINT(seed, &seed_sparsity)) ;
    GRB_TRY(LG_GET_FORMAT_HINT(step, &step_sparsity)) ;

    bool step_is_sparse = step_sparsity == LG_HYPERSPARSE ||
        step_sparsity == LG_SPARSE ;
    bool step_is_ultrasparse = step_nnz <=
        n / RPQ_MATRIX_NAIVE_STEP_DEGREE_DENOMINATOR ;
    *use_naive = seed_sparsity == LG_HYPERSPARSE && step_is_sparse &&
        step_is_ultrasparse && seed_nnz > 0 &&
        seed_nnz <= step_nnz / RPQ_MATRIX_NAIVE_SEED_RATIO ;
    return (GrB_SUCCESS) ;
}

// Compute the fixed point from the accumulated result.  This avoids mask
// overhead for very small hypersparse seeds.
static GrB_Info LAGraph_RPQMatrixNaiveClosure(
    GrB_Matrix *result,
    GrB_Matrix seed,
    GrB_Matrix step,
    bool step_on_left,
    char *msg)
{
    LG_ASSERT(result != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(seed != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(step != NULL, GrB_NULL_POINTER) ;

    GrB_Matrix S = GrB_NULL ;
    GrB_Matrix T = GrB_NULL ;
    GrB_Matrix U = GrB_NULL ;
    GRB_TRY(GrB_Matrix_dup(&S, seed)) ;

    GrB_Index n = 0, nnz = 0, previous_nnz = 0, step_nnz = 0 ;
    GRB_TRY(GrB_Matrix_nrows(&n, seed)) ;
    GRB_TRY(GrB_Matrix_nvals(&previous_nnz, S)) ;
    GRB_TRY(GrB_Matrix_nvals(&step_nnz, step)) ;

    while (previous_nnz > 0 && step_nnz > 0)
    {
        GRB_TRY(GrB_Matrix_new(&T, GrB_BOOL, n, n)) ;
        if (step_on_left)
        {
            GRB_TRY(GrB_mxm(T, GrB_NULL, GrB_NULL, sr,
                            step, S, GrB_NULL)) ;
        }
        else
        {
            GRB_TRY(GrB_mxm(T, GrB_NULL, GrB_NULL, sr,
                            S, step, GrB_NULL)) ;
        }

        GRB_TRY(GrB_Matrix_new(&U, GrB_BOOL, n, n)) ;
        GRB_TRY(GrB_eWiseAdd(U, GrB_NULL, GrB_NULL,
                             GxB_ANY_BOOL, S, T, GrB_NULL)) ;
        GRB_TRY(GrB_Matrix_free(&T)) ;
        GRB_TRY(GrB_Matrix_nvals(&nnz, U)) ;

        GRB_TRY(GrB_Matrix_free(&S)) ;
        S = U ;
        U = GrB_NULL ;
        if (nnz == previous_nnz)
        {
            break ;
        }
        previous_nnz = nnz ;
    }

    *result = S ;
    return (GrB_SUCCESS) ;
}

// Compute the least fixed point starting from seed.  S keeps all discovered
// pairs, while frontier contains only pairs discovered by the previous step.
static GrB_Info LAGraph_RPQMatrixFrontierClosure(
    GrB_Matrix *result,
    GrB_Matrix seed,
    GrB_Matrix step,
    bool step_on_left,
    char *msg)
{
    LG_ASSERT(result != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(seed != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(step != NULL, GrB_NULL_POINTER) ;

    GrB_Matrix S = GrB_NULL ;
    GrB_Matrix frontier = GrB_NULL ;
    GRB_TRY(GrB_Matrix_dup(&S, seed)) ;
    GRB_TRY(GrB_Matrix_dup(&frontier, seed)) ;

    GrB_Index frontier_nnz = 0, step_nnz = 0 ;
    GRB_TRY(GrB_Matrix_nvals(&frontier_nnz, frontier)) ;
    GRB_TRY(GrB_Matrix_nvals(&step_nnz, step)) ;

    while (frontier_nnz > 0 && step_nnz > 0)
    {
        // frontier<!S> = step * frontier or frontier * step.  The replace,
        // structural, complemented mask leaves only newly discovered pairs.
        if (step_on_left)
        {
            GRB_TRY(GrB_mxm(frontier, S, GrB_NULL, sr,
                            step, frontier, GrB_DESC_RSC)) ;
        }
        else
        {
            GRB_TRY(GrB_mxm(frontier, S, GrB_NULL, sr,
                            frontier, step, GrB_DESC_RSC)) ;
        }

        GRB_TRY(GrB_Matrix_nvals(&frontier_nnz, frontier)) ;
        GRB_TRY(GrB_eWiseAdd(S, GrB_NULL, GrB_NULL,
                             GxB_ANY_BOOL, S, frontier, GrB_NULL)) ;
    }

    GRB_TRY(GrB_Matrix_free(&frontier)) ;
    *result = S ;
    return (GrB_SUCCESS) ;
}

static GrB_Info LAGraph_RPQMatrixAdaptiveClosure(
    GrB_Matrix *result,
    GrB_Matrix seed,
    GrB_Matrix step,
    bool step_on_left,
    char *msg)
{
    bool use_naive = false ;
    GRB_TRY(LAGraph_RPQMatrixUseNaiveClosure(
        &use_naive, seed, step, msg)) ;
    if (use_naive)
    {
        return LAGraph_RPQMatrixNaiveClosure(
            result, seed, step, step_on_left, msg) ;
    }
    return LAGraph_RPQMatrixFrontierClosure(
        result, seed, step, step_on_left, msg) ;
}

static GrB_Info LAGraph_RPQMatrixKleene(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG(plan->op == RPQ_MATRIX_OP_KLEENE, GrB_INVALID_VALUE, "operator is not kleene") ;
    LG_ASSERT_MSG(plan->res_mat == NULL, GrB_INVALID_VALUE, "resulting matrix is already set for kleene") ;

    RPQMatrixPlan *lhs = plan->lhs ;
    RPQMatrixPlan *rhs = plan->rhs ;

    // Kleene star should have one child. Always right.
    LG_ASSERT_MSG(lhs == NULL, GrB_INVALID_VALUE, "lhs is expected to be NULL for kleene") ;
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER) ;

    OK(LAGraph_RPQMatrix_solver(rhs, msg)) ;

    GrB_Matrix B = rhs->res_mat ;
    GrB_Matrix S = GrB_NULL ;
    GrB_Matrix U = GrB_NULL ;

    GrB_Index n ;
    GRB_TRY(GrB_Matrix_nrows(&n, B)) ;

    GRB_TRY(LAGraph_RPQMatrixAdaptiveClosure(&S, B, B, false, msg)) ;

    // Add I to recieve B* from B+.
    GrB_Vector v = GrB_NULL ;
    GrB_Matrix I = GrB_NULL ;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n)) ;
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL)) ;
    GRB_TRY(GrB_Matrix_diag(&I, v, 0)) ;
    GRB_TRY(GrB_Vector_free(&v)) ;

    GRB_TRY(GrB_Matrix_new(&U, GrB_BOOL, n, n)) ;
    GRB_TRY(GrB_eWiseAdd(U, GrB_NULL, GrB_NULL, GxB_ANY_BOOL, S, I, GrB_NULL)) ;
    GRB_TRY(GrB_Matrix_free(&S)) ;
    GRB_TRY(GrB_Matrix_free(&I)) ;

    plan->res_mat = U ;
    return (GrB_SUCCESS) ;
}

// this function need to handle special case where some optimization
// are available.
//
// consider following AST:
//    ┌─┐
//    │/|
//    └┬┘
// ┌─┬─┴─┬─┐
// │*│   │b│
// └┬┘   └─┘
// ┌┴┐
// │a│
// └─┘
// If matrix B is sparse and A is dense, then instead of naive
// way:
//
// (I + A + A x A + ...) x B
//
// we can do:
//
// (B + A x B + A x A x B + ...)
//
// and AST should be rewritten in the following way:
//   ┌───┐
//   │L^*│
//   └─┬─┘
// ┌─┬─┴─┬─┐
// │a│   │b│
// └─┘   └─┘
static GrB_Info LAGraph_RPQMatrixKleene_L(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG(plan->op == RPQ_MATRIX_OP_KLEENE_L, GrB_INVALID_VALUE, "different operator is not expected left-kleene") ;
    LG_ASSERT_MSG(plan->res_mat == NULL, GrB_INVALID_VALUE, "resulting matrix is already set for left-kleene") ;

    RPQMatrixPlan *lhs = plan->lhs ; // A
    RPQMatrixPlan *rhs = plan->rhs ; // B

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER) ;

    OK(LAGraph_RPQMatrix_solver(lhs, msg)) ;
    OK(LAGraph_RPQMatrix_solver(rhs, msg)) ;

    GrB_Matrix A = lhs->res_mat ;
    GrB_Matrix B = rhs->res_mat ;

    GrB_Matrix S = GrB_NULL ;
    GRB_TRY(LAGraph_RPQMatrixAdaptiveClosure(&S, B, A, true, msg)) ;

    plan->res_mat = S ;
    return GrB_SUCCESS ;
}
// this function need to handle special case where some optimization
// are available.
// consider following AST:
//    ┌─┐
//    │/|
//    └┬┘
// ┌─┬─┴─┬─┐
// │a│   │*│
// └─┘   └┬┘
//       ┌┴┐
//       │b│
//       └─┘
// If matrix A is sparse and B is dense, then instead of naive
// way:
//
// A x (I + B + B x B + ...)
//
// we can do:
//
// (A + A x B + A x B x B + ...)
//
// and AST should be rewritten in the following way:
//   ┌───┐
//   │R^*│
//   └─┬─┘
// ┌─┬─┴─┬─┐
// │a│   │b│
// └─┘   └─┘

static GrB_Info LAGraph_RPQMatrixKleene_R(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG(plan->op == RPQ_MATRIX_OP_KLEENE_R, GrB_INVALID_VALUE, "different operator is not expected right-kleene") ;
    LG_ASSERT_MSG(plan->res_mat == NULL, GrB_INVALID_VALUE, "resulting matrix is already set for right-kleene") ;

    RPQMatrixPlan *lhs = plan->lhs ; // A
    RPQMatrixPlan *rhs = plan->rhs ; // B

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER) ;

    OK(LAGraph_RPQMatrix_solver(lhs, msg)) ;
    OK(LAGraph_RPQMatrix_solver(rhs, msg)) ;

    GrB_Matrix A = lhs->res_mat ;
    GrB_Matrix B = rhs->res_mat ;

    GrB_Matrix S = GrB_NULL ;
    GRB_TRY(LAGraph_RPQMatrixAdaptiveClosure(&S, A, B, false, msg)) ;

    plan->res_mat = S ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_solver(RPQMatrixPlan *plan, char *msg)
{
    if (plan->res_mat != NULL)
    {
        return (GrB_SUCCESS) ;
    }

    switch (plan->op)
    {
    case RPQ_MATRIX_OP_LABEL:
        LG_ASSERT_MSG(plan->lhs == NULL && plan->rhs == NULL,
                      GrB_INVALID_VALUE, "label node should not have any children nodes") ;
        plan->res_mat = plan->mat ;
        return (GrB_SUCCESS) ;
    case RPQ_MATRIX_OP_LOR:
        return LAGraph_RPQMatrixLor(plan, msg) ;
    case RPQ_MATRIX_OP_CONCAT:
        return LAGraph_RPQMatrixConcat(plan, msg) ;
    case RPQ_MATRIX_OP_KLEENE:
        return LAGraph_RPQMatrixKleene(plan, msg) ;
    case RPQ_MATRIX_OP_KLEENE_L:
        return LAGraph_RPQMatrixKleene_L(plan, msg) ;
    case RPQ_MATRIX_OP_KLEENE_R:
        return LAGraph_RPQMatrixKleene_R(plan, msg) ;
    default:
        LG_ASSERT_MSG(false, GrB_INVALID_VALUE, "invalid graph node type") ;
    }
    return (GrB_SUCCESS) ;
}

GrB_Info LAGraph_RPQMatrix_initialize(void)
{
    if (sr != GrB_NULL)
    {
        return GrB_SUCCESS ;
    }
    sr = GxB_ANY_PAIR_BOOL ;
    op = GxB_ANY_BOOL_MONOID ;
    srand(time(NULL)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix(
    // output:
    GrB_Index *nnz, // number of nonzero values in
                    // result reachability matrix

    // input:
    RPQMatrixPlan *plan, // root of abstarct syntax tree of
                         // regular expression
    char *msg            // LAGraph output message
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER) ;
    GrB_Index dimension = -1 ;
    GrB_Info info = LAGraph_RPQMatrix_check(plan, &dimension, msg) ;
    LG_ASSERT_MSG(info == GrB_SUCCESS, info, msg) ;

    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------

    LAGraph_RPQMatrix_initialize() ;

    //--------------------------------------------------------------------------
    // run solver
    //--------------------------------------------------------------------------

    info = LAGraph_RPQMatrix_solver(plan, msg) ;
    LG_ASSERT_MSG(info == GrB_SUCCESS, info, msg) ;
    GrB_Matrix_nvals(nnz, plan->res_mat) ;
    return GrB_SUCCESS ;
}
