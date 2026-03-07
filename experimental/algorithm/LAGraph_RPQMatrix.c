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
                         GrB_LOR, lhs_mat, rhs_mat, GrB_DESC_R)) ;
    plan->res_mat = res ;

    return (GrB_SUCCESS) ;
}

#define abs(a) ((a) > 0.0f ? (a) : (-(a)))

double error = 0.0 ;
double error2 = 0.0 ;
double c = 0.0 ;

GrB_Info LAGraph_RPQMatrixEstimate_WanderJoin(GrB_Index *estimate, GrB_Matrix lhs, GrB_Matrix rhs) ;
GrB_Info LAGraph_RPQMatrixEstimate_WanderJoin_Kleene(GrB_Index *estimate, GrB_Matrix rhs) ;

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
    // S <- I
    GrB_Matrix S ;

    // Creating identity matrix.
    GrB_Index n ;
    GRB_TRY(GrB_Matrix_nrows(&n, B)) ;
    // GrB_Matrix I ;
    // GRB_TRY(GrB_Matrix_new(&I, GrB_BOOL, n, n)) ;

    GrB_Vector v ;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n)) ;
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL)) ;

    GRB_TRY(GrB_Matrix_diag(&S, v, 0)) ;

    bool changed = true ;
    GrB_Index nnz_S = n, nnz_Sold = 0 ;

    while (changed)
    {
        // S <- S x (B + I)
        GRB_TRY(GrB_mxm(S, S, GrB_NULL,
                        sr, S, B, GrB_DESC_C)) ;

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S)) ;
        if (nnz_S != nnz_Sold)
        {
            changed = true ;
            nnz_Sold = nnz_S ;
        }
        else
        {
            changed = false ;
        }
    }
    GrB_Vector_free(&v) ;
    plan->res_mat = S ;

    // GRB_TRY(GrB_Matrix_free(&I)) ;
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

    // S <- B
    GrB_Matrix S ;
    GRB_TRY(GrB_Matrix_dup(&S, B)) ;

    bool changed = true ;
    GrB_Index nnz_S = 0, nnz_Sold = 0 ;

    while (changed)
    {
        // S <- (A + I) x S
        GRB_TRY(GrB_mxm(S, S, NULL, sr, A, S, GrB_DESC_C)) ;

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S)) ;
        if (nnz_S != nnz_Sold)
        {
            changed = true ;
            nnz_Sold = nnz_S ;
        }
        else
        {
            changed = false ;
        }
    }

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

    // S <- A
    GrB_Matrix S ;
    GRB_TRY(GrB_Matrix_dup(&S, A)) ;

    bool changed = true ;
    GrB_Index nnz_S = 0, nnz_Sold = 0 ;

    while (changed)
    {
        // S <- S x (B + I)
        GRB_TRY(GrB_mxm(S, S, NULL, sr, S, B, GrB_DESC_C)) ;

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S)) ;
        if (nnz_S != nnz_Sold)
        {
            changed = true ;
            nnz_Sold = nnz_S ;
        }
        else
        {
            changed = false ;
        }
    }

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
    sr = LAGraph_any_one_bool ;
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

GrB_Info LAGraph_RPQMatrixOuts (GrB_Vector w, const GrB_Matrix A) {
    return GrB_reduce (w, GrB_NULL, GrB_NULL, op, A, GrB_DESC_T0) ;
}
GrB_Info LAGraph_RPQMatrixIns (GrB_Vector w, const GrB_Matrix A) {
    return GrB_reduce (w, GrB_NULL, GrB_NULL, op, A, GrB_NULL) ;
}

#define max(a, b) ((a) >= (b) ? (a) : (b))
#define min(a, b) ((a) <= (b) ? (a) : (b))

GrB_Info LAGraph_RPQMatrixEstimate (float *estimate,float *estimate2,  GrB_Vector outs, GrB_Vector ins) {
    GrB_Index outc ;
    GrB_Index inc ;

    OK (GrB_Vector_nvals(&outc, outs)) ;
    OK (GrB_Vector_nvals(&inc, ins)) ;
    //printf ("OUTC INC: %llu %llu\n\n", outc, inc) ;

    GrB_Index dimension ;
    GrB_Vector_size(&dimension, ins) ;

    GrB_Vector inter;
    OK (GrB_Vector_new (&inter, GrB_BOOL, dimension)) ;
    GrB_Index interc ;

    OK (GrB_assign(inter, ins, GrB_NULL, outs, GrB_ALL, dimension, GrB_DESC_RS));
    //GxB_print(inter, 2);

    OK (GrB_Vector_nvals(&interc, inter)) ;

    *estimate = ((float) interc) / ((float) max(outc, inc)) / ((float) min(outc, inc));
    *estimate2 = 1 / ((float) max(outc, inc)) ;

    return GrB_SUCCESS ;
}

#define WANDER_JOIN_FACTOR 8

int compare( const void* a, const void* b)
{
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );
     
     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}


GrB_Info LAGraph_RPQMatrixEstimate_WanderJoin (GrB_Index *estimate, GrB_Matrix lhs, GrB_Matrix rhs) {
    // TODO: check if square, of the same dimensions.

    GrB_Index dim ;
    GRB_TRY (GrB_Matrix_nrows(&dim, lhs)) ;
    GrB_Index nvals ;
    GRB_TRY (GrB_Matrix_nvals(&nvals, lhs)) ;

    // nvals = n * n ~> join_factor = n
    // nvals = big ~> join_factor = big
    // nvals = small ~> join_factor = small
    // nvals = 0 ~> join_factor = 1
    GrB_Index join_factor = WANDER_JOIN_FACTOR ;

    GrB_Index sdim = dim / join_factor ;
    GrB_Matrix slhs ;
    GrB_Matrix srhs ;
    GRB_TRY (GrB_Matrix_new(&slhs, GrB_BOOL, sdim, sdim)) ;
    GRB_TRY (GrB_Matrix_new(&srhs, GrB_BOOL, sdim, sdim)) ;

    GrB_Index I[sdim] ;
    for (size_t i = 0; i < sdim; i++) {
        GrB_Index k ;
        while (true) {
            k = rand() % dim ;
            for (size_t j = 0; j < i; j++) {
                if (I[j] == k) {
                    goto skip;
                }
            }
            break;
            skip:;
        }
        I[i] = k ;
    }

    GRB_TRY (GrB_extract(slhs, GrB_NULL, GrB_NULL, lhs, I, sdim, I, sdim, GrB_DESC_R)) ;
    GRB_TRY (GrB_extract(srhs, GrB_NULL, GrB_NULL, rhs, I, sdim, I, sdim, GrB_DESC_R)) ;

    GrB_Matrix sres ;
    GRB_TRY (GrB_Matrix_new(&sres, GrB_BOOL, sdim, sdim)) ;

    GRB_TRY (GrB_mxm(sres, GrB_NULL, GrB_NULL, sr, slhs, srhs, GrB_DESC_R)) ;

    GrB_Index snv ;
    GRB_TRY (GrB_Matrix_nvals(&snv, sres)) ;

    *estimate = snv * join_factor * join_factor * join_factor ;

    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrixEstimate_WanderJoin_Kleene (GrB_Index *estimate, GrB_Matrix rhs) {
    // TODO: check if square, of the same dimensions.

    GrB_Index dim ;

    GRB_TRY (GrB_Matrix_nrows(&dim, rhs)) ;

    GrB_Index sdim = dim / WANDER_JOIN_FACTOR ;
    GrB_Matrix srhs ;
    GRB_TRY (GrB_Matrix_new(&srhs, GrB_BOOL, sdim, sdim)) ;

    GrB_Index I[sdim] ;
    for (size_t i = 0; i < sdim; i++) {
        GrB_Index k ;
        while (true) {
            k = rand() % dim ;
            for (size_t j = 0; j < i; j++) {
                if (I[j] == k) {
                    goto skip;
                }
            }
            break;
            skip:;
        }
        I[i] = k ;
    }

    GRB_TRY (GrB_extract(srhs, GrB_NULL, GrB_NULL, rhs, I, sdim, I, sdim, GrB_DESC_R)) ;

    GrB_Matrix S ;
    GRB_TRY(GrB_Matrix_dup(&S, srhs)) ;
    
    bool changed = true ;
    GrB_Index nnz_S = 0, nnz_Sold = 0 ;

    while (changed)
    {
        // S <- S x (B + I)
        GRB_TRY(GrB_mxm(S, S, GrB_NULL,
                        sr, S, srhs, GrB_DESC_C)) ;

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S)) ;
        if (nnz_S != nnz_Sold)
        {
            changed = true ;
            nnz_Sold = nnz_S ;
        }
        else
        {
            changed = false ;
        }
    }

    GrB_Index snv ;
    GRB_TRY (GrB_Matrix_nvals(&snv, S)) ;

    *estimate = snv * WANDER_JOIN_FACTOR * WANDER_JOIN_FACTOR * WANDER_JOIN_FACTOR;

    return GrB_SUCCESS ;
}

static GrB_Index *I = NULL;

GrB_Info LAGraph_RPQMatrix_Alt (GrB_Matrix lhs, GrB_Matrix rhs, GrB_Matrix *res, uint64_t *nvals) {
    GrB_Index sdim;
    GRB_TRY (GrB_Matrix_nrows(&sdim, rhs)) ;

    GRB_TRY (GrB_Matrix_new(res, GrB_BOOL, sdim, sdim)) ;

    GRB_TRY(GrB_eWiseAdd(*res, GrB_NULL, GrB_NULL,
                         op, lhs, rhs, GrB_DESC_R)) ;

    GRB_TRY (GrB_Matrix_nvals(nvals, *res)) ;
}
GrB_Info LAGraph_RPQMatrix_Seq (GrB_Matrix lhs, GrB_Matrix rhs, GrB_Matrix *res, uint64_t *nvals) {
    GrB_Index sdim;
    GRB_TRY (GrB_Matrix_nrows(&sdim, rhs)) ;

    GRB_TRY (GrB_Matrix_new(res, GrB_BOOL, sdim, sdim)) ;

    GRB_TRY (GrB_mxm(*res, GrB_NULL, GrB_NULL, sr, lhs, rhs, GrB_DESC_R)) ;

    GRB_TRY (GrB_Matrix_nvals(nvals, *res)) ;
}
GrB_Info LAGraph_RPQMatrix_ExtractRandom (GrB_Matrix rhs, GrB_Matrix *srhs, uint64_t seed) {
    // TODO: check if square, of the same dimensions.
    GrB_Index dim ;

    GRB_TRY (GrB_Matrix_nrows(&dim, rhs)) ;

    GrB_Index sdim = dim / WANDER_JOIN_FACTOR ;
    GRB_TRY (GrB_Matrix_new(srhs, GrB_BOOL, sdim, sdim)) ;


    if (I == NULL) {
        srand(seed);
        I = calloc(sizeof(GrB_Index), sdim);

        for (size_t i = 0; i < sdim; i++) {
            size_t u = rand();
            size_t k = u % WANDER_JOIN_FACTOR + WANDER_JOIN_FACTOR * i;
            I[i] = k ;
        }
    }

    GRB_TRY (GrB_extract(*srhs, GrB_NULL, GrB_NULL, rhs, I, sdim, I, sdim, GrB_DESC_R)) ;
    return GrB_SUCCESS ;
}

/*GrB_Info LAGraph_RPQMatrixEstimateSum(GrB_Index *estimate, RPQMatrix *lhs, RPQMatrix *rhs, char *msg) {
    GrB_Index lhs_nvals, rhs_nvals ;
    OK (GrB_Matrix_nvals(&lhs_nvals, lhs->m)) ;
    OK (GrB_Matrix_nvals(&rhs_nvals, lhs->m)) ;

    *estimate = lhs_nvals + rhs_nvals ;

    return GrB_SUCCESS ;
}*/
