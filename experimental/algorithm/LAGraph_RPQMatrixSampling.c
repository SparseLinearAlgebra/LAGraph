#include "LG_internal.h"
#include "LAGraphX.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

static char msg[LAGRAPH_MSG_LEN] ;

GrB_Info LAGraph_RPQMatrix_sample_submatrix(GrB_Matrix *result, GrB_Matrix source, const GrB_Index *vertices, GrB_Index count)
{
    LG_ASSERT(result != NULL && source != GrB_NULL && (vertices != NULL || count == 0), GrB_NULL_POINTER) ;
    *result = GrB_NULL ;
    GrB_Index rows, cols ;
    OK(GrB_Matrix_nrows(&rows, source)) ;
    OK(GrB_Matrix_ncols(&cols, source)) ;
    if (rows != cols)
    {
        return GrB_DIMENSION_MISMATCH ;
    }
    OK(GrB_Matrix_new(result, GrB_BOOL, count, count)) ;
    GrB_Info info = GrB_Matrix_extract(*result, GrB_NULL, GrB_NULL, source, vertices, count, vertices, count, GrB_NULL) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free(result) ;
    }
    return info ;
}

GrB_Info LAGraph_RPQMatrix_sample_identity(GrB_Matrix *result, GrB_Index n)
{
    LG_ASSERT(result != NULL, GrB_NULL_POINTER) ;
    *result = GrB_NULL ;
    GrB_Vector diagonal = GrB_NULL ;
    OK(GrB_Vector_new(&diagonal, GrB_BOOL, n)) ;
    GrB_Info info = GrB_Vector_assign_BOOL(diagonal, GrB_NULL, GrB_NULL, true, GrB_ALL, n, GrB_NULL) ;
    if (info == GrB_SUCCESS)
    {
        info = GrB_Matrix_diag(result, diagonal, 0) ;
    }
    GrB_Vector_free(&diagonal) ;
    return info ;
}

GrB_Info LAGraph_RPQMatrix_sample_apply(GrB_Matrix *result, GrB_Matrix lhs, GrB_Matrix rhs)
{
    LG_ASSERT(result != NULL && lhs != GrB_NULL && rhs != GrB_NULL, GrB_NULL_POINTER) ;
    *result = GrB_NULL ;
    GrB_Index lhs_rows, lhs_cols, rhs_rows, rhs_cols ;
    OK(GrB_Matrix_nrows(&lhs_rows, lhs)) ;
    OK(GrB_Matrix_ncols(&lhs_cols, lhs)) ;
    OK(GrB_Matrix_nrows(&rhs_rows, rhs)) ;
    OK(GrB_Matrix_ncols(&rhs_cols, rhs)) ;
    if (lhs_cols != rhs_rows)
    {
        return GrB_DIMENSION_MISMATCH ;
    }
    OK(GrB_Matrix_new(result, GrB_BOOL, lhs_rows, rhs_cols)) ;
    GrB_Info info = GrB_mxm(*result, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, lhs, rhs, GrB_NULL) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free(result) ;
    }
    return info ;
}

GrB_Info LAGraph_RPQMatrix_sample_union(GrB_Matrix *result, GrB_Matrix lhs, GrB_Matrix rhs)
{
    LG_ASSERT(result != NULL && lhs != GrB_NULL && rhs != GrB_NULL, GrB_NULL_POINTER) ;
    *result = GrB_NULL ;
    GrB_Index lhs_rows, lhs_cols, rhs_rows, rhs_cols ;
    OK(GrB_Matrix_nrows(&lhs_rows, lhs)) ;
    OK(GrB_Matrix_ncols(&lhs_cols, lhs)) ;
    OK(GrB_Matrix_nrows(&rhs_rows, rhs)) ;
    OK(GrB_Matrix_ncols(&rhs_cols, rhs)) ;
    if (lhs_rows != rhs_rows || lhs_cols != rhs_cols)
    {
        return GrB_DIMENSION_MISMATCH ;
    }
    OK(GrB_Matrix_new(result, GrB_BOOL, lhs_rows, lhs_cols)) ;
    GrB_Info info = GrB_eWiseAdd(*result, GrB_NULL, GrB_NULL, GxB_ANY_BOOL, lhs, rhs, GrB_NULL) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free(result) ;
    }
    return info ;
}

GrB_Info LAGraph_RPQMatrix_sample_stats(GrB_Index *nvals, GrB_Index *active_rows, GrB_Index *active_cols, GrB_Index *diagonal_nvals, GrB_Matrix sample)
{
    LG_ASSERT(nvals != NULL && active_rows != NULL && active_cols != NULL && diagonal_nvals != NULL && sample != GrB_NULL, GrB_NULL_POINTER) ;
    OK(GrB_Matrix_nvals(nvals, sample)) ;
    OK(LAGraph_RPQMatrix_reduce(active_rows, sample, 0)) ;
    OK(LAGraph_RPQMatrix_reduce(active_cols, sample, 1)) ;
    GrB_Index n ;
    GrB_Vector diagonal = GrB_NULL ;
    OK(GrB_Matrix_nrows(&n, sample)) ;
    OK(GrB_Vector_new(&diagonal, GrB_BOOL, n)) ;
    GrB_Info info = GxB_Vector_diag(diagonal, sample, 0, GrB_NULL) ;
    if (info == GrB_SUCCESS)
    {
        info = GrB_Vector_nvals(diagonal_nvals, diagonal) ;
    }
    GrB_Vector_free(&diagonal) ;
    return info ;
}

GrB_Info LAGraph_RPQMatrix_reduce_count_vector(GrB_Vector *res, GrB_Matrix mat, uint8_t reduce_type)
{
    LG_ASSERT(res != NULL && mat != GrB_NULL, GrB_NULL_POINTER) ;
    if (reduce_type > 1)
    {
        return GrB_INVALID_VALUE ;
    }
    GrB_Index rows, cols ;
    OK(GrB_Matrix_nrows(&rows, mat)) ;
    OK(GrB_Matrix_ncols(&cols, mat)) ;
    OK(GrB_Vector_new(res, GrB_UINT64, reduce_type == 0 ? rows : cols)) ;
    GrB_Info info = GrB_reduce(*res, GrB_NULL, GrB_NULL, GrB_PLUS_MONOID_UINT64, mat, reduce_type == 0 ? GrB_NULL : GrB_DESC_T0) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Vector_free(res) ;
    }
    return info ;
}

static GrB_Info LAGraph_RPQMatrix_extended_count_vectors_work(GrB_Vector *row_extended, GrB_Vector *col_extended, GrB_Vector *singleton_rows, GrB_Vector *singleton_cols, GrB_Matrix mat, GrB_Vector row_counts, GrB_Vector col_counts, GrB_Index rows, GrB_Index cols)
{
    OK(GrB_Vector_new(singleton_rows, GrB_BOOL, rows)) ;
    OK(GrB_Vector_new(singleton_cols, GrB_BOOL, cols)) ;
    OK(GrB_select(*singleton_rows, GrB_NULL, GrB_NULL, GrB_VALUEEQ_UINT64, row_counts, 1UL, GrB_NULL)) ;
    OK(GrB_select(*singleton_cols, GrB_NULL, GrB_NULL, GrB_VALUEEQ_UINT64, col_counts, 1UL, GrB_NULL)) ;
    OK(GrB_Vector_new(row_extended, GrB_UINT64, rows)) ;
    OK(GrB_mxv(*row_extended, GrB_NULL, GrB_NULL, GxB_PLUS_PAIR_UINT64, mat, *singleton_cols, GrB_NULL)) ;
    OK(GrB_Vector_new(col_extended, GrB_UINT64, cols)) ;
    OK(GrB_vxm(*col_extended, GrB_NULL, GrB_NULL, GxB_PLUS_PAIR_UINT64, *singleton_rows, mat, GrB_NULL)) ;
    return GrB_SUCCESS ;
}

// h^er counts entries in singleton columns; h^ec counts entries in singleton rows.
GrB_Info LAGraph_RPQMatrix_extended_count_vectors(GrB_Vector *row_extended, GrB_Vector *col_extended, GrB_Matrix mat, GrB_Vector row_counts, GrB_Vector col_counts)
{
    LG_ASSERT(row_extended != NULL && col_extended != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(row_extended != col_extended, GrB_INVALID_VALUE) ;
    LG_ASSERT(mat != GrB_NULL && row_counts != GrB_NULL && col_counts != GrB_NULL, GrB_NULL_POINTER) ;

    *row_extended = GrB_NULL ;
    *col_extended = GrB_NULL ;
    GrB_Vector singleton_rows = GrB_NULL, singleton_cols = GrB_NULL ;
    GrB_Index rows, cols, row_size, col_size ;
    OK(GrB_Matrix_nrows(&rows, mat)) ;
    OK(GrB_Matrix_ncols(&cols, mat)) ;
    OK(GrB_Vector_size(&row_size, row_counts)) ;
    OK(GrB_Vector_size(&col_size, col_counts)) ;
    if (rows != row_size || cols != col_size)
    {
        return GrB_DIMENSION_MISMATCH ;
    }

    GrB_Info info = LAGraph_RPQMatrix_extended_count_vectors_work(row_extended, col_extended, &singleton_rows, &singleton_cols, mat, row_counts, col_counts, rows, cols) ;
    GrB_Vector_free(&singleton_rows) ;
    GrB_Vector_free(&singleton_cols) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Vector_free(row_extended) ;
        GrB_Vector_free(col_extended) ;
    }
    return info ;
}

GrB_Info LAGraph_RPQMatrix_count_vector_dot(double *res, GrB_Vector lhs, GrB_Vector rhs)
{
    LG_ASSERT(res != NULL && lhs != GrB_NULL && rhs != GrB_NULL, GrB_NULL_POINTER) ;
    GrB_Index ls, rs ;
    OK(GrB_Vector_size(&ls, lhs)) ;
    OK(GrB_Vector_size(&rs, rhs)) ;
    if (ls != rs)
    {
        return GrB_DIMENSION_MISMATCH ;
    }
    GrB_Vector product = GrB_NULL ;
    OK(GrB_Vector_new(&product, GrB_FP64, ls)) ;
    GrB_Info info = GrB_eWiseMult(product, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, lhs, rhs, GrB_NULL) ;
    if (info == GrB_SUCCESS)
    {
        info = GrB_reduce(res, GrB_NULL, GrB_PLUS_MONOID_FP64, product, GrB_NULL) ;
    }
    GrB_Vector_free(&product) ;
    return info ;
}

static GrB_Info LAGraph_RPQMatrix_count_vector_count_where(GrB_Index *res, GrB_Vector vector, double threshold, int equal)
{
    LG_ASSERT(res != NULL && vector != GrB_NULL, GrB_NULL_POINTER) ;
    GrB_Vector selected = GrB_NULL ;
    GrB_Index size ;
    GrB_Type type ;
    OK(GrB_Vector_size(&size, vector)) ;
    OK(GxB_Vector_type(&type, vector)) ;
    if (type != GrB_UINT64 && type != GrB_FP64)
    {
        return GrB_DOMAIN_MISMATCH ;
    }
    OK(GrB_Vector_new(&selected, GrB_BOOL, size)) ;
    GrB_Info info ;
    if (type == GrB_UINT64)
    {
        if (equal)
        {
            info = GrB_select(selected, GrB_NULL, GrB_NULL, GrB_VALUEEQ_UINT64, vector, (uint64_t) 1, GrB_NULL) ;
        }
        else
        {
            info = GrB_select(selected, GrB_NULL, GrB_NULL, GrB_VALUEGT_UINT64, vector, (uint64_t) threshold, GrB_NULL) ;
        }
    }
    else
    {
        info = GrB_select(selected, GrB_NULL, GrB_NULL, equal ? GrB_VALUEEQ_FP64 : GrB_VALUEGT_FP64, vector, equal ? 1.0 : threshold, GrB_NULL) ;
    }
    if (info == GrB_SUCCESS)
    {
        info = GrB_Vector_nvals(res, selected) ;
    }
    GrB_Vector_free(&selected) ;
    return info ;
}

static GrB_Info LAGraph_RPQMatrix_mnc_generic_product(GrB_Vector *product, GrB_Index *nvals, GrB_Vector lhs_cols, GrB_Vector rhs_rows, GrB_Index size)
{
    OK(GrB_Vector_new(product, GrB_FP64, size)) ;
    OK(GrB_eWiseMult(*product, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, lhs_cols, rhs_rows, GrB_NULL)) ;
    OK(GrB_Vector_nvals(nvals, *product)) ;
    return GrB_SUCCESS ;
}

static GrB_Info LAGraph_RPQMatrix_mnc_generic_nnz(double *res, GrB_Vector lhs_cols, GrB_Vector rhs_rows, double p)
{
    LG_ASSERT(res != NULL && lhs_cols != GrB_NULL && rhs_rows != GrB_NULL, GrB_NULL_POINTER) ;
    *res = 0.0 ;
    if (p <= 0.0)
    {
        return GrB_SUCCESS ;
    }

    GrB_Vector product = GrB_NULL ;
    GrB_Index lhs_n, product_nvals ;
    OK(GrB_Vector_size(&lhs_n, lhs_cols)) ;
    GrB_Info info = LAGraph_RPQMatrix_mnc_generic_product(&product, &product_nvals, lhs_cols, rhs_rows, lhs_n) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Vector_free(&product) ;
        return info ;
    }
    if (product_nvals == 0)
    {
        GrB_Vector_free(&product) ;
        return GrB_SUCCESS ;
    }

    double *values = malloc(product_nvals * sizeof(double)) ;
    if (values == NULL)
    {
        GrB_Vector_free(&product) ;
        return GrB_OUT_OF_MEMORY ;
    }
    GrB_Index extracted = product_nvals ;
    info = GrB_Vector_extractTuples_FP64(NULL, values, &extracted, product) ;
    GrB_Vector_free(&product) ;
    if (info != GrB_SUCCESS)
    {
        free(values) ;
        return info ;
    }

    double log_zero_probability = 0.0 ;
    for (GrB_Index i = 0 ; i < extracted ; i++)
    {
        double probability = fmin(fmax(values[i] / p, 0.0), 1.0) ;
        if (probability >= 1.0)
        {
            log_zero_probability = -INFINITY ;
            break ;
        }
        log_zero_probability += log1p(-probability) ;
    }
    *res = -p * expm1(log_zero_probability) ;
    free(values) ;
    return GrB_SUCCESS ;
}

static GrB_Info LAGraph_RPQMatrix_mnc_left_exact(double *exact_nnz, GrB_Index *remaining_rows, GrB_Vector *lhs_residual, GrB_Vector lhs_rows, GrB_Vector lhs_cols, GrB_Vector rhs_rows, GrB_Vector lhs_col_extended, GrB_Index size)
{
    GrB_Index singleton_rows ;
    OK(LAGraph_RPQMatrix_count_vector_count_where(&singleton_rows, lhs_rows, 1.0, 1)) ;
    *remaining_rows -= singleton_rows ;
    OK(LAGraph_RPQMatrix_count_vector_dot(exact_nnz, lhs_col_extended, rhs_rows)) ;
    OK(GrB_Vector_new(lhs_residual, GrB_FP64, size)) ;
    OK(GrB_eWiseAdd(*lhs_residual, GrB_NULL, GrB_NULL, GrB_MINUS_FP64, lhs_cols, lhs_col_extended, GrB_NULL)) ;
    return GrB_SUCCESS ;
}

static GrB_Info LAGraph_RPQMatrix_mnc_right_exact(double *exact_nnz, GrB_Index *remaining_cols, GrB_Vector *rhs_residual, GrB_Vector generic_lhs, GrB_Vector rhs_rows, GrB_Vector rhs_cols, GrB_Vector rhs_row_extended, GrB_Index size)
{
    GrB_Index singleton_cols ;
    OK(LAGraph_RPQMatrix_count_vector_count_where(&singleton_cols, rhs_cols, 1.0, 1)) ;
    *remaining_cols -= singleton_cols ;
    double additional_exact = 0.0 ;
    OK(LAGraph_RPQMatrix_count_vector_dot(&additional_exact, generic_lhs, rhs_row_extended)) ;
    *exact_nnz += additional_exact ;
    OK(GrB_Vector_new(rhs_residual, GrB_FP64, size)) ;
    OK(GrB_eWiseAdd(*rhs_residual, GrB_NULL, GrB_NULL, GrB_MINUS_FP64, rhs_rows, rhs_row_extended, GrB_NULL)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_count_vector_mnc_matmul_nnz(double *res, GrB_Vector lhs_rows, GrB_Vector lhs_cols, GrB_Vector rhs_rows, GrB_Vector rhs_cols, GrB_Vector lhs_col_extended, GrB_Vector rhs_row_extended)
{
    LG_ASSERT(res != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(lhs_rows != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(lhs_cols != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs_rows != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs_cols != GrB_NULL, GrB_NULL_POINTER) ;

    GrB_Index lhs_m ;
    GrB_Index lhs_n ;
    GrB_Index rhs_n ;
    GrB_Index rhs_l ;
    OK(GrB_Vector_size(&lhs_m, lhs_rows)) ;
    OK(GrB_Vector_size(&lhs_n, lhs_cols)) ;
    OK(GrB_Vector_size(&rhs_n, rhs_rows)) ;
    OK(GrB_Vector_size(&rhs_l, rhs_cols)) ;
    if (lhs_n != rhs_n)
    {
        return GrB_DIMENSION_MISMATCH ;
    }
    GrB_Index extended_size ;
    if (lhs_col_extended != GrB_NULL)
    {
        OK(GrB_Vector_size(&extended_size, lhs_col_extended)) ;
        if (extended_size != lhs_n)
        {
            return GrB_DIMENSION_MISMATCH ;
        }
    }
    if (rhs_row_extended != GrB_NULL)
    {
        OK(GrB_Vector_size(&extended_size, rhs_row_extended)) ;
        if (extended_size != lhs_n)
        {
            return GrB_DIMENSION_MISMATCH ;
        }
    }

    GrB_Type lhs_row_type, rhs_col_type ;
    OK(GxB_Vector_type(&lhs_row_type, lhs_rows)) ;
    OK(GxB_Vector_type(&rhs_col_type, rhs_cols)) ;
    double lhs_row_max = 0 ;
    double rhs_col_max = 0 ;
    if (lhs_row_type == GrB_UINT64)
    {
        OK(GrB_reduce(&lhs_row_max, GrB_NULL, GrB_MAX_MONOID_FP64, lhs_rows, GrB_NULL)) ;
    }
    if (rhs_col_type == GrB_UINT64)
    {
        OK(GrB_reduce(&rhs_col_max, GrB_NULL, GrB_MAX_MONOID_FP64, rhs_cols, GrB_NULL)) ;
    }

    double output_cells = (double) lhs_m * (double) rhs_l ;
    // Only UINT64 leaf counts certify the exact structural condition in Theorem 3.1.
    if ((lhs_row_type == GrB_UINT64 && lhs_row_max <= 1.0) || (rhs_col_type == GrB_UINT64 && rhs_col_max <= 1.0))
    {
        double dot = 0 ;
        OK(LAGraph_RPQMatrix_count_vector_dot(&dot, lhs_cols, rhs_rows)) ;
        *res = fmin(fmax(dot, 0.0), output_cells) ;
        return (GrB_SUCCESS) ;
    }

    GrB_Index lhs_nonempty_rows ;
    GrB_Index rhs_nonempty_cols ;
    OK(GrB_Vector_nvals(&lhs_nonempty_rows, lhs_rows)) ;
    OK(GrB_Vector_nvals(&rhs_nonempty_cols, rhs_cols)) ;
    GrB_Index remaining_rows = lhs_nonempty_rows ;
    GrB_Index remaining_cols = rhs_nonempty_cols ;
    double exact_nnz = 0.0, generic_nnz = 0.0 ;
    GrB_Vector lhs_residual = GrB_NULL, rhs_residual = GrB_NULL ;
    GrB_Vector generic_lhs = lhs_cols, generic_rhs = rhs_rows ;
    if (lhs_col_extended != GrB_NULL)
    {
        GrB_Info info = LAGraph_RPQMatrix_mnc_left_exact(&exact_nnz, &remaining_rows, &lhs_residual, lhs_rows, lhs_cols, rhs_rows, lhs_col_extended, lhs_n) ;
        if (info != GrB_SUCCESS)
        {
            GrB_Vector_free(&lhs_residual) ;
            return info ;
        }
        generic_lhs = lhs_residual ;
    }
    if (rhs_row_extended != GrB_NULL)
    {
        GrB_Info info = LAGraph_RPQMatrix_mnc_right_exact(&exact_nnz, &remaining_cols, &rhs_residual, generic_lhs, rhs_rows, rhs_cols, rhs_row_extended, rhs_n) ;
        if (info != GrB_SUCCESS)
        {
            GrB_Vector_free(&lhs_residual) ;
            GrB_Vector_free(&rhs_residual) ;
            return info ;
        }
        generic_rhs = rhs_residual ;
    }
    double p = (double) remaining_rows * (double) remaining_cols ;
    GrB_Info info = LAGraph_RPQMatrix_mnc_generic_nnz(&generic_nnz, generic_lhs, generic_rhs, p) ;
    GrB_Vector_free(&lhs_residual) ;
    GrB_Vector_free(&rhs_residual) ;
    if (info != GrB_SUCCESS)
    {
        return info ;
    }

    GrB_Index lhs_half_rows ;
    GrB_Index rhs_half_cols ;
    OK(LAGraph_RPQMatrix_count_vector_count_where(&lhs_half_rows, lhs_rows, (double) lhs_n / 2.0, 0)) ;
    OK(LAGraph_RPQMatrix_count_vector_count_where(&rhs_half_cols, rhs_cols, (double) lhs_n / 2.0, 0)) ;
    double lower_bound = (double) lhs_half_rows * (double) rhs_half_cols ;

    double support_cells = (double) lhs_nonempty_rows * (double) rhs_nonempty_cols ;
    *res = fmin(fmax(fmax(exact_nnz + generic_nnz, lower_bound), 0.0), fmin(support_cells, output_cells)) ;
    return (GrB_SUCCESS) ;
}

GrB_Info LAGraph_RPQMatrix_count_vector_sum(double *res, GrB_Vector vector)
{
    LG_ASSERT(res != NULL && vector != GrB_NULL, GrB_NULL_POINTER) ;
    return GrB_reduce(res, GrB_NULL, GrB_PLUS_MONOID_FP64, vector, GrB_NULL) ;
}

GrB_Info LAGraph_RPQMatrix_count_vector_scale(GrB_Vector *res, GrB_Vector vector, double scale, double cap)
{
    LG_ASSERT(res != NULL && vector != GrB_NULL, GrB_NULL_POINTER) ;
    GrB_Index size ;
    OK(GrB_Vector_size(&size, vector)) ;
    OK(GrB_Vector_new(res, GrB_FP64, size)) ;
    if (scale == 0.0)
    {
        return GrB_SUCCESS ;
    }
    GrB_Info info = GrB_apply(*res, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, scale, vector, GrB_NULL) ;
    if (info == GrB_SUCCESS)
    {
        info = GrB_apply(*res, GrB_NULL, GrB_NULL, GrB_MIN_FP64, cap, *res, GrB_NULL) ;
    }
    if (info != GrB_SUCCESS)
    {
        GrB_Vector_free(res) ;
    }
    return info ;
}

static GrB_Info LAGraph_RPQMatrix_count_vector_mnc_add_work(GrB_Vector *res, GrB_Vector *product, GrB_Vector lhs, GrB_Vector rhs, double lambda, double cap, GrB_Index size)
{
    OK(GrB_Vector_new(product, GrB_FP64, size)) ;
    OK(GrB_eWiseMult(*product, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, lhs, rhs, GrB_NULL)) ;
    OK(GrB_apply(*product, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, lambda, *product, GrB_NULL)) ;
    OK(GrB_Vector_new(res, GrB_FP64, size)) ;
    OK(GrB_eWiseAdd(*res, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, lhs, rhs, GrB_NULL)) ;
    OK(GrB_eWiseAdd(*res, GrB_NULL, GrB_NULL, GrB_MINUS_FP64, *res, *product, GrB_NULL)) ;
    OK(GrB_apply(*res, GrB_NULL, GrB_NULL, GrB_MAX_FP64, 0.0, *res, GrB_NULL)) ;
    OK(GrB_apply(*res, GrB_NULL, GrB_NULL, GrB_MIN_FP64, cap, *res, GrB_NULL)) ;
    OK(GrB_select(*res, GrB_NULL, GrB_NULL, GrB_VALUEGT_FP64, *res, 0.0, GrB_NULL)) ;
    return GrB_SUCCESS ;
}

GrB_Info LAGraph_RPQMatrix_count_vector_mnc_add(GrB_Vector *res, GrB_Vector lhs, GrB_Vector rhs, double lambda, double cap)
{
    LG_ASSERT(res != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(lhs != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT(rhs != GrB_NULL, GrB_NULL_POINTER) ;

    GrB_Index lhs_size ;
    GrB_Index rhs_size ;
    OK(GrB_Vector_size(&lhs_size, lhs)) ;
    OK(GrB_Vector_size(&rhs_size, rhs)) ;
    if (lhs_size != rhs_size)
    {
        return GrB_DIMENSION_MISMATCH ;
    }

    GrB_Vector product = GrB_NULL ;
    *res = GrB_NULL ;
    GrB_Info info = LAGraph_RPQMatrix_count_vector_mnc_add_work(res, &product, lhs, rhs, lambda, cap, lhs_size) ;
    GrB_Vector_free(&product) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Vector_free(res) ;
    }
    return info ;
}
