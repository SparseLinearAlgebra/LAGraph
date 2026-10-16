#define LG_FREE_WORK                                                         \
{                                                                            \
    GrB_free (&v_diag);                                                      \
    GrB_free (&identity_matrix);                                             \
    GrB_free (&start_matrix);                                                \
    GrB_free (&final_column);                                                \
    GrB_free (&Rblock);                                                      \
    GrB_free (&delta);                                                       \
    GrB_free (&DeltaM);                                                      \
    GrB_free (&M_work);                                                      \
    GrB_free (&Reach);                                                       \
    GrB_free (&Frontier);                                                    \
    GrB_free (&NextFrontier);                                                \
    if (Sel != NULL)                                                         \
    {                                                                        \
        for (int64_t k = 0; k < nt_count; ++k)                               \
        {                                                                    \
            GrB_free (&Sel [k]);                                             \
        }                                                                    \
    }                                                                        \
    LAGraph_Free ((void **) &Sel, NULL);                                     \
}

#define LG_FREE_ALL                                                          \
{                                                                            \
    LG_FREE_WORK;                                                            \
    for (int64_t k = 0; outputs != NULL && k < nt_count; ++k)                \
    {                                                                        \
        GrB_free (&outputs [k]);                                             \
    }                                                                        \
    if (M_intersect_out != NULL) GrB_free (M_intersect_out);                 \
}

#include "LG_internal.h"
#include <LAGraphX.h>

GrB_Info LAGraph_CFL_AllPaths_Kronecker
(
    GrB_Matrix *M_intersect_out,                                        
    GrB_Matrix *outputs,                                                                                  
    const GrB_Matrix *adj_matrices,                                    
    int64_t terms_count,            
    RSM *rsm,                       
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    GrB_Vector v_diag = NULL;
    GrB_Matrix identity_matrix = NULL, start_matrix = NULL, final_column = NULL;
    GrB_Matrix M_work = NULL, DeltaM = NULL;
    GrB_Matrix Reach = NULL;
    GrB_Matrix Frontier = NULL, NextFrontier = NULL;
    GrB_Matrix Rblock = NULL, delta = NULL;
    GrB_Matrix *Sel = NULL;
    int64_t nt_count = 0;
    LG_CLEAR_MSG;

    LG_ASSERT_MSG (M_intersect_out != NULL, GrB_NULL_POINTER,
        "The output matrix pointer cannot be null.");
    *M_intersect_out = NULL;
    LG_ASSERT_MSG (outputs != NULL, GrB_NULL_POINTER,
        "The outputs array cannot be null.");
    LG_ASSERT_MSG (adj_matrices != NULL, GrB_NULL_POINTER,
        "The adjacency matrices array cannot be null.");
    LG_ASSERT_MSG (rsm != NULL, GrB_NULL_POINTER, "The RSM cannot be null.");
    LG_ASSERT_MSG (terms_count > 0, GrB_INVALID_VALUE,
        "The number of terminals must be greater than zero.");
    LG_ASSERT_MSG (rsm->nonterminal_count > 0, GrB_INVALID_VALUE,
        "The number of non-terminals must be greater than zero.");
    for (int64_t t = 0; t < terms_count; ++t)
    {
        LG_ASSERT_MSGF (adj_matrices [t] != NULL, GrB_NULL_POINTER,
            "Adjacency matrix with index %" PRId64 " is null.", t);
    }

    nt_count = rsm->nonterminal_count;
    for (int64_t nt = 0; nt < nt_count; ++nt)
    {
        outputs [nt] = NULL;
    }

    GrB_Index g_dim;
    GRB_TRY (GrB_Matrix_ncols (&g_dim, adj_matrices [0]));
    GrB_Index r_dim = rsm->state_count;
    GrB_Index kronecker_dim = r_dim * g_dim;
    GrB_Index sources_count = nt_count * g_dim;

    GRB_TRY (GrB_Vector_new (&v_diag, GrB_BOOL, g_dim));
    GRB_TRY (GrB_Vector_assign_BOOL (v_diag, NULL, NULL, true, GrB_ALL, g_dim,
        NULL));
    GRB_TRY (GrB_Matrix_diag (&identity_matrix, v_diag, 0));
    GRB_TRY (GrB_free (&v_diag));

    GRB_TRY (GrB_Matrix_new (&M_work, GrB_BOOL, kronecker_dim, kronecker_dim));
    GRB_TRY (GrB_Matrix_new (&DeltaM, GrB_BOOL, kronecker_dim, kronecker_dim));
    GRB_TRY (GrB_Matrix_new (&Frontier, GrB_BOOL, sources_count,
        kronecker_dim));
    GRB_TRY (GrB_Matrix_new (&NextFrontier, GrB_BOOL, sources_count,
        kronecker_dim));
    GRB_TRY (GrB_Matrix_new (&Rblock, GrB_BOOL, g_dim, kronecker_dim));
    GRB_TRY (GrB_Matrix_new (&delta, GrB_BOOL, g_dim, g_dim));
    GRB_TRY (GrB_Matrix_new (&start_matrix, GrB_BOOL, nt_count, r_dim));
    GRB_TRY (GrB_Matrix_new (&final_column, GrB_BOOL, r_dim, 1));
    LG_TRY (LAGraph_Calloc ((void **) &Sel, nt_count, sizeof (GrB_Matrix), msg));

    for (int64_t t = 0; t < terms_count; ++t)
    {
        GRB_TRY (GrB_kronecker (M_work, NULL, GxB_ANY_BOOL, GxB_PAIR_BOOL,
            rsm->terminal_matrices [t], adj_matrices [t], NULL));
    }

    for (int64_t nt = 0; nt < nt_count; ++nt)
    {
        GrB_Index s = rsm->start_states [nt];
        GRB_TRY (GrB_Matrix_setElement_BOOL (start_matrix, true, nt, s));

        bool is_final = false;
        GRB_TRY (GrB_Vector_extractElement_BOOL (&is_final,
            rsm->final_states [nt], s));

        if (is_final)
        {
            GRB_TRY (GrB_Matrix_dup (&outputs [nt], identity_matrix));
            GRB_TRY (GrB_kronecker (M_work, NULL, GxB_ANY_BOOL, GxB_PAIR_BOOL,
                rsm->nonterminal_matrices [nt], identity_matrix, NULL));
        }
        else
        {
            GRB_TRY (GrB_Matrix_new (&outputs [nt], GrB_BOOL, g_dim, g_dim));
        }

        GRB_TRY (GrB_Col_assign (final_column, NULL, NULL,
            rsm->final_states [nt], GrB_ALL, r_dim, 0, NULL));
        GRB_TRY (GrB_Matrix_new (&Sel [nt], GrB_BOOL, kronecker_dim, g_dim));
        GRB_TRY (GrB_kronecker (Sel [nt], NULL, NULL, GxB_PAIR_BOOL,
            final_column, identity_matrix, NULL));
    }

    GRB_TRY (GrB_kronecker (Frontier, NULL, NULL, GxB_PAIR_BOOL, start_matrix,
        identity_matrix, NULL));
    GRB_TRY (GrB_Matrix_dup (&Reach, Frontier));


    GrB_Index frontier_nvals;
    GRB_TRY (GrB_Matrix_nvals (&frontier_nvals, Frontier));
    while (frontier_nvals > 0)
    {
        GRB_TRY (GrB_Matrix_clear (DeltaM));
        for (int64_t nt = 0; nt < nt_count; ++nt)
        {
            GrB_Index rows [2] = { nt * g_dim, (nt + 1) * g_dim - 1 };
            GRB_TRY (GrB_Matrix_extract (Rblock, NULL, NULL, Frontier, rows,
                GxB_RANGE, GrB_ALL, kronecker_dim, NULL));

            GRB_TRY (GrB_mxm (delta, outputs [nt], NULL, GxB_ANY_PAIR_BOOL,
                Rblock, Sel [nt], GrB_DESC_RSC));

            GrB_Index delta_nvals;
            GRB_TRY (GrB_Matrix_nvals (&delta_nvals, delta));
            if (delta_nvals > 0)
            {
                GRB_TRY (GrB_eWiseAdd (outputs [nt], NULL, NULL, GxB_ANY_BOOL,
                    outputs [nt], delta, NULL));
                GRB_TRY (GrB_kronecker (DeltaM, NULL, GxB_ANY_BOOL,
                    GxB_PAIR_BOOL, rsm->nonterminal_matrices [nt], delta,
                    NULL));
            }
        }

        GRB_TRY (GrB_mxm (NextFrontier, Reach, NULL, GxB_ANY_PAIR_BOOL,
            Frontier, M_work, GrB_DESC_RSC));

        GrB_Index new_edges;
        GRB_TRY (GrB_Matrix_nvals (&new_edges, DeltaM));
        if (new_edges > 0)
        {
            GRB_TRY (GrB_mxm (NextFrontier, Reach, GxB_ANY_BOOL,
                GxB_ANY_PAIR_BOOL, Reach, DeltaM, GrB_DESC_SC));
            GRB_TRY (GrB_eWiseAdd (M_work, NULL, NULL, GxB_ANY_BOOL, M_work,
                DeltaM, NULL));
        }

        GRB_TRY (GrB_eWiseAdd (Reach, NULL, NULL, GxB_ANY_BOOL, Reach,
            NextFrontier, NULL));

        GrB_Matrix temp = Frontier;
        Frontier = NextFrontier;
        NextFrontier = temp;
        GRB_TRY (GrB_Matrix_nvals (&frontier_nvals, Frontier));
    }

    *M_intersect_out = M_work;
    M_work = NULL;
    LG_FREE_WORK;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
