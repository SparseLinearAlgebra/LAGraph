#include "LAGraph.h"
#include "LG_internal.h"
#include <GraphBLAS.h>

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define LG_FREE_WORK                                                                               \
    {                                                                                              \
        GrB_free(&M);                                                                              \
        GrB_free(&P);                                                                              \
        GrB_free(&M_term);                                                                         \
        GrB_free(&M_nonterm);                                                                      \
        GrB_free(&M_call);                                                                         \
        GrB_free(&M_return);                                                                       \
        GrB_free(&M_new);                                                                          \
        GrB_free(&temp1);                                                                          \
        GrB_free(&temp2);                                                                          \
        GrB_free(&mask_call);                                                                      \
        GrB_free(&frame);                                                                          \
        GrB_free(&new_edges);                                                                      \
        GrB_free(&diag_entry);                                                                     \
        GrB_free(&outer_result);                                                                   \
        GrB_free(&entry_vec);                                                                      \
        if (Stack != NULL)                                                                         \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterminals; ++_i)                                       \
                GrB_free(&Stack[_i]);                                                              \
            LAGraph_Free((void **)&Stack, GrB_NULL);                                               \
        }                                                                                          \
        if (graph_nt != NULL)                                                                      \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterminals; ++_i)                                       \
                GrB_free(&graph_nt[_i]);                                                           \
            LAGraph_Free((void **)&graph_nt, GrB_NULL);                                            \
        }                                                                                          \
    }
#define LG_FREE_ALL                                                                                \
    {                                                                                              \
        LG_FREE_WORK;                                                                              \
        GrB_free(reachable);                                                                       \
    }

#define H_TRY(expr)                                                                                \
    {                                                                                              \
        GrB_Info _hi = (expr);                                                                     \
        if (_hi != GrB_SUCCESS)                                                                    \
            return _hi;                                                                            \
    }

static GrB_Info s_mxm_chain(GrB_Matrix C,     // output accumulator  |Q|*|V*V|
                            GrB_Matrix A,     // RSM matrix          |Q|*|Q|
                            GrB_Matrix M,     // frontier            |Q|*|V*V|
                            GrB_Matrix B,     // graph matrix        |V|*|V|
                            GrB_Matrix temp1, // scratch             |Q|*|V*V|
                            GrB_Matrix temp2, // scratch             |Q|*|V*V|
                            GrB_Index Q, GrB_Index V)
{
    H_TRY(GrB_Matrix_clear(temp1));
    H_TRY(GrB_Matrix_clear(temp2));

    // temp1 = A^T * M -> |Q|*|V*V|
    H_TRY(GrB_mxm(temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, A, M, GrB_DESC_T0));

    // reshape: |Q| × |V*V| ->  |Q*V|*|V|
    H_TRY(GxB_Matrix_reshape(temp1, false, Q * V, V, GrB_NULL));
    H_TRY(GxB_Matrix_reshape(temp2, false, Q * V, V, GrB_NULL));

    // temp2 = temp1 * B -> |Q*V|*|V|
    H_TRY(GrB_mxm(temp2, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, temp1, B, GrB_NULL));

    // reshape: |Q*V|*|V| -> |Q|*|V*V|
    H_TRY(GxB_Matrix_reshape(temp1, false, Q, V * V, GrB_NULL));
    H_TRY(GxB_Matrix_reshape(temp2, false, Q, V * V, GrB_NULL));

    // C |= temp2
    H_TRY(GrB_eWiseAdd(C, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, C, temp2, GrB_NULL));

    return GrB_SUCCESS;
}
#undef H_TRY

int LAGraph_CFPQ_RSM(
    // output:
    GrB_Vector *reachable,

    // RSM terminal input:
    size_t num_terminals, GrB_Matrix *rsm_term, GrB_Matrix *graph_term,

    // RSM non-terminal input:
    size_t num_nonterminals, GrB_Matrix *rsm_nonterm, GrB_Matrix *rsm_call, GrB_Vector *rsm_start,
    GrB_Vector *rsm_final,

    size_t start_nonterm, GrB_Index start_vertex,

    GrB_Index Q, GrB_Index V,

    char *msg)
{
    LG_CLEAR_MSG;

    // working matrices / vectors
    GrB_Matrix M = GrB_NULL;
    GrB_Matrix P = GrB_NULL;
    GrB_Matrix M_term = GrB_NULL;
    GrB_Matrix M_nonterm = GrB_NULL;
    GrB_Matrix M_call = GrB_NULL;
    GrB_Matrix M_return = GrB_NULL;
    GrB_Matrix M_new = GrB_NULL;
    GrB_Matrix temp1 = GrB_NULL;
    GrB_Matrix temp2 = GrB_NULL;
    GrB_Matrix mask_call = GrB_NULL;
    GrB_Matrix frame = GrB_NULL;
    GrB_Matrix new_edges = GrB_NULL;
    GrB_Matrix diag_entry = GrB_NULL;   // add diagonal matrix for entry_vec
    GrB_Matrix outer_result = GrB_NULL; // additional matrix
    GrB_Vector entry_vec = GrB_NULL;

    GrB_Matrix *Stack = NULL;
    GrB_Matrix *graph_nt = NULL;

    LG_ASSERT(reachable != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_term != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_nonterm != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_call != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_start != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_final != NULL, GrB_NULL_POINTER);
    LG_ASSERT(start_nonterm < num_nonterminals, GrB_INVALID_VALUE);
    LG_ASSERT(start_vertex < V, GrB_INVALID_VALUE);

    *reachable = GrB_NULL;

    GrB_Index VV = V * V;

    // allocate per-nonterminal arrays
    LG_TRY(LAGraph_Malloc((void **)&Stack, num_nonterminals, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Malloc((void **)&graph_nt, num_nonterminals, sizeof(GrB_Matrix), msg));

    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        Stack[i] = GrB_NULL;
        graph_nt[i] = GrB_NULL;
    }

    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        GRB_TRY(GrB_Matrix_new(&Stack[i], GrB_BOOL, Q, VV));
        GRB_TRY(GrB_Matrix_new(&graph_nt[i], GrB_BOOL, V, V));
    }

    // allocate working matrices
    GRB_TRY(GrB_Matrix_new(&M, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&P, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&M_term, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&M_nonterm, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&M_call, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&M_return, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&M_new, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&temp1, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&temp2, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&mask_call, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&frame, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Matrix_new(&new_edges, GrB_BOOL, 1, VV));
    // GRB_TRY(GrB_Matrix_new(&diag_entry, GrB_BOOL, V, V));
    GRB_TRY(GrB_Matrix_new(&outer_result, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Vector_new(&entry_vec, GrB_BOOL, V));
    GRB_TRY(GrB_Vector_new(reachable, GrB_BOOL, V));

    // seed initial frontier M and Stack
    {
        GrB_Index seed_col = start_vertex * V + start_vertex;

        GrB_Index nstart = 0;
        GRB_TRY(GrB_Vector_nvals(&nstart, rsm_start[start_nonterm]));

        GrB_Index *start_states = NULL;
        bool *start_vals = NULL;
        LG_TRY(LAGraph_Malloc((void **)&start_states, nstart, sizeof(GrB_Index), msg));
        LG_TRY(LAGraph_Malloc((void **)&start_vals, nstart, sizeof(bool), msg));

        GrB_Info info = GrB_Vector_extractTuples_BOOL(start_states, start_vals, &nstart,
                                                      rsm_start[start_nonterm]);
        LAGraph_Free((void **)&start_vals, GrB_NULL);

        if (info != GrB_SUCCESS)
        {
            LAGraph_Free((void **)&start_states, GrB_NULL);
            LG_FREE_ALL;
            return info;
        }

        for (GrB_Index k = 0; k < nstart; ++k)
        {
            GrB_Index q = start_states[k];
            GRB_TRY(GrB_Matrix_setElement_BOOL(M, true, q, seed_col));
            GRB_TRY(GrB_Matrix_setElement_BOOL(Stack[start_nonterm], true, q, seed_col));
        }

        LAGraph_Free((void **)&start_states, GrB_NULL);
    }

    // main fixed-point loop
    GrB_Index m_nvals = 0;
    GRB_TRY(GrB_Matrix_nvals(&m_nvals, M));

    while (m_nvals > 0)
    {

        // PHASE 1: TERMINAL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_term));
        for (size_t i = 0; i < num_terminals; ++i)
        {
            if (rsm_term[i] == GrB_NULL || graph_term[i] == GrB_NULL)
                continue;
            LG_TRY(s_mxm_chain(M_term, rsm_term[i], M, graph_term[i], temp1, temp2, Q, V));
        }

        // PHASE 2: NON-TERMINAL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_nonterm));
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (rsm_nonterm[i] == GrB_NULL)
                continue;
            LG_TRY(s_mxm_chain(M_nonterm, rsm_nonterm[i], M, graph_nt[i], temp1, temp2, Q, V));
        }

        // PHASE 3: CALL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_call));
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (rsm_call[i] == GrB_NULL || rsm_nonterm[i] == GrB_NULL)
                continue;

            // mask_call = rsm_call[i]^T * M -> |Q|*|V*V|
            GRB_TRY(GrB_Matrix_clear(mask_call));
            GRB_TRY(GrB_mxm(mask_call, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, rsm_call[i],
                            M, GrB_DESC_T0));

            GrB_Index mc_nvals = 0;
            GRB_TRY(GrB_Matrix_nvals(&mc_nvals, mask_call));
            if (mc_nvals == 0)
                continue;

            // frame = rsm_nonterm[i]^T * M, accumulated into Stack[i]
            GRB_TRY(GrB_Matrix_clear(frame));
            GRB_TRY(GrB_mxm(frame, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, rsm_nonterm[i], M,
                            GrB_DESC_T0));
            GRB_TRY(GrB_eWiseAdd(Stack[i], GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, Stack[i],
                                 frame, GrB_NULL));

            // entry_vec = column-OR of mask_call reshaped to |Q*V|*|V|
            GRB_TRY(GrB_Matrix_clear(temp1));
            GRB_TRY(GrB_eWiseAdd(temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, temp1,
                                 mask_call, GrB_NULL));
            GRB_TRY(GxB_Matrix_reshape(temp1, false, Q * V, V, GrB_NULL));

            GRB_TRY(GrB_Vector_clear(entry_vec));
            // column OR = row OR of transpose
            GRB_TRY(
                GrB_reduce(entry_vec, GrB_NULL, GrB_LOR, GrB_LOR_MONOID_BOOL, temp1, GrB_DESC_T0));

            GRB_TRY(GxB_Matrix_reshape(temp1, false, Q, VV, GrB_NULL));

            // Build diag(entry_vec) reshaped to 1*|V*V|
            diag_entry = GrB_NULL;
            GRB_TRY(GrB_Matrix_diag(&diag_entry, entry_vec, 0)); // |V|*|V|
            GRB_TRY(GxB_Matrix_reshape(diag_entry, false, 1, VV, GrB_NULL));

            GRB_TRY(GrB_Matrix_clear(outer_result));
            // |Q|*1 * 1*|V*V| -> |Q|*|V*V|
            GRB_TRY(GrB_mxm(outer_result, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL,
                            (GrB_Matrix)rsm_start[i], diag_entry, GrB_NULL));

            // M_call |= outer_result
            GRB_TRY(GrB_eWiseAdd(M_call, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_call,
                                 outer_result, GrB_NULL));

            GRB_TRY(GrB_Matrix_free(&diag_entry));
        }

        // PHASE 4: RETURN TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_return));
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (rsm_final[i] == GrB_NULL)
                continue;

            // 1*|Q| * |Q|*|V*V| -> 1*|V*V|
            GRB_TRY(GrB_Matrix_clear(new_edges));
            GRB_TRY(GrB_mxm(new_edges, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                            (GrB_Matrix)rsm_final[i], M, GrB_DESC_T0));

            GRB_TRY(GxB_Matrix_reshape(new_edges, false, V, V, GrB_NULL));

            GRB_TRY(GrB_eWiseAdd(graph_nt[i], GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 graph_nt[i], new_edges, GrB_NULL));

            GRB_TRY(GrB_Matrix_clear(temp1));
            GRB_TRY(GrB_eWiseAdd(temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, temp1,
                                 Stack[i], GrB_NULL));
            GRB_TRY(GxB_Matrix_reshape(temp1, false, Q * V, V, GrB_NULL));

            GRB_TRY(GrB_Matrix_clear(temp2));
            GRB_TRY(GxB_Matrix_reshape(temp2, false, Q * V, V, GrB_NULL));
            // |Q*V|*|V| * |V|*|V| -> |Q*V|*|V|
            GRB_TRY(GrB_mxm(temp2, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, temp1, new_edges,
                            GrB_NULL));

            GRB_TRY(GxB_Matrix_reshape(temp1, false, Q, VV, GrB_NULL));
            GRB_TRY(GxB_Matrix_reshape(temp2, false, Q, VV, GrB_NULL));

            GRB_TRY(GrB_eWiseAdd(M_return, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_return,
                                 temp2, GrB_NULL));
            GRB_TRY(GxB_Matrix_reshape(new_edges, false, 1, VV, GrB_NULL));
        }

        // PHASE 5: MERGE & MARK VISITED
        // GrB_assign has big overhead
        // M_new = (M_term | M_nonterm | M_call | M_return) & ~P
        GRB_TRY(GrB_Matrix_clear(M_new));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_new, M_term,
                             GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_new, M_nonterm,
                             GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_new, M_call,
                             GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, M_new, M_return,
                             GrB_NULL));

        // M_new &= ~P
        GRB_TRY(GrB_assign(M_new, P, GrB_NULL, M_new, GrB_ALL, Q, GrB_ALL, VV, GrB_DESC_RSC));

        // P |= M_new
        GRB_TRY(GrB_eWiseAdd(P, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, P, M_new, GrB_NULL));

        GrB_Matrix swap = M;
        M = M_new;
        M_new = swap;

        GRB_TRY(GrB_Matrix_nvals(&m_nvals, M));
    }

    // GrB_DESC_T0 transposes graph_nt so Col_extract gives row start_vertex
    GRB_TRY(GrB_Col_extract(*reachable, GrB_NULL, GrB_NULL, graph_nt[start_nonterm], GrB_ALL, V,
                            start_vertex, GrB_DESC_T0));

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
