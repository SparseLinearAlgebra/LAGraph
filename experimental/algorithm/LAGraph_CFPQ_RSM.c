#include "LAGraph.h"
#include "LG_internal.h"
#include <GraphBLAS.h>

#include <stddef.h>
#include <stdint.h>
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
        GrB_free(&source_mask);                                                                    \
        if (graph_nt != NULL)                                                                      \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterm; ++_i)                                            \
                GrB_free(&graph_nt[_i]);                                                           \
            LAGraph_Free((void **)&graph_nt, msg);                                                 \
        }                                                                                          \
        if (rsm_start != NULL)                                                                     \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterm; ++_i)                                            \
                GrB_free(&rsm_start[_i]);                                                          \
            LAGraph_Free((void **)&rsm_start, msg);                                                \
        }                                                                                          \
        if (rsm_call != NULL)                                                                      \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterm; ++_i)                                            \
                GrB_free(&rsm_call[_i]);                                                           \
            LAGraph_Free((void **)&rsm_call, msg);                                                 \
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

static inline GrB_Info reshape( // input/output:
    GrB_Matrix *C,              // input/output matrix, reshaped in place
    bool by_col,                // true if reshape by column, false if by row
    GrB_Index nrows_new,        // new number of rows of C
    GrB_Index ncols_new,        // new number of columns of C
    const GrB_Descriptor desc   // to control # of threads used
)
{
    if (C == NULL || *C == NULL)
        return GrB_NULL_POINTER;

    GrB_Matrix temp = GrB_NULL;
    H_TRY(GxB_Matrix_reshapeDup(&temp, *C, by_col, nrows_new, ncols_new, desc));
    H_TRY(GrB_Matrix_free(C));

    *C = temp;
    return GrB_SUCCESS;
}

static GrB_Info s_mxm_chain(GrB_Matrix C,      // output accumulator  |Q|*|V*V|
                            GrB_Matrix A,      // RSM matrix          |Q|*|Q|
                            GrB_Matrix M,      // frontier            |Q|*|V*V|
                            GrB_Matrix B,      // graph matrix        |V|*|V|
                            GrB_Matrix *temp1, // scratch             |Q|*|V*V|
                            GrB_Matrix *temp2, // scratch             |Q|*|V*V|
                            GrB_Index Q, GrB_Index V)
{
    H_TRY(GrB_Matrix_clear(*temp1));
    H_TRY(GrB_Matrix_clear(*temp2));

    // temp1 = A^T * M -> |Q|*|V*V|
    H_TRY(GrB_mxm(*temp1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, A, M, GrB_DESC_T0));

    // reshape: |Q| × |V*V| ->  |Q*V|*|V|
    H_TRY(reshape(temp1, false, Q * V, V, GrB_NULL));
    H_TRY(reshape(temp2, false, Q * V, V, GrB_NULL));

    // temp2 = temp1 * B -> |Q*V|*|V|
    H_TRY(GrB_mxm(*temp2, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, *temp1, B, GrB_NULL));

    // reshape: |Q*V|*|V| -> |Q|*|V*V|
    H_TRY(reshape(temp1, false, Q, V * V, GrB_NULL));
    H_TRY(reshape(temp2, false, Q, V * V, GrB_NULL));

    // C |= temp2
    H_TRY(GrB_eWiseAdd(C, GrB_NULL, GrB_NULL, GrB_LOR, C, *temp2, GrB_NULL));

    return GrB_SUCCESS;
}

static GrB_Info build_call_matrix(GrB_Matrix *call, GrB_Matrix rsm_nt, GrB_Vector rsm_start,
                                  GrB_Index Q)
{
    GrB_Vector call_mask = GrB_NULL;

    H_TRY(GrB_Matrix_new(call, GrB_BOOL, Q, Q));
    H_TRY(GrB_Vector_new(&call_mask, GrB_BOOL, Q));

    // [q_call, q_ret] -> [q_call]
    H_TRY(GrB_reduce(call_mask, GrB_NULL, GrB_NULL, GrB_LOR_MONOID_BOOL, rsm_nt, GrB_NULL));
    // |Q|*1 * 1*|Q| -> |Q|*|Q|
    H_TRY(GrB_mxm(*call, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, (GrB_Matrix)call_mask,
                  (GrB_Matrix)rsm_start, GrB_DESC_T1));

    H_TRY(GrB_Vector_free(&call_mask));

    return GrB_SUCCESS;
}

#undef H_TRY
typedef struct
{
    GrB_Index state_count;            // Q
    GrB_Index terminal_count;         // num_terminals
    GrB_Index nonterminal_count;      // num_nonterminals
    GrB_Index start_nonterminal;      // start_nonterm
    GrB_Matrix *terminal_matrices;    // rsm_term
    GrB_Matrix *nonterminal_matrices; // rsm_nonterm
    GrB_Index *start_states;          // rsm_start
    GrB_Vector *final_states;         // rsm_final
} RSM;

int LAGraph_CFPQ_RSM(GrB_Vector *reachable, const RSM *rsm, const GrB_Matrix *graph_term,
                     const GrB_Index *sources, size_t num_sources, GrB_Index V, char *msg)
{
    LG_CLEAR_MSG;

    // shortcuts
    GrB_Index VV = V * V;
    GrB_Index Q = rsm->state_count;
    GrB_Index num_term = rsm->terminal_count;
    GrB_Index num_nonterm = rsm->nonterminal_count;
    GrB_Index start_nonterm = rsm->start_nonterminal;
    GrB_Matrix *rsm_term = rsm->terminal_matrices;
    GrB_Matrix *rsm_nonterm = rsm->nonterminal_matrices;
    GrB_Index *start_states = rsm->start_states;
    GrB_Vector *final_states = rsm->final_states;

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
    GrB_Matrix diag_entry = GrB_NULL;
    GrB_Matrix outer_result = GrB_NULL;
    GrB_Vector entry_vec = GrB_NULL;
    GrB_Vector source_mask = GrB_NULL;

    GrB_Matrix *graph_nt = NULL;
    GrB_Matrix *rsm_call = NULL;
    GrB_Vector *rsm_start = NULL; // rsm_start[S][i] = true, iff start_states[S] = i

    LG_ASSERT(reachable != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_term != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rsm_nonterm != NULL, GrB_NULL_POINTER);
    LG_ASSERT(start_states != NULL, GrB_NULL_POINTER);
    LG_ASSERT(final_states != NULL, GrB_NULL_POINTER);
    LG_ASSERT(start_nonterm < num_nonterm, GrB_INVALID_VALUE);
    LG_ASSERT(num_sources <= V, GrB_INVALID_VALUE);

    *reachable = GrB_NULL;

    // allocate per-nonterminal arrays
    LG_TRY(LAGraph_Malloc((void **)&graph_nt, rsm->nonterminal_count, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Malloc((void **)&rsm_start, rsm->nonterminal_count, sizeof(GrB_Vector), msg));
    LG_TRY(LAGraph_Malloc((void **)&rsm_call, rsm->nonterminal_count, sizeof(GrB_Matrix), msg));

    for (size_t i = 0; i < num_nonterm; ++i)
    {
        graph_nt[i] = GrB_NULL;
        rsm_start[i] = GrB_NULL;
        rsm_call[i] = GrB_NULL;
    }

    for (size_t i = 0; i < num_nonterm; ++i)
    {
        GRB_TRY(GrB_Matrix_new(&graph_nt[i], GrB_BOOL, V, V));
        GRB_TRY(GrB_Vector_new(&rsm_start[i], GrB_BOOL, Q));
        GrB_Vector_setElement_BOOL(rsm_start[i], true, start_states[i]);
        GRB_TRY(build_call_matrix(&rsm_call[i], rsm_nonterm[i], rsm_start[i], Q));
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
    // there is no need to create matrix for diag_entry
    GRB_TRY(GrB_Matrix_new(&outer_result, GrB_BOOL, Q, VV));
    GRB_TRY(GrB_Vector_new(&entry_vec, GrB_BOOL, V));
    GRB_TRY(GrB_Vector_new(&source_mask, GrB_BOOL, V));
    GRB_TRY(GrB_Vector_new(reachable, GrB_BOOL, V));

    // seed initial frontier M and visited P
    {

        GrB_Index nstart = 0;
        GRB_TRY(GrB_Vector_nvals(&nstart, rsm_start[start_nonterm]));

        GrB_Index *seed_states = NULL;
        bool *seed_vals = NULL;
        LG_TRY(LAGraph_Malloc((void **)&seed_states, nstart, sizeof(GrB_Index), msg));
        LG_TRY(LAGraph_Malloc((void **)&seed_vals, nstart, sizeof(bool), msg));

        GrB_Info info = GrB_Vector_extractTuples_BOOL(seed_states, seed_vals, &nstart,
                                                      rsm_start[start_nonterm]);
        LAGraph_Free((void **)&seed_vals, GrB_NULL);

        if (info != GrB_SUCCESS)
        {
            LAGraph_Free((void **)&seed_states, GrB_NULL);
            LG_FREE_ALL;
            return info;
        }

        for (size_t i = 0; i < num_sources; ++i)
        {
            GrB_Index seed_col = sources[i] * V + sources[i];
            for (GrB_Index k = 0; k < nstart; ++k)
            {
                GrB_Index q = seed_states[k];
                GRB_TRY(GrB_Matrix_setElement_BOOL(M, true, q, seed_col));
                GRB_TRY(GrB_Matrix_setElement_BOOL(P, true, q, seed_col));
            }
        }

        LAGraph_Free((void **)&seed_states, GrB_NULL);
    }

    // main fixed-point loop
    GrB_Index m_nvals = 0;
    GRB_TRY(GrB_Matrix_nvals(&m_nvals, M));
    while (m_nvals > 0)
    {
        // PHASE 1: TERMINAL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_term));
        for (size_t i = 0; i < num_term; ++i)
        {
            if (rsm_term[i] == GrB_NULL || graph_term[i] == GrB_NULL)
                continue;
            LG_TRY(s_mxm_chain(M_term, rsm_term[i], M, graph_term[i], &temp1, &temp2, Q, V));
        }

        // PHASE 2: NON-TERMINAL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_nonterm));
        for (size_t i = 0; i < num_nonterm; ++i)
        {
            if (rsm_nonterm[i] == GrB_NULL)
                continue;

            GrB_Index gnt_nvals = 0;
            GRB_TRY(GrB_Matrix_nvals(&gnt_nvals, graph_nt[i]));
            if (gnt_nvals == 0)
                continue;

            LG_TRY(s_mxm_chain(M_nonterm, rsm_nonterm[i], M, graph_nt[i], &temp1, &temp2, Q, V));
        }

        // PHASE 3: CALL TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_call));
        for (size_t i = 0; i < num_nonterm; ++i)
        {
            if (rsm_call[i] == GrB_NULL || rsm_nonterm[i] == GrB_NULL)
                continue;

            // mask_call = rsm_call[i]^T * M -> |Q|*|V*V|
            GRB_TRY(GrB_Matrix_clear(mask_call));
            GRB_TRY(GrB_mxm(mask_call, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, rsm_call[i],
                            M, GrB_DESC_T0));

            GrB_Index mc_nvals = 0;
            GRB_TRY(GrB_Matrix_nvals(&mc_nvals, mask_call));
            if (mc_nvals == 0)
                continue;

            // entry_vec = column-OR of mask_call reshaped to |Q*V|*|V|
            GRB_TRY(GrB_Matrix_clear(temp1));
            GRB_TRY(GrB_eWiseAdd(temp1, GrB_NULL, GrB_NULL, GrB_LOR, temp1, mask_call, GrB_NULL));
            GRB_TRY(reshape(&temp1, false, Q * V, V, GrB_NULL));

            GRB_TRY(GrB_Vector_clear(entry_vec));
            // column OR = row OR of transpose
            GRB_TRY(
                GrB_reduce(entry_vec, GrB_NULL, GrB_NULL, GrB_LOR_MONOID_BOOL, temp1, GrB_DESC_T0));

            GRB_TRY(reshape(&temp1, false, Q, VV, GrB_NULL));

            GrB_Index ev_nvals = 0;
            GRB_TRY(GrB_Vector_nvals(&ev_nvals, entry_vec));
            if (ev_nvals == 0)
                continue;

            // Build diag(entry_vec) and flatten to 1*|V*V|
            GRB_TRY(GrB_Matrix_free(&diag_entry));
            diag_entry = GrB_NULL;
            GRB_TRY(GrB_Matrix_diag(&diag_entry, entry_vec, 0)); // |V|*|V|
            GRB_TRY(reshape(&diag_entry, false, 1, VV, GrB_NULL));

            GRB_TRY(GrB_Matrix_clear(outer_result));
            // |Q|*1 * 1*|V*V| -> |Q|*|V*V|
            GRB_TRY(GrB_mxm(outer_result, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL,
                            (GrB_Matrix)rsm_start[i], diag_entry, GrB_NULL));

            // M_call |= outer_result
            GRB_TRY(
                GrB_eWiseAdd(M_call, GrB_NULL, GrB_NULL, GrB_LOR, M_call, outer_result, GrB_NULL));
        }

        // PHASE 4: RETURN TRANSITIONS
        GRB_TRY(GrB_Matrix_clear(M_return));
        for (size_t i = 0; i < num_nonterm; ++i)
        {
            if (final_states[i] == GrB_NULL || rsm_nonterm[i] == GrB_NULL)
                continue;

            // 1*|Q| * |Q|*|V*V| -> 1*|V*V|
            GRB_TRY(GrB_Matrix_clear(new_edges));
            GRB_TRY(GrB_mxm(new_edges, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                            (GrB_Matrix)final_states[i], M, GrB_DESC_T0));

            GrB_Index ne_vals = 0;
            GRB_TRY(GrB_Matrix_nvals(&ne_vals, new_edges));
            if (ne_vals == 0)
                continue;

            GRB_TRY(reshape(&new_edges, false, V, V, GrB_NULL));
            GRB_TRY(GrB_eWiseAdd(graph_nt[i], GrB_NULL, GrB_NULL, GrB_LOR, graph_nt[i], new_edges,
                                 GrB_NULL));

            LG_TRY(s_mxm_chain(M_return, rsm_nonterm[i], P, new_edges, &temp1, &temp2, Q, V));

            GRB_TRY(reshape(&new_edges, false, 1, VV, GrB_NULL));
        }

        // PHASE 5: MERGE & MARK VISITED
        // GrB_assign has big overhead
        // M_new = (M_term | M_nonterm | M_call | M_return) & ~P
        GRB_TRY(GrB_Matrix_clear(M_new));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_NULL, GrB_LOR, M_new, M_term, GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_NULL, GrB_LOR, M_new, M_nonterm, GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_NULL, GrB_LOR, M_new, M_call, GrB_NULL));
        GRB_TRY(GrB_eWiseAdd(M_new, GrB_NULL, GrB_NULL, GrB_LOR, M_new, M_return, GrB_NULL));

        // M_new &= ~P
        GRB_TRY(GrB_assign(M_new, P, GrB_NULL, M_new, GrB_ALL, Q, GrB_ALL, VV, GrB_DESC_RSC));

        // P |= M_new
        GRB_TRY(GrB_eWiseAdd(P, GrB_NULL, GrB_NULL, GrB_LOR, P, M_new, GrB_NULL));

        GrB_Matrix swap = M;
        M = M_new;
        M_new = swap;

        GRB_TRY(GrB_Matrix_nvals(&m_nvals, M));
    }

    // Extract all reachable vertices
    {
        bool *vals = NULL;
        LG_TRY(LAGraph_Malloc((void **)&vals, num_sources, sizeof(bool), msg));
        for (size_t i = 0; i < num_sources; ++i)
            vals[i] = true;

        GrB_Info _bi =
            GrB_Vector_build_BOOL(source_mask, sources, vals, num_sources, GrB_SECOND_BOOL);
        LAGraph_Free((void **)&vals, GrB_NULL);
        GRB_TRY(_bi);
    }

    GRB_TRY(GrB_mxv(*reachable, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                    graph_nt[start_nonterm], source_mask, GrB_DESC_T0));

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
