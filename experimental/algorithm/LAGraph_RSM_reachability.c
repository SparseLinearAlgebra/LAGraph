//------------------------------------------------------------------------------
// LAGraph_RSM_reachability.c
//------------------------------------------------------------------------------

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "LAGraph.h"
#include "LG_internal.h"

#define LG_FREE_WORK                                                                               \
    {                                                                                              \
        GrB_free(&frontier);                                                                       \
        GrB_free(&next_frontier);                                                                  \
        GrB_free(&front_temp1);                                                                    \
        GrB_free(&front_temp2);                                                                    \
        GrB_free(&P);                                                                              \
        GrB_free(&K);                                                                              \
        GrB_free(&final_reducer);                                                                  \
        GrB_free(&M_ret);                                                                          \
        GrB_free(&M_call);                                                                         \
        GrB_free(&M_call_new);                                                                     \
        GrB_free(&G_S_new);                                                                        \
        GrB_free(&step1);                                                                          \
        LAGraph_Free((void **)&A_a, GrB_NULL);                                                     \
        LAGraph_Free((void **)&B_a, GrB_NULL);                                                     \
        LAGraph_Free((void **)&B_S, GrB_NULL);                                                     \
        LAGraph_Free((void **)&B_call, GrB_NULL);                                                  \
        LAGraph_Free((void **)&B_ret, GrB_NULL);                                                   \
        LAGraph_Free((void **)&B_S, GrB_NULL);                                                     \
        if (A_S != GrB_NULL)                                                                       \
        {                                                                                          \
            for (size_t _i = 0; _i < num_nonterminals; ++_i)                                       \
                GrB_free(&A_S[_i]);                                                                \
            LAGraph_Free((void **)&A_S, GrB_NULL);                                                 \
        }                                                                                          \
    }

#define LG_FREE_ALL                                                                                \
    {                                                                                              \
        LG_FREE_WORK;                                                                              \
        GrB_free(reachable);                                                                       \
    }

int LAGraph_RSM_reachability(
    // output:
    GrB_Vector *reachable,

    // input:
    size_t num_terminals, // # terminal labels
    LAGraph_Graph *G_a,   // terminal transitions
    LAGraph_Graph *N_a,   // terminal transitions

    size_t num_nonterminals, // # nonterminal labels
    LAGraph_Graph *N_S,      // nonterminal RSM transitions
    LAGraph_Graph *Call_S,   // call transition matrix
    LAGraph_Graph *Ret_S,    // return transition matrix

    const GrB_Index *QS, // staring states in RSM
    size_t nqs,          // # starting states
    const GrB_Index *QF, // final states in RSM
    size_t nqf,          // # final states

    const GrB_Index *Source, // sources vertices
    size_t ns,               // # source vertices

    char *msg

)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    // Working matrices
    GrB_Matrix frontier = GrB_NULL; // traversal frontier
    GrB_Matrix next_frontier = GrB_NULL;
    GrB_Matrix front_temp1 = GrB_NULL;
    GrB_Matrix front_temp2 = GrB_NULL;
    GrB_Matrix P = GrB_NULL;             // visited pairs (state, vertex)
    GrB_Matrix K = GrB_NULL;             // visited pairs in second loop
    GrB_Vector final_reducer = GrB_NULL; // vector for reducing the visited matrix

    // Backward-phase working matrices
    GrB_Matrix M_ret = GrB_NULL;
    GrB_Matrix M_call = GrB_NULL;
    GrB_Matrix M_call_new = GrB_NULL;
    GrB_Matrix G_S_new = GrB_NULL;
    GrB_Matrix step1 = GrB_NULL;

    GrB_Index ng = 0;      // # nodes in the graph
    GrB_Index nr = 0;      // # states in the RSM
    GrB_Index states = ns; // # pairs in current frontier

    GrB_Index rows = 0, cols = 0, vals = 0;

    GrB_Matrix *A_a = GrB_NULL; // for G_a
    GrB_Matrix *B_a = GrB_NULL; // for N_a

    GrB_Matrix *A_S = GrB_NULL;    // Derived graph edge matrices
    GrB_Matrix *B_S = GrB_NULL;    // for N_S
    GrB_Matrix *B_call = GrB_NULL; // for Call_S
    GrB_Matrix *B_ret = GrB_NULL;  // for Ret_S

    LG_ASSERT(reachable != GrB_NULL, GrB_NULL_POINTER);

    LG_ASSERT(G_a != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(N_a != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(N_S != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(Call_S != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(Ret_S != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(QS != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(QF != GrB_NULL, GrB_NULL_POINTER);
    LG_ASSERT(Source != GrB_NULL, GrB_NULL_POINTER);

    (*reachable) = GrB_NULL;

    //--------------------------------------------------------------------------
    // validate input graphs
    //--------------------------------------------------------------------------

    for (size_t i = 0; i < num_terminals; ++i)
        if (G_a[i] != GrB_NULL)
            LG_TRY(LAGraph_CheckGraph(G_a[i], msg))

    for (size_t i = 0; i < num_terminals; ++i)
        if (N_a[i] != GrB_NULL)
            LG_TRY(LAGraph_CheckGraph(N_a[i], msg))

    for (size_t i = 0; i < num_nonterminals; ++i)
        if (N_S[i] != GrB_NULL)
            LG_TRY(LAGraph_CheckGraph(N_S[i], msg))

    for (size_t i = 0; i < num_nonterminals; ++i)
        if (Call_S[i] != GrB_NULL)
            LG_TRY(LAGraph_CheckGraph(Call_S[i], msg))

    for (size_t i = 0; i < num_nonterminals; ++i)
        if (Ret_S[i] != GrB_NULL)
            LG_TRY(LAGraph_CheckGraph(Ret_S[i], msg))

    LG_TRY(LAGraph_Malloc((void **)&A_a, num_terminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (G_a[i] == GrB_NULL)
            A_a[i] = GrB_NULL;
        else
            A_a[i] = G_a[i]->A;
    }
    LG_TRY(LAGraph_Malloc((void **)&B_a, num_terminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (N_a[i] == GrB_NULL)
            B_a[i] = GrB_NULL;
        else
            B_a[i] = N_a[i]->A;
    }
    LG_TRY(LAGraph_Malloc((void **)&B_S, num_nonterminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (N_S[i] == GrB_NULL)
            B_S[i] = GrB_NULL;
        else
            B_S[i] = N_S[i]->A;
    }
    LG_TRY(LAGraph_Malloc((void **)&B_call, num_nonterminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (Call_S[i] == GrB_NULL)
            B_call[i] = GrB_NULL;
        else
            B_call[i] = Call_S[i]->A;
    }
    LG_TRY(LAGraph_Malloc((void **)&B_ret, num_nonterminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (Ret_S[i] == GrB_NULL)
            B_ret[i] = GrB_NULL;
        else
            B_ret[i] = Ret_S[i]->A;
    }

    //--------------------------------------------------------------------------
    // determine ng (graph size) and nr (RSM state count)
    //--------------------------------------------------------------------------

    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (A_a[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&ng, A_a[i]));
        break;
    }

    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (B_a[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&nr, B_a[i]));
        break;
    }

    //--------------------------------------------------------------------------
    // dimension checks
    //--------------------------------------------------------------------------

    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (A_a[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&rows, A_a[i]));
        GRB_TRY(GrB_Matrix_ncols(&cols, A_a[i]));

        LG_ASSERT_MSG(rows == ng && cols == ng, LAGRAPH_NOT_CACHED,
                      "all the matrices in the graph adjacency matrix decomposition "
                      "should have the same dimensions and be square");
    }
    for (size_t i = 0; i < num_terminals; ++i)
    {
        if (B_a[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&rows, B_a[i]));
        GRB_TRY(GrB_Matrix_ncols(&cols, B_a[i]));

        LG_ASSERT_MSG(rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
                      "all the matrices in the RSM adjacency matrix decomposition "
                      "should have the same dimensions and be square")
    }
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (B_S[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&rows, B_S[i]));
        GRB_TRY(GrB_Matrix_ncols(&cols, B_S[i]));

        LG_ASSERT_MSG(rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
                      "all the matrices in the RSM nonterminal matrix decomposition "
                      "should have the same dimensions and be square")
    }
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (B_call[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&rows, B_call[i]));
        GRB_TRY(GrB_Matrix_ncols(&cols, B_call[i]));

        LG_ASSERT_MSG(rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
                      "all the matrices in the RSM call matrix decomposition "
                      "should have the same dimensions and be square")
    }
    for (size_t i = 0; i < num_nonterminals; ++i)
    {
        if (B_ret[i] == GrB_NULL)
            continue;

        GRB_TRY(GrB_Matrix_nrows(&rows, B_ret[i]));
        GRB_TRY(GrB_Matrix_ncols(&cols, B_ret[i]));

        LG_ASSERT_MSG(rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
                      "all the matrices in the RSM return matrix decomposition "
                      "should have the same dimensions and be square")
    }

    // Check source nodes in the graph
    for (size_t i = 0; i < ns; ++i)
    {
        GrB_Index s = Source[i];
        LG_ASSERT_MSG(s < ng, GrB_INVALID_INDEX, "invalid graph source node");
    }
    // Check starting states of the RSM
    for (size_t i = 0; i < nqs; ++i)
    {
        GrB_Index qs = QS[i];
        LG_ASSERT_MSG(qs < nr, GrB_INVALID_INDEX, "invalid RSM starting state");
    }
    // Check final states of the RSM
    for (size_t i = 0; i < nqf; ++i)
    {
        GrB_Index qf = QF[i];
        LG_ASSERT_MSG(qf < nr, GrB_INVALID_INDEX, "invalid RSM final state");
    }

    //--------------------------------------------------------------------------
    // allocate derived graph-edge matrices A_S[s]  (ng x ng, initially empty)
    //--------------------------------------------------------------------------

    LG_TRY(LAGraph_Malloc((void **)&A_S, num_nonterminals, sizeof(GrB_Matrix), msg))
    for (size_t i = 0; i < num_nonterminals; ++i)
        GRB_TRY(GrB_Matrix_new(&A_S[i], GrB_BOOL, ng, ng));

    // -------------------------------------------------------------------------
    // initialization
    // -------------------------------------------------------------------------

    GRB_TRY(GrB_Vector_new(reachable, GrB_BOOL, ng));
    GRB_TRY(GrB_Vector_new(&final_reducer, GrB_BOOL, nr));

    // final_reducer[QF] = true
    GrB_assign(final_reducer, GrB_NULL, GrB_NULL, true, QF, nqf, GrB_NULL);

    // frontier / visited matrices  (nr x ng)
    GRB_TRY(GrB_Matrix_new(&frontier, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&next_frontier, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&front_temp1, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&front_temp2, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&P, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&K, GrB_BOOL, nr, ng))

    GRB_TRY(GrB_Matrix_new(&M_ret, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&M_call, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&M_call_new, GrB_BOOL, nr, ng));
    GRB_TRY(GrB_Matrix_new(&G_S_new, GrB_BOOL, ng, ng));
    GRB_TRY(GrB_Matrix_new(&step1, GrB_BOOL, ng, nr));

    // Seed next_frontier, P, K from (start-states × source-vertices)
    GrB_assign(next_frontier, GrB_NULL, GrB_NULL, true, QS, nqs, Source, ns, GrB_NULL);
    GrB_assign(P, GrB_NULL, GrB_NULL, true, QS, nqs, Source, ns, GrB_NULL);
    GrB_assign(K, GrB_NULL, GrB_NULL, true, QS, nqs, Source, ns, GrB_NULL);

    GRB_TRY(GrB_Matrix_nvals(&states, next_frontier));

    int iteration = 0;

    // Main loop
    while (states != 0)
    {
        iteration++;

        GrB_Matrix old_frontier = frontier;
        frontier = next_frontier;
        next_frontier = old_frontier;
        GRB_TRY(GrB_Matrix_clear(next_frontier));

        //----------------------------------------------------------------------
        // FORWARD PHASE
        //----------------------------------------------------------------------

        // Terminals
        for (size_t i = 0; i < num_terminals; ++i)
        {
            if (A_a[i] == GrB_NULL || B_a[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(GrB_Matrix_clear(front_temp2));

            // front_temp1  = N_a[a]^T * M
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, B_a[i],
                            frontier, GrB_DESC_T0));

            // front_temp2 = front_temp1  * G_a[a]
            GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, front_temp1,
                            A_a[i], GrB_NULL));

            // M_new |= front_temp2 & ~P
            GRB_TRY(GrB_eWiseAdd(next_frontier, P, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 next_frontier, front_temp2, GrB_DESC_SC));
        }

        // Nonterminals (derived A_S edges)
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (B_S[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(GrB_Matrix_clear(front_temp2));

            // front_temp1 = N_S[i]^T * frontier
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_S[i],
                            frontier, GrB_DESC_T0));

            // front_temp2 = front_temp1 * A_S[i]
            GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                            front_temp1, A_S[i], GrB_NULL));

            // next_frontier |= symbol_frontier & ~P
            GRB_TRY(GrB_eWiseAdd(next_frontier, P, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 next_frontier, front_temp2, GrB_DESC_SC));
        }

        // Call edges
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (B_call[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(front_temp1));

            // front_temp1 = B_call[s]^T * frontier
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_call[i],
                            frontier, GrB_DESC_T0));

            // next_frontier |= front_temp1 & ~P
            GRB_TRY(GrB_eWiseAdd(next_frontier, P, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 next_frontier, front_temp1, GrB_DESC_SC));
        }

        //----------------------------------------------------------------------
        // BACKWARD PHASE – match return edges back to their call sites
        //----------------------------------------------------------------------

        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (B_ret[i] == GrB_NULL || B_call[i] == GrB_NULL || B_S[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(M_ret));
            GRB_TRY(GrB_Matrix_clear(M_call));

            // M_ret = B_ret[s]^T * frontier  (positions that hit a return state)
            GRB_TRY(GrB_mxm(M_ret, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_ret[i],
                            frontier, GrB_DESC_T0));

            GrB_Index ret_nvals;
            GRB_TRY(GrB_Matrix_nvals(&ret_nvals, M_ret));
            if (ret_nvals == 0)
                continue;

            // M_call = (B_ret[s] * M_ret) & frontier
            GRB_TRY(GrB_mxm(M_call, frontier, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_ret[i], M_ret,
                            GrB_DESC_S));

            // walk backwards through visited set P to find the
            // matching call-site positions.

            GrB_Index call_nvals;
            GRB_TRY(GrB_Matrix_nvals(&call_nvals, M_call));
            while (call_nvals > 0)
            {
                GRB_TRY(GrB_Matrix_clear(M_call_new));

                // Backward terminal step: M_call_new |= (N_a[i] * M_call * G_a[i]^T) & P
                for (size_t i = 0; i < num_terminals; ++i)
                {
                    if (A_a[i] == GrB_NULL || B_a[i] == GrB_NULL)
                        continue;

                    GRB_TRY(GrB_Matrix_clear(front_temp1));
                    GRB_TRY(GrB_Matrix_clear(front_temp2));

                    // front_temp1  = N_a[a] * M_call
                    GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                    B_a[i], M_call, GrB_NULL));

                    // front_temp2 = front_temp1  * G_a[a]^T
                    GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                    front_temp1, A_a[i], GrB_DESC_T1));

                    // M_call_new |= front_temp2 & P
                    GRB_TRY(GrB_eWiseAdd(M_call_new, P, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                         M_call_new, front_temp2, GrB_DESC_S));
                }
                // Backward nonterminal step: M_call_new |= (N_S[i] * M_call * G_S[i]^T) & P
                for (size_t i = 0; i < num_nonterminals; ++i)
                {
                    if (A_S[i] == GrB_NULL)
                        continue;

                    GRB_TRY(GrB_Matrix_clear(front_temp1));
                    GRB_TRY(GrB_Matrix_clear(front_temp2));

                    // front_temp1 = N_S[S] * M_call
                    GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                    B_S[i], M_call, GrB_NULL));

                    // front_temp2 = front_temp1 * G_S[S]^T
                    GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                                    front_temp1, A_S[i], GrB_DESC_T1));

                    // M_call_new |= front_temp2 & P
                    GRB_TRY(GrB_eWiseAdd(M_call_new, P, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                         M_call_new, front_temp2, GrB_DESC_S));
                }

                GRB_TRY(GrB_Matrix_nvals(&call_nvals, M_call_new));
                if (call_nvals == 0)
                    break;

                GRB_TRY(GrB_Matrix_clear(M_call));
                GRB_TRY(GrB_eWiseAdd(M_call, GrB_NULL, GrB_NULL, GrB_LOR, M_call, M_call_new,
                                     GrB_NULL));
            }

            // need to check if we reach call pos
            // M_call = (B_call[s] * M_call) & P
            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_call[i],
                            M_call, GrB_NULL));
            GRB_TRY(GrB_Matrix_clear(M_call));
            GRB_TRY(GrB_eWiseAdd(M_call, P, GrB_NULL, GrB_LOR, M_call, front_temp1, GrB_DESC_S));

            GRB_TRY(GrB_Matrix_nvals(&call_nvals, M_call));
            if (call_nvals == 0)
                continue;

            // we found some call positions
            // but still not sure if they are correct ones
            // G_S_new = (M_call.T @ N_S[S] @ M_ret) & ~G_S[S]
            GRB_TRY(GrB_Matrix_clear(step1));
            GRB_TRY(GrB_Matrix_clear(front_temp2));
            GRB_TRY(GrB_Matrix_clear(G_S_new));

            GRB_TRY(GrB_mxm(step1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, M_call, B_S[i],
                            GrB_DESC_T0)); // step1 is ng*nr

            GRB_TRY(GrB_Matrix_clear(G_S_new));
            GRB_TRY(GrB_mxm(G_S_new, A_S[i], GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, step1, M_ret,
                            GrB_DESC_SC));

            GrB_Index new_edge_nvals;
            GRB_TRY(GrB_Matrix_nvals(&new_edge_nvals, G_S_new));
            if (new_edge_nvals == 0)
                continue;

            /*
            call positions are correct
            we derived some edges for graph
            now we need to traverse this part again */

            // G_S[S] |= G_S_new
            GRB_TRY(GrB_eWiseAdd(A_S[i], GrB_NULL, GrB_NULL, GrB_LOR, A_S[i], G_S_new, GrB_NULL));

            // P &= ~M_ret
            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(
                GrB_eWiseAdd(front_temp1, GrB_NULL, GrB_NULL, GrB_LOR, front_temp1, P, GrB_NULL));
            GRB_TRY(GrB_Matrix_clear(P));
            GRB_TRY(GrB_eWiseAdd(P, M_ret, GrB_NULL, GrB_LOR, P, front_temp1, GrB_DESC_SC));

            // next_frontier |= M_call
            GRB_TRY(GrB_eWiseAdd(next_frontier, GrB_NULL, GrB_NULL, GrB_LOR, next_frontier, M_call,
                                 GrB_NULL));
        }

        //----------------------------------------------------------------------
        // ADVANCE: P |= next_frontier
        //----------------------------------------------------------------------

        GRB_TRY(GrB_eWiseAdd(P, GrB_NULL, GrB_NULL, GrB_LOR, P, next_frontier, GrB_NULL));
        GRB_TRY(GrB_Matrix_nvals(&states, next_frontier));

        if (iteration > 3000)
        {
            printf("Warning: Maximum iterations reached\n");
            break;
        }
    }

    GRB_TRY(GrB_Matrix_clear(next_frontier));
    GrB_assign(next_frontier, GrB_NULL, GrB_NULL, true, QS, nqs, Source, ns, GrB_NULL);

    GRB_TRY(GrB_Matrix_nvals(&states, next_frontier));
    // second loop
    while (states != 0)
    {
        iteration++;
        GrB_Matrix old_frontier = frontier;
        frontier = next_frontier;
        next_frontier = old_frontier;
        GRB_TRY(GrB_Matrix_clear(next_frontier));

        // Terminals
        for (size_t i = 0; i < num_terminals; ++i)
        {
            if (A_a[i] == GrB_NULL || B_a[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(GrB_Matrix_clear(front_temp2));

            // front_temp1  = N_a[a]^T * M
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, B_a[i],
                            frontier, GrB_DESC_T0));

            // front_temp2 = front_temp1  * G_a[a]
            GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, front_temp1,
                            A_a[i], GrB_NULL));

            // M_new |= front_temp2 & ~K
            GRB_TRY(GrB_eWiseAdd(next_frontier, K, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 next_frontier, front_temp2, GrB_DESC_SC));
        }

        // Nonterminals (derived A_S edges)
        for (size_t i = 0; i < num_nonterminals; ++i)
        {
            if (B_S[i] == GrB_NULL)
                continue;

            GRB_TRY(GrB_Matrix_clear(front_temp1));
            GRB_TRY(GrB_Matrix_clear(front_temp2));

            // front_temp1 = N_S[i]^T * frontier
            GRB_TRY(GrB_mxm(front_temp1, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, B_S[i],
                            frontier, GrB_DESC_T0));

            // front_temp2 = front_temp1 * G_S[i]
            GRB_TRY(GrB_mxm(front_temp2, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL,
                            front_temp1, A_S[i], GrB_NULL));

            // next_frontier |= symbol_frontier & ~K
            GRB_TRY(GrB_eWiseAdd(next_frontier, K, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL,
                                 next_frontier, front_temp2, GrB_DESC_SC));
        }

        //----------------------------------------------------------------------
        // ADVANCE: P |= next_frontier
        //----------------------------------------------------------------------

        GRB_TRY(GrB_eWiseAdd(K, GrB_NULL, GrB_NULL, GrB_LOR, K, next_frontier, GrB_NULL));
        GRB_TRY(GrB_Matrix_nvals(&states, next_frontier));

        if (iteration > 3000)
        {
            printf("Warning: Maximum iterations reached\n");
            break;
        }
    }

    GRB_TRY(GrB_vxm(*reachable, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, final_reducer, K,
                    GrB_NULL))

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
