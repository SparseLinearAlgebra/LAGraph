//------------------------------------------------------------------------------
// LAGraph_RSM_reachability.c
//------------------------------------------------------------------------------

// For and edge-labelled directed graph the algorithm computes the set of nodes
// for which these conditions are held:
// * The node is reachable by a path from one of the source nodes.
// * The concatenation of the labels over this path's edges is a word from the
//   specified context-free language
//
// The context-free constraints are specified by a recursive state machine (RSM)
// over a subset of the graph edge labels. The algorithm assumes the labels are
// enumerated from 0 to some nl-1. The input graph and the RSM are defined by
// adjacency matrix decomposition. They are represented by arrays of graphs G and
// R both of length

#include "GraphBLAS.h"
#define LG_FREE_WORK \
{                    \
    GrB_free (&frontier) ; \
    GrB_free (&next_frontier) ; \
    GrB_free (&symbol_frontier) ; \
    GrB_free (&visited) ; \
    GrB_free (&final_reducer) ; \
    LAGraph_Free ((void **) &A, GrB_NULL) ; \
    LAGraph_Free ((void **) &B, GrB_NULL) ; \
    LAGraph_Free ((void **) &BT, GrB_NULL) ; \
    LAGraph_Free ((void **) &B_call, GrB_NULL) ; \
    LAGraph_Free ((void **) &BT_call, GrB_NULL) ; \
}

#define LG_FREE_ALL        \
{                          \
    LG_FREE_WORK ;         \
    GrB_free (reachable) ; \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_RSM_reachability
(
    // output:
    GrB_Vector *reachable, // reachable(i) = true if node i is reachable
                           // from one of the starting nodes by a path
                           // satisfying regular constraints

    LAGraph_Graph *R,      // input recursive state machine's
                           // adjacency matrix decomposition
    size_t nl,             // total label count, # of matrices graph and
                           // RSM adjacency matrix decomposition
    const GrB_Index *QS,   // starting states in RSM
    size_t nqs,            // number of starting states in RSM
    const GrB_Index *QF,   // final states in RSM
    size_t nqf,            // number of final states in RSM

    LAGraph_Graph *R_call, // recursive state machine's call/return nr*nr
                           // matrix, where (i, j) = 1 if there is transition
                           // from i state of one automata to j starting 
                           // state of another automata
    LAGraph_Graph *G,      // input graph adjacency matrix decomposition
    const GrB_Index *S,    // source vertices to start searching paths
    size_t ns,             // number of source vertices
    char *msg              // LAGraph output message
) 
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Matrix frontier = GrB_NULL ;        // traversal frontier representing
                                            // correspondence between RSM states
                                            // and graph vertices
    GrB_Matrix symbol_frontier = GrB_NULL ; // part of the new frontier for the 
                                            // specific label
    GrB_Matrix next_frontier = GrB_NULL ;   // frontier value on the next
                                            // traversal step
    GrB_Matrix visited = GrB_NULL ;         // visited pairs (state, vertex)

    GrB_Vector final_reducer = GrB_NULL ;   // auxiliary vector for reducing the
                                            // visited matrix to an answer

    GrB_Index ng = 0 ;      // # nodes in the graph
    GrB_Index nr = 0 ;      // # states in the RSM
    GrB_Index states = ns ; // # pairs in the current
                            // correspondence between the graph and
                            // the RSM

    GrB_Index rows = 0 ;    // utility matrix row count
    GrB_Index cols = 0 ;    // utility matrix column count
    GrB_Index vals = 0 ;    // utility matrix valiue count

    GrB_Matrix *A = GrB_NULL ;
    GrB_Matrix *B = GrB_NULL ;
    GrB_Matrix *BT = GrB_NULL ; 

    GrB_Matrix *B_call = GrB_NULL ;
    GrB_Matrix *BT_call = GrB_NULL ;

    LG_ASSERT (reachable != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (R != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (G != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (S != GrB_NULL, GrB_NULL_POINTER) ;

    LG_ASSERT (R_call != GrB_NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (*R_call != GrB_NULL, GrB_NULL_POINTER) ;

    (*reachable) = GrB_NULL ;

    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (G[i] == GrB_NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (G[i], msg)) ;
    }

    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (R[i] == GrB_NULL) continue;
        LG_TRY (LAGraph_CheckGraph (R[i], msg)) ;
    }

    LG_TRY (LAGraph_Malloc ((void **) &A, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (G[i] == GrB_NULL)
        {
            A[i] = GrB_NULL ;
            continue ;
        }

        A[i] = G[i]->A ;
    }

    LG_TRY (LAGraph_Malloc((void **) &B, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc((void **) &BT, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; ++i)
    {
        BT[i] = GrB_NULL ;

        if (R[i] == GrB_NULL)
        {
            B[i] = GrB_NULL ;
            continue ;
        }

        B[i] = R[i]->A ;
        if (R[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            BT[i] = B[i] ;
        } else 
        {
            // BT[i] could be NULL and the matrix will be transposed by a
            // descriptor
            BT[i] = R[i]->AT ;
        }
    }
    


    LG_TRY (LAGraph_Malloc((void **) &B_call, 1, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc((void **) &BT_call, 1, sizeof (GrB_Matrix), msg)) ;

    *B_call = (*R_call)->A ;
    *BT_call = (*R_call)->AT ;



    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (A[i] == GrB_NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&ng, A[i])) ;
        break ;
    } 

    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (B[i] == GrB_NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&nr, B[i])) ;
        break ;
    }

    // Check all the matrices in graph adjacency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (A[i] == GrB_NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&rows, A[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, A[i])) ;

        LG_ASSERT_MSG (rows == ng && cols == ng, LAGRAPH_NOT_CACHED,
            "all the matrices in the graph adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check all the matrices in RSM adjancency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; ++i)
    {
        if (B[i] == GrB_NULL) continue ;
        
        GRB_TRY (GrB_Matrix_nrows (&rows, B[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, B[i])) ;

        LG_ASSERT_MSG (rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
            "all the matrices in the RSM adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check size of RSM call/return matrix

    GRB_TRY (GrB_Matrix_nrows (&rows, *B_call)) ;
    GRB_TRY (GrB_Matrix_ncols (&cols, *B_call)) ;

    LG_ASSERT_MSG (rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
        "RSM's call/return matrix should have nr*nr size") ;

    // Check source nodes in the graph
    for (size_t i = 0 ; i < ns ; ++i)
    {
        GrB_Index s = S[i] ;
        LG_ASSERT_MSG(s < ng, GrB_INVALID_INDEX, "invalid graph source node") ;
    }

    // Check starting states of the RSM
    for (size_t i = 0 ; i < nqs ; ++i)
    {
        GrB_Index qs = QS[i] ;
        LG_ASSERT_MSG (qs < nr, GrB_INVALID_INDEX,
            "invalid RSM starting state") ;
    }

    // Check final states of the RSM
    for (size_t i = 0 ; i < nqf ; ++i)
    {
        GrB_Index qf = QF[i] ;
        LG_ASSERT_MSG (qf < nr, GrB_INVALID_INDEX, "invalid RSM final state") ;
    }

    // -------------------------------------------------------------------------
    // initialization
    // -------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (reachable, GrB_BOOL, ng)) ;
    GRB_TRY (GrB_Vector_new (&final_reducer, GrB_BOOL, nr)) ;

    // Initialize matrix for reducing the result
    GrB_assign (final_reducer, GrB_NULL, GrB_NULL, true, QF, nqf, GrB_NULL) ;

    GRB_TRY (GrB_Matrix_new (&frontier, GrB_BOOL, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&next_frontier, GrB_BOOL, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&symbol_frontier, GrB_BOOL, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&visited, GrB_BOOL, nr, ng)) ;

    // Initialize frontier with the source nodes
    GrB_assign (next_frontier, GrB_NULL, GrB_NULL, true, QS, nqs, S, ns, GrB_NULL) ;
    GrB_assign (visited, GrB_NULL, GrB_NULL, true, QS, nqs, S, ns, GrB_NULL) ;
    
    // Main loop
    while (states != 0)
    {
        GrB_Matrix old_frontier = frontier ;
        frontier = next_frontier ;
        next_frontier = old_frontier ;

        GRB_TRY (GrB_Matrix_clear (next_frontier)) ;

        // Obtain a new relation between the RSM states and the graph nodes
        for (size_t i = 0 ; i < nl ; ++i)
        {
            if (A[i] == GrB_NULL || B[i] == GrB_NULL) continue ;

            // Traverse the RSM
            // Try to use a provided transposed matrix or use the descriptor

            GRB_TRY (GrB_Matrix_clear (symbol_frontier)) ;

            if (BT[i] != GrB_NULL)
            {
                GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_LOR,
                    GrB_LOR_LAND_SEMIRING_BOOL, BT[i], frontier, GrB_NULL)) ;
            }
            else
            {
                GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_LOR,
                    GrB_LOR_LAND_SEMIRING_BOOL, B[i], frontier, GrB_DESC_T0)) ;
            }

            if ((*BT_call) != GrB_NULL) 
            {
                GRB_TRY (GrB_mxm (next_frontier, visited, GrB_LOR,
                    GrB_LOR_LAND_SEMIRING_BOOL, *BT_call, frontier, GrB_DESC_SC)) ; // не уверен по поводу дескриптора
            }
            else
            {
                GRB_TRY (GrB_mxm (next_frontier, visited, GrB_LOR,
                    GrB_LOR_LAND_SEMIRING_BOOL, *B_call, frontier, GrB_DESC_SCT0)) ; // не уверен по поводу дескриптора
            }

            // я реально не понимаю как это работает, но окей
            // Traverse the graph
            GRB_TRY (GrB_mxm (next_frontier, visited, GrB_LOR,
                GrB_LOR_LAND_SEMIRING_BOOL, symbol_frontier, A[i], GrB_DESC_SC)) ;
        }

        // я реально не понимаю как это работает, но окей
        // Accumulate the new state <-> node correspondence
        GRB_TRY (GrB_assign (visited, visited, GrB_NULL, next_frontier,
            GrB_ALL, nr, GrB_ALL, ng, GrB_DESC_SC)) ;

        GRB_TRY (GrB_Matrix_nvals (&states, next_frontier)) ;
    }

    // Extract the nodes matching the final RSM states
    GRB_TRY (GrB_vxm (*reachable, GrB_NULL, GrB_NULL,
        GrB_LOR_LAND_SEMIRING_BOOL, final_reducer, visited, GrB_NULL)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
