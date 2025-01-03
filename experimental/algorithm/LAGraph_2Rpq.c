//
// Different RPQ semantics
//

#define LG_FREE_WORK                            \
{                                               \
}
#define LG_FREE_ALL                         \
{                                           \
        LG_FREE_WORK ;                          \
}

#include "LG_internal.h"
#include "LAGraphX.h"
#include <assert.h>
#include <limits.h>

#define PATH_LIMIT 100000

typedef struct {
    Path paths[QUICK_PATH_COUNT];
    size_t path_count;
    Path *extra_paths;
} MultiplePaths ;

MultiplePaths multiple_paths_identity ;

void Path_print (const Path *x)
{
    if (x->vertex_count == 0)
    {
        printf ("empty path \n") ;
        return ;
    }

    for (size_t i = 0 ; i < x->vertex_count ; i++)
    {
        // Increase the vertex by 1 since usually user expects the same
        // numbering as in the input determined by MTX file in which the
        // entries are enumerated starting from 1.
        printf ("(%ld)", (i < QUICK_PATH_LENGTH ? x->vertices[i] : x->extra_vertices[i - QUICK_PATH_LENGTH]) + 1) ;

        if (i != x->vertex_count - 1)
        {
            printf ("-") ;
        }
    }

    printf ("\n") ;
}

static void MultiplePaths_print (const MultiplePaths *x)
{
    printf("Multiple paths:\n") ;
    for (size_t i = 0 ; i < x->path_count ; i++)
    {

        printf("\t Path %ld: ", i) ;
        Path_print (&x->paths[i]) ;
    }
    printf("\n") ;
}

GrB_Type multiple_paths ;
GrB_BinaryOp combine_multiple_paths_op ;
GrB_Monoid combine_multiple_paths ;
GrB_BinaryOp first_multiple_paths ;
GrB_BinaryOp second_multiple_paths ;
GrB_Semiring first_combine_multiple_paths ;
GrB_Semiring second_combine_multiple_paths ;
GrB_IndexUnaryOp extend_multiple_paths ;
GrB_IndexUnaryOp extend_multiple_simple ;
GrB_IndexUnaryOp extend_multiple_trails ;

void first_multiple_paths_f(MultiplePaths *z, MultiplePaths *x, bool *_y)
{
    *z = *x;
}
void second_multiple_paths_f(MultiplePaths *z, bool *_x, MultiplePaths *y)
{
    *z = *y;
}

void combine_multiple_paths_f(MultiplePaths *z, const MultiplePaths *x, const MultiplePaths *y)
{
    z->path_count = x->path_count + y->path_count ;
    assert (z->path_count < QUICK_PATH_COUNT) ;

    for (size_t i = 0 ; i < x->path_count ; i++)
    {
        z->paths[i] = x->paths[i] ;
    }

    for (size_t i = 0 ; i < y->path_count ; i++)
    {
        z->paths[x->path_count + i] = y->paths[i] ;
    }

    // TODO: Support more than QUICK_PATH_COUNT paths.
}

static inline void path_extend(Path *path, Vertex vertex)
{
    if (path->vertex_count == 0)
    {
        return ;
    }

    if (path->vertex_count < QUICK_PATH_LENGTH)
    {
        path->vertices[path->vertex_count++] = vertex ;
    }
    else
    {
        if (path->extra_vertices == NULL)
        {
            LG_TRY (LAGraph_Calloc ((void **) &path->extra_vertices, 64, sizeof (Vertex), NULL)) ;
        }
        path->extra_vertices [(path->vertex_count++) - QUICK_PATH_LENGTH] = vertex ;
    }

    // TODO: Support more than QUICK_PATH_LENGTH vertices.
}

static inline bool path_is_empty(Path *path)
{
    return path->vertex_count == 0;
}


static inline void multiple_paths_append(MultiplePaths *multiple_paths, const Path *path)
{
    multiple_paths->paths[multiple_paths->path_count++] = *path ;

    // TODO: Support more than QUICK_PATH_COUNT paths.
}

//
// ALL PATHS.
//

// NB: Using this semantic without a length limit makes the code behave like a
// procedure for searching all paths satisfying the constraints.
// It means it may not finish if there is loops.

void extend_multiple_paths_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *_y)
{
    /*if (z != x)
        for (size_t i = 0 ; i < x->path_count ; i++)
        {
            multiple_paths_append(z, &x->paths[i]) ;
            path_extend (&z->paths[i], col) ;
        }
    {*/
        for (size_t i = 0 ; i < z->path_count ; i++)
        {
            Path *path = &z->paths[i] ;
            path_extend (&z->paths[i], col) ;
        }
    //}
}

//
// ALL SIMPLE
//

static inline bool path_extending_will_add_repeated_non_starting_vertex(const Path *path, Vertex vertex)
{
    if (path->vertex_count <= 1) 
    {
        return false ;
    }

    for (size_t i = 1 ; i < path->vertex_count ; i++)
    {
        if (path->vertices[i] == vertex)
        {
            return true ;
        }
    }

    Vertex last_vertex = path->vertices[path->vertex_count - 1] ;

    return path->vertices[0] == last_vertex;
}

void extend_multiple_simple_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *_y)
{
    /*if (z != x)
    {
        for (size_t i = 0 ; i < x->path_count ; i++)
        {
            const Path *path = &x->paths[i] ;
            if (path_has_loop_at_end (path))
            {
                continue;
            }

            multiple_paths_append(z, path) ;
            path_extend (&z->paths[i], col) ;
        }
    }
    else
    {*/
        for (size_t i = 0 ; i < z->path_count ; i++)
        {
            Path *path = &z->paths[i] ;
            if (path_extending_will_add_repeated_non_starting_vertex (path, col))
            {
                path->vertex_count = 0 ;
                continue ;
            }

            path_extend (&z->paths[i], col) ;
        }
    //}
}

//
// ALL TRAILS
//

static inline bool path_extending_will_add_repeated_edge(const Path *path, Vertex vertex_2)
{
    if (path->vertex_count == 0)
    {
        return false ;
    }

    // We identify edges as pairs of vertices.
    Vertex vertex_1 = path->vertices[path->vertex_count - 1] ;

    for (size_t i = 0 ; i < path->vertex_count - 1; i++)
    {
        if (path->vertices[i] == vertex_1 && path->vertices[i + 1] == vertex_2)
        {
            return true ;
        }
    }

    return false ;
}
void extend_multiple_trails_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *y)
{
    /*if (z != x)
    {
        z->path_count = x->path_count ;

        for (size_t i = 0 ; i < x->path_count ; i++)
        {
            const Path *path = &x->paths[i] ;
            if (path_extending_will_add_repeated_edge (path, col))
            {
                continue ;
            }

            multiple_paths_append(z, path) ;
            path_extend (&z->paths[i], col) ;
        }
    }
    else
    {*/
        for (size_t i = 0 ; i < z->path_count ; i++)
        {
            Path *path = &z->paths[i] ;
            if (path_extending_will_add_repeated_edge (path, col))
            {
                path->vertex_count = 0 ;
                continue ;
            }

            path_extend (&z->paths[i], col) ;
        }
    //}
}

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&next_frontier) ;                 \
    GrB_free (&symbol_frontier) ;               \
    GrB_free (&final_reducer) ;                 \
    LAGraph_Free ((void **) &A, NULL) ;         \
    LAGraph_Free ((void **) &B, NULL) ;         \
    LAGraph_Free ((void **) &BT, NULL) ;        \
}

#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    LAGraph_Free ((void **) paths, NULL) ;  \
}

static int LAGraph_2Rpq
(
    // output:
    Path **paths,               // simple paths from one of the starting
                                // nodes satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    uint64_t limit,             // maximum path count
    char *msg,                  // LAGraph output message
    GrB_IndexUnaryOp op         // index unary op for a specific semantic
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Matrix frontier = NULL ;         // traversal frontier representing
                                         // correspondence between NFA states
                                         // and graph vertices
    GrB_Matrix symbol_frontier = NULL ;  // part of the new frontier for the
                                         // specific label
    GrB_Matrix next_frontier = NULL ;    // frontier value on the next
                                         // traversal step
    GrB_Vector final_reducer = NULL ;    // auxiliary vector for reducing the
                                         // visited matrix to an answer

    GrB_Index ng = 0 ;                   // # nodes in the graph
    GrB_Index nr = 0 ;                   // # states in the NFA
    GrB_Index nv = 0 ;                   // # pair count in the frontier
    GrB_Index states = ns ;              // # pairs in the current
                                         // correspondence between the graph and
                                         // the NFA

    GrB_Index rows = 0 ;                 // utility matrix row count
    GrB_Index cols = 0 ;                 // utility matrix column count
    GrB_Index vals = 0 ;                 // utility matrix value count

    // TODO: This names might be too short.
    GrB_Semiring sr1 = first_combine_multiple_paths ;
    GrB_Semiring sr2 = second_combine_multiple_paths ;
    GrB_BinaryOp acc = combine_multiple_paths_op ;

    GrB_Matrix *A = NULL ;
    GrB_Matrix *AT = NULL ;
    GrB_Matrix *B = NULL ;
    GrB_Matrix *BT = NULL ;

    LG_ASSERT (paths != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (path_count != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (G != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (R != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (S != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (op != NULL, GrB_NULL_POINTER) ;

    (*paths) = NULL ;
    (*path_count) = 0 ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (G[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (G[i], msg)) ;
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (R[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (R[i], msg)) ;
    }

    LG_TRY (LAGraph_Malloc ((void **) &A, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &AT, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (G[i] == NULL)
        {
            A[i] = NULL ;
            AT[i] = NULL ;
            continue ;
        }

        A[i] = G[i]->A ;
        if (G[i]->kind == LAGraph_ADJACENCY_UNDIRECTED ||
            G[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            AT[i] = A[i] ;
        }
        else
        {
            // AT[i] could be NULL and the matrix will be transposed by a
            // descriptor
            AT[i] = G[i]->AT ;
        }
    }

    LG_TRY (LAGraph_Malloc ((void **) &B, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &BT, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        BT[i] = NULL ;

        if (R[i] == NULL)
        {
            B[i] = NULL ;
            BT[i] = NULL ;
            continue ;
        }

        B[i] = R[i]->A ;
        if (R[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            BT[i] = B[i] ;
        }
        else
        {
            // BT[i] could be NULL and the matrix will be transposed by a
            // descriptor
            BT[i] = R[i]->AT ;
        }
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&ng, A[i])) ;
        break ;
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (B[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&nr, B[i])) ;
        break ;
    }

    // Check all the matrices in graph adjacency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; i++)
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&rows, A[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, A[i])) ;

        LG_ASSERT_MSG (rows == ng && cols == ng, LAGRAPH_NOT_CACHED,
            "all the matrices in the graph adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check all the matrices in NFA adjacency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; i++)
    {
        if (B[i] == NULL) continue ;

        GrB_Index rows = 0 ;
        GrB_Index cols = 0 ;

        GRB_TRY (GrB_Matrix_nrows (&rows, B[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, B[i])) ;

        LG_ASSERT_MSG (rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
            "all the matrices in the NFA adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check source nodes in the graph
    for (size_t i = 0 ; i < ns ; i++)
    {
        GrB_Index s = S [i] ;
        LG_ASSERT_MSG (s < ng, GrB_INVALID_INDEX, "invalid graph source node") ;
    }

    // Check starting states of the NFA
    for (size_t i = 0 ; i < nqs ; i++)
    {
        GrB_Index qs = QS [i] ;
        LG_ASSERT_MSG (qs < nr, GrB_INVALID_INDEX,
            "invalid NFA starting state") ;
    }

    // Check final states of the NFA
    for (size_t i = 0 ; i < nqf ; i++)
    {
        GrB_Index qf = QF [i] ;
        LG_ASSERT_MSG (qf < nr, GrB_INVALID_INDEX, "invalid NFA final state") ;
    }

    // -------------------------------------------------------------------------
    // initialization
    // -------------------------------------------------------------------------

    LG_TRY (LAGraph_Calloc ((void **) paths, PATH_LIMIT, sizeof (Path), msg)) ;

    GRB_TRY (GrB_Vector_new (&final_reducer, GrB_BOOL, nr)) ;

    // Initialize matrix for reducing the result
    GRB_TRY (GrB_assign (final_reducer, NULL, NULL, true, QF, nqf, NULL)) ;

    GRB_TRY (GrB_Matrix_new (&next_frontier, multiple_paths, nr, ng)) ;

    // Initialize frontier with the source nodes

    for (size_t i = 0 ; i < ns ; i++)
    {
        GrB_Index s = S[i] ;
        MultiplePaths value = {
            .paths = {
                {
                    .vertices = { s },
                    .vertex_count = 1
                }
            },
            .path_count = 1
        };

        for (size_t j = 0 ; j < nqs ; j++)
        {
            GrB_Index qs = QS[j] ;

            GRB_TRY (GrB_Matrix_setElement_UDT (next_frontier, &value, qs, s)) ;
        }
    }

    // Initialize a few utility matrices
    GRB_TRY (GrB_Matrix_new (&frontier, multiple_paths, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&symbol_frontier, multiple_paths, nr, ng)) ;

    // Main loop
    while (true)
    {
        //printf("Iteration\n");
        GrB_Index nvals = 0 ;
        GRB_TRY (GrB_Matrix_nvals (&nvals, next_frontier)) ;

        MultiplePaths *X ;
        GrB_Index *I ;
        bool had_non_empty_path = false ;

        //MultiplePaths *X;
        LG_TRY (LAGraph_Calloc ((void **) &X, nvals, sizeof (MultiplePaths), msg)) ;
        LG_TRY (LAGraph_Calloc ((void **) &I, nvals, sizeof (GrB_Index), msg)) ;

        // TODO: Change to a generic call.
        GRB_TRY (GrB_Matrix_extractTuples_UDT (I, GrB_NULL, (void**) X, &nvals, next_frontier)) ;
        //printf("Next frontier with %d entries\n", nvals);

        for (size_t i = 0 ; i < nvals ; i++)
        {
            for (size_t j = 0 ; j < X[i].path_count ; j++)
            {
                if (!path_is_empty(&X[i].paths[j]))
                {
                    had_non_empty_path = true;
                    break;
                }
            }

            //MultiplePaths_print (&X[i]) ;
            bool final = false ;
            for (size_t j = 0 ; j < nqf ; j++)
            {
                if (I[i] == QF[j])
                {
                    final = true ;
                    break ;
                }
            }
            //printf("Path at %ld final is %b", I[i], final) ;

            if (!final)
            {
                continue ;
            }

            //printf("Found final paths!\n");
            for (size_t j = 0 ; j < X[i].path_count && (*path_count) < limit ; j++)
            {
                const Path *path = &X[i].paths[j] ;
                if (!path_is_empty(path))
                {
                    (*paths)[(*path_count)++] = *path ;
                }
            }
        }

        if (!had_non_empty_path || (*path_count) == limit)
        {
            //printf("breaking\n");
            break;
        }

        GrB_Matrix old_frontier = frontier ;
        frontier = next_frontier ;
        next_frontier = old_frontier ;

        GRB_TRY (GrB_Matrix_clear(next_frontier)) ;

        // Obtain a new relation between the NFA states and the graph nodes
        for (size_t i = 0 ; i < nl ; i++)
        {
            if (A[i] == NULL || B[i] == NULL) continue ;

            // Traverse the NFA
            // Try to use a provided transposed matrix or use the descriptor
            if (!inverse) {
                if (BT[i] != NULL)
                {
                    GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL,
                        sr2, BT[i], frontier, GrB_DESC_R)) ;
                }
                else
                {
                    GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL,
                        sr2, B[i], frontier, GrB_DESC_RT0)) ;
                }
            } else {
                GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL, sr2, B[i], frontier, GrB_DESC_R )) ;
            }

            // TODO: Skip the iteration if symbol_frontier is already empty.

            // Traverse the graph
            if (!inverse_labels[i]) {
                if (!inverse) {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, A[i], GrB_NULL)) ;
                } else if (AT[i]) {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, AT[i], GrB_NULL)) ;
                } else {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, A[i], GrB_DESC_T1)) ;
                }
            } else {
                if (!inverse && AT[i]) {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, AT[i], GrB_NULL)) ;
                } else if (!inverse) {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, A[i], GrB_DESC_T1)) ;
                } else {
                    GRB_TRY (GrB_mxm (next_frontier, GrB_NULL, acc, sr1, symbol_frontier, A[i], GrB_NULL)) ;
                }
            }
        }

        GRB_TRY (GrB_apply (next_frontier, GrB_NULL, GrB_NULL, op, next_frontier, false, GrB_NULL)) ;

    }

    //LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}


int LAGraph_2Rpq_AllSimple      // All simple paths satisfying regular
                                // expression. Simple paths are paths without
                                // loops or the ones with the same starting
                                // and final nodes.
(
    // output:
    Path **paths,               // simple paths from one of the starting
                                // nodes satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    char *msg                   // LAGraph output message
)
{
    return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, ULLONG_MAX, msg, extend_multiple_simple) ;
}

LAGRAPHX_PUBLIC
int LAGraph_2Rpq_AllTrails      // All trails satisfying regular expression.
                                // Trails are paths without repeated edges.
(
    // output:
    Path **paths,                // trails from one of the starting nodes
                                // satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    char *msg                   // LAGraph output message
)
{
    return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, ULLONG_MAX, msg, extend_multiple_trails) ;
}

int LAGraph_2Rpq_AllPaths       // All paths satisfying regular expression
(
    // output:
    Path **paths,               // paths from one of the starting nodes
                                // satisfying regular constraints
    size_t *path_count,         // resulting path count
                                // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    uint64_t limit,             // maximum path count
    char *msg                   // LAGraph output message
    )
{
        return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, limit, msg, extend_multiple_paths) ;
}

#define LG_FREE_WORK                            \
{                                               \
}
#define LG_FREE_ALL                         \
{                                           \
        LG_FREE_WORK ;                          \
}

int LAGraph_Rpq_initialize(char *msg)
{
        GRB_TRY (GrB_Type_new (&multiple_paths, sizeof(MultiplePaths))) ;

        GRB_TRY (GrB_BinaryOp_new (&combine_multiple_paths_op, (GxB_binary_function) &combine_multiple_paths_f, multiple_paths, multiple_paths, multiple_paths)) ;
        GRB_TRY (GrB_BinaryOp_new (&first_multiple_paths, (GxB_binary_function) &first_multiple_paths_f, multiple_paths, multiple_paths, GrB_BOOL)) ;
        GRB_TRY (GrB_BinaryOp_new (&second_multiple_paths, (GxB_binary_function) &second_multiple_paths_f, multiple_paths, GrB_BOOL, multiple_paths)) ;
        GRB_TRY (GrB_Monoid_new (&combine_multiple_paths, combine_multiple_paths_op, (void*) &multiple_paths_identity)) ;

        GRB_TRY (GrB_Semiring_new (&first_combine_multiple_paths, combine_multiple_paths, first_multiple_paths)) ;
        GRB_TRY (GrB_Semiring_new (&second_combine_multiple_paths, combine_multiple_paths, second_multiple_paths)) ;

        GRB_TRY (GrB_IndexUnaryOp_new (&extend_multiple_paths, (GxB_index_unary_function) &extend_multiple_paths_f, multiple_paths, multiple_paths, GrB_BOOL)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&extend_multiple_simple, (GxB_index_unary_function) &extend_multiple_simple_f, multiple_paths, multiple_paths, GrB_BOOL)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&extend_multiple_trails, (GxB_index_unary_function) &extend_multiple_trails_f, multiple_paths, multiple_paths, GrB_BOOL)) ;
}
