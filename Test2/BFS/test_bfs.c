//------------------------------------------------------------------------------
// LAGraph/Test2/BFS/test_bfs.c: test LAGraph_BreadthFirstSearch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph2.h"

#define NTHREAD_LIST 2
// #define NTHREAD_LIST 1
#define THREAD_LIST 0

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGraph_Delete (&G, msg) ;      \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&parent) ;            \
    GrB_free (&level) ;             \
    GrB_free (&SourceNodes) ;       \
}

#define LAGraph_CATCH(status)                                               \
{                                                                           \
    printf ("LAGraph error: %s line: %d, status: %d: %s\n", __FILE__,       \
        __LINE__, status, msg) ;                                            \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}

#define GrB_CATCH(info)                                                     \
{                                                                           \
    printf ("GraphBLAS error: %s line: %d, info: %d: %s\n", __FILE__,       \
        __LINE__, info, msg) ;                                              \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}

int main (int argc, char **argv)
{

    printf ("%s v%d.%d.%d [%s]\n",
        GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR,
        GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB,
        GxB_IMPLEMENTATION_DATE) ;

    GrB_Info info ;
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector level = NULL ;
    GrB_Vector parent = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    LAGraph_TRY (LAGraph_Init (msg)) ;
    GrB_TRY (GxB_set (GxB_BURBLE, false)) ;

    uint64_t seed = 1 ;
    FILE *f ;
    int nthreads ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ; 
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    double tpl [nthreads_max+1] ;
    double tp [nthreads_max+1] ;
    double tl [nthreads_max+1] ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;

    double t_save ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file
    //--------------------------------------------------------------------------

    if (argc > 1)
    {
        // Usage:
        //      ./test_bfs matrixfile.mtx sources.mtx
        //      ./test_bfs matrixfile.grb sources.mtx

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            LAGraph_TRY (LAGraph_BinRead (&A, filename, msg)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGraph_TRY (LAGraph_MMRead (&A, f, msg)) ;
            fclose (f) ;
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2)
        {
            filename = argv [2] ;
            printf ("sources: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Source node file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGraph_TRY (LAGraph_MMRead (&SourceNodes, f, msg)) ;
            fclose (f) ;
        }
    }
    else
    {

        // Usage:  ./test_bfs < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGraph_TRY (LAGraph_MMRead (&A, stdin, msg)) ;
    }

    // convert to boolean, pattern-only
    LAGraph_TRY (LAGraph_Pattern (&Abool, A, msg)) ;
    // LAGraph_mmwrite (Abool, stderr) ;
    GrB_free (&A) ;
    A = Abool ;
    Abool = NULL ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    GrB_Index n = nrows ;
    if (nrows != ncols) { printf ("A must be square\n") ; abort ( ) ; }
    double t_read ;
    LAGraph_TRY (LAGraph_Toc (&t_read, tic, msg)) ;
    printf ("read time: %g\n", t_read) ;

    //--------------------------------------------------------------------------
    // construct the graph
    //--------------------------------------------------------------------------

    bool A_is_symmetric =
        (nrows == 134217726 ||  // HACK for kron
         nrows == 134217728) ;  // HACK for urand

    if (A_is_symmetric)
    {
        // A is known to be symmetric
        // TODO: LAGraph_New should set G->A_pattern_is_symmetric if
        // the G->kind is LAGRAPH_ADJACENCY_UNDIRECTED
        LAGraph_TRY (LAGraph_New (&G, &A, LAGRAPH_ADJACENCY_UNDIRECTED, false,
            msg)) ;
        G->A_pattern_is_symmetric = true ;
    }
    else
    {
        // compute G->AT and determine if A has a symmetric pattern
        LAGraph_TRY (LAGraph_New (&G, &A, LAGRAPH_ADJACENCY_DIRECTED, false,
            msg)) ;
        LAGraph_TRY (LAGraph_Property_ASymmetricPattern (G, msg)) ;
        if (G->A_pattern_is_symmetric)
        {
            // if G->A has a symmetric pattern, declare the graph undirected
            // and free G->AT since it isn't needed.  The BFS only looks at
            // the pattern of A anyway.
            G->kind = LAGRAPH_ADJACENCY_UNDIRECTED ;
            GrB_TRY (GrB_Matrix_free (&(G->AT))) ;
        }
    }

    // compute G->rowdegree
    LAGraph_TRY (LAGraph_Property_RowDegree (G, msg)) ;

    // compute G->coldegree, just to test it (not needed for BFS)
    LAGraph_TRY (LAGraph_Property_ColDegree (G, msg)) ;

    LAGraph_TRY (LAGraph_DisplayGraph (G, 0, msg)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    #define NSOURCES 64

    if (SourceNodes == NULL)
    {
        GrB_TRY (GrB_Matrix_new (&SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // SourceNodes [k] = i 
            GrB_TRY (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
    }

    int64_t ntrials ;
    GrB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;

    // HACK
    // ntrials = 1 ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

    char filename [1024] ;

    //--------------------------------------------------------------------------
    // BFS
    //--------------------------------------------------------------------------

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;
        tp [nthreads] = 0 ;
        tl [nthreads] = 0 ;
        tpl [nthreads] = 0 ;
        printf ("\n------------------------------- threads: %2d\n", nthreads) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            GrB_TRY (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            double ttrial ;

            //------------------------------------------------------------------
            // BFS to compute just parent
            //------------------------------------------------------------------

            GrB_free (&parent) ;
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (NULL, &parent,
                G, s, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tp [nthreads] += ttrial ;
            printf ("parent only  trial: %2d threads: %2d src: %9ld "
                "%10.4f sec\n",
                trial, nthreads, s, ttrial) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (parent, 2)) ;
            GrB_free (&parent) ;

            //------------------------------------------------------------------
            // BFS to compute just level
            //------------------------------------------------------------------

            GrB_free (&level) ;
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, NULL,
                G, s, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tl [nthreads] += ttrial ;

            int32_t maxlevel ;
            GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32, level,
                NULL)) ;

            printf ("level only   trial: %2d threads: %2d src: %9ld "
                "%10.4f sec maxlevel %d\n",
                trial, nthreads, s, ttrial, maxlevel) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (level, 2)) ;

            GrB_free (&level) ;

            //------------------------------------------------------------------
            // BFS to compute both parent and level
            //------------------------------------------------------------------

            GrB_free (&parent) ;
            GrB_free (&level) ;
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, &parent,
                G, s, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tpl [nthreads] += ttrial ;

            GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32, level,
                NULL)) ;
            printf ("parent+level trial: %2d threads: %2d src: %9ld "
                "%10.4f sec maxlevel %d\n",
                trial, nthreads, s, ttrial, maxlevel) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (parent, 2)) ;
            // GrB_TRY (GxB_print (level, 2)) ;

            GrB_free (&parent) ;
            GrB_free (&level) ;
        }

        tp [nthreads] = tp [nthreads] / ntrials ;
        tl [nthreads] = tl [nthreads] / ntrials ;
        tpl [nthreads] = tpl [nthreads] / ntrials ;

        fprintf (stderr, "Avg: BFS parent only  %3d: %10.3f sec: %s\n",
             nthreads, tp [nthreads], matrix_name) ;

        fprintf (stderr, "Avg: BFS level only   %3d: %10.3f sec: %s\n",
             nthreads, tl [nthreads], matrix_name) ;

        fprintf (stderr, "Avg: BFS level+parent %3d: %10.3f sec: %s\n",
             nthreads, tpl [nthreads], matrix_name) ;

    }
    // restore default
    LAGraph_TRY (LAGraph_SetNumThreads (nthreads_max, msg)) ;
    printf ("\n") ;

    GrB_free (&parent) ;
    GrB_free (&level) ;

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;

    return (0) ;
}

