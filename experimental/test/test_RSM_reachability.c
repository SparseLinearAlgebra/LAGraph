#include "GraphBLAS.h"
#include "LAGraph.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

#define LEN 512
#define MAX_LABELS 3
#define MAX_RESULTS 2000000

LAGraph_Graph G[MAX_LABELS] ;
LAGraph_Graph R[MAX_LABELS] ;
LAGraph_Graph R_call;

GrB_Matrix A ;

char testcase_name[LEN+1] ;
char filename[LEN+1] ;
char msg[LAGRAPH_MSG_LEN] ;

typedef struct
{
    const char* name ;
    const char* graphs[MAX_LABELS] ;
    const char* rsm[MAX_LABELS] ;
    const char* rsm_call;
    const char* rsm_meta ;
    const char* sources ;
    const GrB_Index expected[MAX_RESULTS] ;
    const size_t expected_count ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {"simple 1",
     {"rsm_data/a.mtx", "rsm_data/b.mtx", NULL},
     {"rsm_data/1_a.mtx", NULL},
     "rsm_data/1_call.mtx",
     "rsm_data/1_meta.txt",
     "rsm_data/1_sources.txt",
     {2},
     1},
    {NULL, NULL, NULL, NULL},
} ;


void test_RSM_reachability (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; ++k) 
    {
        if (files[k].sources == NULL) break ;

        snprintf (testcase_name, LEN, "basic context-free path query %s", files[k].name) ;
        TEST_CASE (testcase_name) ;

        for (int check_symmetry = 0 ; check_symmetry < 2  ; ++check_symmetry)
        {
            // Load graph from MTX files represention its adjacency matrix
            // decomposition
            for (int i = 0 ; ; ++i)
            {
                const char *name = files[k].graphs[i] ;

                if (name == NULL) break ;
                if (strlen(name) == 0) continue ;

                snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
                FILE *f = fopen (filename, "r") ;
                TEST_CHECK (f != NULL) ;
                OK (LAGraph_MMRead (&A, f, msg)) ;
                OK (fclose (f)) ;

                OK (LAGraph_New (&(G[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

                TEST_CHECK (A == NULL) ;
            }

            // Load RSM from MTX files representing its adjacency matrix
            // decomposition
            for (int i = 0 ; ; ++i)
            {
                const char *name = files[k].rsm[i] ;

                if (name == NULL) break ;
                if (strlen(name) == 0) continue ;

                snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
                FILE *f = fopen (filename, "r") ;
                TEST_CHECK (f != NULL) ;
                OK (LAGraph_MMRead (&A, f, msg)) ;
                OK (fclose (f)) ;

                OK (LAGraph_New (&(R[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

                if (check_symmetry)
                {
                    // Check if the pattern is symmetric - if it isn't make it.
                    // Note this also computes R[i]->AT
                    OK (LAGraph_Cached_IsSymmetricStructure (R[i], msg)) ;
                }

                TEST_CHECK (A == NULL) ;
            }

            {
                const char *name = files[k].rsm_call ;
            
                snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
                FILE *f = fopen (filename, "r") ;
                TEST_CHECK (f != NULL) ;
                OK (LAGraph_MMRead (&A, f, msg)) ;
                OK (fclose (f)) ;

                OK (LAGraph_New (&(R_call), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

                if (check_symmetry)
                {
                    OK (LAGraph_Cached_IsSymmetricStructure(R_call, msg)) ;
                }

                TEST_CHECK (A == NULL) ; // зачем нам это проверять?, почему оно NULL
            }

        }

        // Note the matrix rows/cols are enumerated from 0 to n-1.
        // Meanwhile, in MTX format they are enumerated from 1 to n. Thus,
        // when loading/comparing the results these values should be
        // decremented/incremented correspondingly.

        // Load graph source nodes from the sources file
        GrB_Index s ;
        GrB_Index S[16] ;
        size_t ns = 0 ;

        const char *name = files[k].sources ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;

        while (fscanf(f, "%" PRIu64, &s) != EOF)
        {
            S[ns++] = s - 1 ;
        }

        OK (fclose(f)) ;

        // Load RSM starting states from the meta file
        GrB_Index qs ;
        GrB_Index QS[16] ; // почему именно 16, мб личный выбор
        size_t nqs = 0 ;

        name = files[k].rsm_meta ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ; ;

        uint64_t nqs64 = 0 ;
        TEST_CHECK (fscanf (f, "%" PRIu64, &nqs64) != EOF) ;
        nqs = (size_t) nqs64 ;

        for (size_t i = 0 ; i < nqs ; ++i)
        {
            TEST_CHECK(fscanf (f, "%" PRIu64, &qs) != EOF) ;
            QS[i] = qs - 1 ;
        }

        // Load RSM final states from the same file
        uint64_t qf ;
        uint64_t QF[16] ;
        size_t nqf = 0 ;
        uint64_t nqf64 = 0 ;

        TEST_CHECK (fscanf (f, "%" PRIu64, &nqf64) != EOF) ;
        nqf = (size_t) nqf64 ;

        for (size_t i = 0 ; i < nqf ; ++i)
        {
            TEST_CHECK (fscanf (f, "%" PRIu64, &qf) != EOF) ;
            QF[i] = qf - 1 ;
        }

        OK (fclose (f)) ;

        // Evaluate the algorithm
        GrB_Vector r = NULL ;

        OK (LAGraph_RSM_reachability (&r, R, MAX_LABELS, QS, nqs,
                                        QF, nqf, &R_call, G, S, ns, msg)) ;
        
        // Extract results from the output vector
        GrB_Index *reachable ;
        bool *values ;
        
        GrB_Index nvals ;
        GrB_Vector_nvals (&nvals, r) ;

        OK (LAGraph_Malloc ((void **) &reachable, MAX_RESULTS, sizeof (GrB_Index), msg)) ;;
        OK (LAGraph_Malloc ((void **) &values, MAX_RESULTS, sizeof (bool), msg)) ;

        GrB_Vector_extractTuples (reachable, values, &nvals, r) ;
        
        TEST_MSG("returned %lu values:\n", nvals);
        for (uint64_t i = 0; i < nvals; i++)
            TEST_MSG("  %lu\n", reachable[i] + 1);

        // Compare the results with expected values
        TEST_CHECK (nvals == files[k].expected_count) ;
        for (uint64_t i = 0 ; i < nvals ; ++i)
            TEST_CHECK (reachable[i] + 1 == files[k].expected[i]) ;

        // Cleanup
        OK (LAGraph_Free ((void **) &values, NULL)) ;
        OK (LAGraph_Free ((void **) &reachable, NULL)) ;

        OK (GrB_free (&r)) ;

        for (uint64_t i = 0 ; i < MAX_LABELS ; ++i)
        {
            if (G[i] == NULL) continue ;
            OK (LAGraph_Delete (&(G[i]), msg)) ;
        }

        for (uint64_t i = 0; i < MAX_LABELS ; ++i)
        {
            if (R[i] == NULL) continue ;
            OK (LAGraph_Delete (&(R[i]), msg)) ;
        }
    }

    LAGraph_Finalize (msg) ;
}

TEST_LIST = {
    {"RSM_reachability", test_RSM_reachability},
    {NULL, NULL}
} ;
