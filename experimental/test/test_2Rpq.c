#include <stdio.h>
#include <acutest.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

#define LEN 512
#define MAX_LABELS 3
#define MAX_RESULTS 2000000

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G[MAX_LABELS] ;
LAGraph_Graph R[MAX_LABELS] ;
GrB_Matrix A ;

char testcase_name [LEN+1] ;
char filename [LEN+1] ;

typedef struct
{
    const char* name ;
    const char* graphs[MAX_LABELS] ;
    const char* fas[MAX_LABELS] ;
    const char* fa_meta ;
    const char* sources ;
    const GrB_Index expected[MAX_RESULTS] ;
    const size_t expected_count ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {"simple 1 or more",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/1_a.mtx", NULL },                    // Regex: a+
     "rpq_data/1_meta.txt",
     "rpq_data/1_sources.txt",
     {2, 4, 6, 7}, 4},
    {"simple kleene star",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/2_a.mtx", "rpq_data/2_b.mtx", NULL}, // Regex: (a b)*
     "rpq_data/2_meta.txt",
     "rpq_data/2_sources.txt",
     {2, 6, 8}, 3},
    {"kleene star of the conjunction",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"rpq_data/3_a.mtx", "rpq_data/3_b.mtx", NULL}, // Regex: (a | b)*
     "rpq_data/3_meta.txt",
     "rpq_data/3_sources.txt",
     {3, 6}, 2},
    {"simple repeat from n to m times",
     {"rpq_data/a.mtx",   "rpq_data/b.mtx",   NULL},
     {"",                 "rpq_data/4_b.mtx", NULL}, // Regex: b b b (b b)?
     "rpq_data/4_meta.txt",
     "rpq_data/4_sources.txt",
     {3, 4, 6}, 3},
    {NULL, NULL, NULL, NULL},
} ;

//****************************************************************************
static void load_testcase
(
    const char* const graphs [ ],
    const char* const fas [ ],
    const char *fa_meta,
    const char *sources,
    GrB_Index S [ ],
    size_t *ns,
    GrB_Index QS [ ],
    size_t *nqs,
    GrB_Index QF [ ],
    size_t *nqf
)
{
    // Load graph from MTX files representing its adjacency matrix
    // decomposition
    for (int i = 0 ; ; i++)
    {
        const char *name = graphs[i] ;

        if (name == NULL) break ;
        if (strlen(name) == 0) continue ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f));

        OK (LAGraph_New (&(G[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

        TEST_CHECK (A == NULL) ;
    }

    // Load NFA from MTX files representing its adjacency matrix
    // decomposition
    for (int i = 0 ; ; i++)
    {
        const char *name = fas[i] ;

        if (name == NULL) break ;
        if (strlen(name) == 0) continue ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", name) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;

        OK (LAGraph_New (&(R[i]), &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        OK (LAGraph_Cached_AT (R[i], msg)) ;

        TEST_CHECK (A == NULL) ;
    }

    // Note the matrix rows/cols are enumerated from 0 to n-1. Meanwhile, in
    // MTX format they are enumerated from 1 to n. Thus, when
    // loading/comparing the results these values should be
    // decremented/incremented correspondingly.

    // Load graph source nodes from the sources file
    GrB_Index s ;
    *ns = 0 ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", sources) ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;

    while (fscanf(f, "%ld", &s) != EOF)
        S[(*ns)++] = s - 1 ;

    OK (fclose(f)) ;

    // Load NFA starting states from the meta file
    GrB_Index qs ;
    *nqs = 0 ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", fa_meta) ;
    f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;

    TEST_CHECK (fscanf(f, "%ld", nqs) != EOF) ;

    for (uint64_t i = 0; i < *nqs; i++) {
        TEST_CHECK (fscanf(f, "%ld", &qs) != EOF) ;
        QS[i] = qs - 1 ;
    }

    // Load NFA final states from the same file
    uint64_t qf ;
    *nqf = 0 ;

    TEST_CHECK (fscanf(f, "%ld", nqf) != EOF) ;

    for (uint64_t i = 0; i < *nqf; i++) {
        TEST_CHECK (fscanf(f, "%ld", &qf) != EOF) ;
        QF[i] = qf - 1 ;
    }

    OK (fclose(f)) ;
}

//****************************************************************************
static void free_testcase (void)
{
    for (uint64_t i = 0 ; i < MAX_LABELS ; i++)
    {
        if (G[i] == NULL) continue ;
        OK (LAGraph_Delete (&(G[i]), msg)) ;
    }

    for (uint64_t i = 0 ; i < MAX_LABELS ; i++ )
    {
        if (R[i] == NULL) continue ;
        OK (LAGraph_Delete (&(R[i]), msg)) ;
    }
}

//****************************************************************************
void test_Rpq_Simple (void)
{
    LAGraph_Init (msg) ;
    LAGraph_Rpq_initialize (msg) ;

    for (int k = 0 ; ; k++)
    {
        if (files[k].sources == NULL) break ;

        GrB_Index S[16] ;
        size_t ns ;
        GrB_Index QS[16] ;
        size_t nqs ;
        GrB_Index QF[16] ;
        size_t nqf ;

        load_testcase (files[k].graphs, files[k].fas, files[k].fa_meta,
            files[k].sources, S, &ns, QS, &nqs, QF, &nqf) ;

        bool inverse_labels[] = {false, false, false, false, false, false, false, false, false, false, false, false, false};
        bool inverse = false;

        // The operators carry JIT definitions, so run every semantics with
        // both generic and JIT kernels.
        for (int jit = 0 ; jit <= 1 ; jit++)
        {
            snprintf (testcase_name, LEN,
                "basic regular path query %s (jit: %d)", files[k].name, jit) ;
            TEST_CASE (testcase_name) ;
            OK (LG_SET_JIT (jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;
            printf ("jit: %d\n", jit) ;

            Path *paths ;
            size_t path_count ;

            OK (LAGraph_2Rpq_AllSimple (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, msg)) ;

            printf("ALL SIMPLE:\n");
            for (size_t i = 0 ; i < path_count ; i++)
            {
                Path_print (&paths[i]);
            }
            printf("\n");

            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            OK (LAGraph_2Rpq_AllTrails (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, msg)) ;

            printf("ALL TRAILS:\n");
            for (size_t i = 0 ; i < path_count ; i++)
            {
                Path_print (&paths[i]);
            }
            printf("\n");

            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            OK (LAGraph_2Rpq_AllPaths (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, 10, msg)) ;

            printf("ALL PATHS (LIMIT = 10):\n");
            for (size_t i = 0 ; i < path_count ; i++)
            {
                Path_print (&paths[i]);
            }
            printf("\n");

            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            // All shortest paths are searched from exactly one source vertex.
            printf("ALL SHORTEST PATHS:\n");
            for (size_t s = 0 ; s < ns ; s++)
            {
                OK (LAGraph_2Rpq_AllShortestPaths (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, &S[s], 1, inverse, 100,
                                        msg)) ;

                for (size_t i = 0 ; i < path_count ; i++)
                {
                    Path_print (&paths[i]);
                }

                OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;
            }
            printf("\n");
        }

        free_testcase () ;
    }

    OK (LG_SET_JIT (GxB_JIT_ON)) ;
    LAGraph_Finalize (msg) ;
}


//****************************************************************************
#define MAX_PATH_VERTICES 4
#define MAX_EXPECTED_PATHS 4

typedef struct
{
    const size_t vertex_count ;
    const GrB_Index vertices[MAX_PATH_VERTICES] ;
    const GrB_Index labels[MAX_PATH_VERTICES] ;
}
path_info ;

typedef struct
{
    const char* name ;
    const char* graphs[MAX_LABELS] ;
    const char* fas[MAX_LABELS] ;
    const char* fa_meta ;
    const char* sources ;
    const path_info simple[MAX_EXPECTED_PATHS] ;
    const size_t simple_count ;
    const path_info trails[MAX_EXPECTED_PATHS] ;
    const size_t trails_count ;
    const path_info all[MAX_EXPECTED_PATHS] ;
    const size_t all_count ;
    const path_info shortest[MAX_EXPECTED_PATHS] ;
    const size_t shortest_count ;
}
labelled_info ;

// The graph joins vertices 1 and 2 by both labels, so a path is determined by
// its labels and not by its vertices alone. Expected vertices and labels are
// written in the MTX numbering, starting from 1.
const labelled_info labelled_files [ ] =
{
    {"parallel edges with distinct labels",
     {"rpq_data/labels_a.mtx", "rpq_data/labels_b.mtx", NULL},
     {"rpq_data/5_a.mtx",      "rpq_data/5_b.mtx",      NULL}, // Regex: a | b
     "rpq_data/5_meta.txt",
     "rpq_data/5_sources.txt",
     {{2, {1, 2}, {1}}, {2, {1, 2}, {2}}}, 2,
     {{2, {1, 2}, {1}}, {2, {1, 2}, {2}}}, 2,
     {{2, {1, 2}, {1}}, {2, {1, 2}, {2}}}, 2,
     {{2, {1, 2}, {1}}, {2, {1, 2}, {2}}}, 2},
    {"trail revisiting a vertex pair under another label",
     {"rpq_data/labels_a.mtx", "rpq_data/labels_b.mtx", NULL},
     {"rpq_data/6_a.mtx",      "rpq_data/6_b.mtx",      NULL}, // Regex: a a b
     "rpq_data/6_meta.txt",
     "rpq_data/6_sources.txt",
     {{0}}, 0,
     {{4, {1, 2, 1, 2}, {1, 1, 2}}}, 1,
     {{4, {1, 2, 1, 2}, {1, 1, 2}}}, 1,
     {{4, {1, 2, 1, 2}, {1, 1, 2}}}, 1},
    {NULL, NULL, NULL, NULL},
} ;

static bool path_matches (const Path *path, const path_info *expected)
{
    if (path->vertex_count != expected->vertex_count) return false ;

    for (size_t i = 0 ; i < path->vertex_count ; i++)
    {
        if (path->vertices[2 * i] + 1 != expected->vertices[i]) return false ;
    }

    for (size_t i = 0 ; i + 1 < path->vertex_count ; i++)
    {
        if (path->vertices[2 * i + 1] + 1 != expected->labels[i]) return false ;
    }

    return true ;
}

static void check_paths
(
    const char *semantics,
    const Path *paths,
    size_t path_count,
    const path_info *expected,
    size_t expected_count
)
{
    TEST_CHECK (path_count == expected_count) ;
    TEST_MSG ("%s: got %zu paths, expected %zu", semantics, path_count,
        expected_count) ;

    for (size_t i = 0 ; i < expected_count ; i++)
    {
        bool found = false ;

        for (size_t j = 0 ; j < path_count ; j++)
        {
            if (path_matches (&paths[j], &expected[i])) found = true ;
        }

        TEST_CHECK (found) ;
        TEST_MSG ("%s: expected path %zu is missing", semantics, i) ;
    }

    for (size_t j = 0 ; j < path_count ; j++)
    {
        bool found = false ;

        for (size_t i = 0 ; i < expected_count ; i++)
        {
            if (path_matches (&paths[j], &expected[i])) found = true ;
        }

        TEST_CHECK (found) ;
        TEST_MSG ("%s: unexpected path %zu", semantics, j) ;
        if (!found) Path_print (&paths[j]) ;
    }
}

//****************************************************************************
void test_Rpq_Labels (void)
{
    LAGraph_Init (msg) ;
    LAGraph_Rpq_initialize (msg) ;

    for (int k = 0 ; ; k++)
    {
        if (labelled_files[k].sources == NULL) break ;

        GrB_Index S[16] ;
        size_t ns ;
        GrB_Index QS[16] ;
        size_t nqs ;
        GrB_Index QF[16] ;
        size_t nqf ;

        load_testcase (labelled_files[k].graphs, labelled_files[k].fas,
            labelled_files[k].fa_meta, labelled_files[k].sources, S, &ns,
            QS, &nqs, QF, &nqf) ;

        bool inverse_labels[MAX_LABELS] = {false, false, false} ;
        bool inverse = false ;

        for (int jit = 0 ; jit <= 1 ; jit++)
        {
            snprintf (testcase_name, LEN,
                "labelled regular path query %s (jit: %d)",
                labelled_files[k].name, jit) ;
            TEST_CASE (testcase_name) ;
            OK (LG_SET_JIT (jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;

            Path *paths ;
            size_t path_count ;

            OK (LAGraph_2Rpq_AllSimple (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, msg)) ;
            check_paths ("ALL SIMPLE", paths, path_count,
                labelled_files[k].simple, labelled_files[k].simple_count) ;
            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            OK (LAGraph_2Rpq_AllTrails (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, msg)) ;
            check_paths ("ALL TRAILS", paths, path_count,
                labelled_files[k].trails, labelled_files[k].trails_count) ;
            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            OK (LAGraph_2Rpq_AllPaths (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, 10, msg)) ;
            check_paths ("ALL PATHS", paths, path_count,
                labelled_files[k].all, labelled_files[k].all_count) ;
            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;

            OK (LAGraph_2Rpq_AllShortestPaths (&paths, &path_count, R,
                                        inverse_labels, MAX_LABELS, QS, nqs,
                                        QF, nqf, G, S, ns, inverse, 100,
                                        msg)) ;
            check_paths ("ALL SHORTEST PATHS", paths, path_count,
                labelled_files[k].shortest, labelled_files[k].shortest_count) ;
            OK (LAGraph_2Rpq_FreePaths (&paths, path_count, msg)) ;
        }

        free_testcase () ;
    }

    OK (LG_SET_JIT (GxB_JIT_ON)) ;
    LAGraph_Finalize (msg) ;
}

TEST_LIST = {
    {"Rpq_Simple", test_Rpq_Simple},
    {"Rpq_Labels", test_Rpq_Labels},
    {NULL, NULL}
};
