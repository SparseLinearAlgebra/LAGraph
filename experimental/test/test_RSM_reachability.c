#include <GraphBLAS.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "LAGraph.h"

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <acutest.h>
#include <stdint.h>

#define LEN 512
#define MAX_LABELS 3
#define MAX_RESULTS 2000000

#define CHECK(info)                                                                                \
    do                                                                                             \
    {                                                                                              \
        GrB_Info _i = (info);                                                                      \
        if (_i != GrB_SUCCESS && _i != GrB_NO_VALUE)                                               \
        {                                                                                          \
            fprintf(stderr, "GraphBLAS error %d at %s:%d\n", _i, __FILE__, __LINE__);              \
            exit(EXIT_FAILURE);                                                                    \
        }                                                                                          \
    } while (0)

char msg[LAGRAPH_MSG_LEN];

static void setup(void)
{
    LAGraph_Init(msg);
}

static void teardown(void)
{
    LAGraph_Finalize(msg);
}

static GrB_Matrix make_bool_matrix(GrB_Index rows, GrB_Index cols, const GrB_Index *I,
                                   const GrB_Index *J, GrB_Index nvals)
{
    GrB_Matrix A = GrB_NULL;
    TEST_CHECK(GrB_SUCCESS == GrB_Matrix_new(&A, GrB_BOOL, rows, cols));

    for (GrB_Index k = 0; k < nvals; k++)
        TEST_CHECK(GrB_SUCCESS == GrB_Matrix_setElement_BOOL(A, true, I[k], J[k]));
    return A;
}

static LAGraph_Graph make_graph(GrB_Matrix *pA)
{
    LAGraph_Graph G = GrB_NULL;
    int retval = LAGraph_New(&G, pA, LAGraph_ADJACENCY_DIRECTED, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("LAGraph_New failed: %s", msg);
    return G;
}

static void free_graphs(LAGraph_Graph *arr, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        if (arr[i] != GrB_NULL)
            LAGraph_Delete(&arr[i], msg);
}

static void check_result(GrB_Vector result, GrB_Index V, const bool *expected,
                         const char *test_name)
{
    GrB_Index n = 0;
    TEST_CHECK(GrB_SUCCESS == GrB_Vector_size(&n, result));
    TEST_CHECK_(n == V, "%s: result vector size %llu != %llu", test_name, (unsigned long long)n,
                (unsigned long long)V);

    for (GrB_Index v = 0; v < V; v++)
    {
        bool val = false;
        GrB_Info info = GrB_Vector_extractElement_BOOL(&val, result, v);

        if (expected[v])
        {
            TEST_CHECK_(info == GrB_SUCCESS && val,
                        "%s: vertex %llu should be reachable but is not", test_name,
                        (unsigned long long)v);
        }
        else
        {
            bool absent = (info == GrB_NO_VALUE) || (!val);
            TEST_CHECK_(absent, "%s: vertex %llu should not be reachable but is", test_name,
                        (unsigned long long)v);
        }
    }
}

//==============================================================================
// TEST CASES
//==============================================================================

void test_based(void)
{
    setup();

    const GrB_Index Q = 4, V = 6;

    // Graph: 0 -a-> 1 -a-> 2 -b-> 3 -b-> 4 -b-> 5
    GrB_Index Ga_I[] = {0, 1}, Ga_J[] = {1, 2};
    GrB_Index Gb_I[] = {2, 3, 4}, Gb_J[] = {3, 4, 5};

    // RSM: S: q0 -a-> q1 -S-> q2 -b-> q3

    GrB_Index Na_I[] = {0}, Na_J[] = {1};
    GrB_Index Nb_I[] = {1, 2}, Nb_J[] = {3, 3}; // q1->q3, q2->q3

    GrB_Index Call_I[] = {1}, Call_J[] = {0};
    GrB_Index Ret_I[] = {3}, Ret_J[] = {2};
    GrB_Index Ns_I[] = {1}, Ns_J[] = {2}; // N_S: q1->q2

    GrB_Matrix Ga_raw = make_bool_matrix(V, V, Ga_I, Ga_J, 2);
    GrB_Matrix Gb_raw = make_bool_matrix(V, V, Gb_I, Gb_J, 3);
    GrB_Matrix Na_raw = make_bool_matrix(Q, Q, Na_I, Na_J, 1);
    GrB_Matrix Nb_raw = make_bool_matrix(Q, Q, Nb_I, Nb_J, 2);
    GrB_Matrix Cs_raw = make_bool_matrix(Q, Q, Call_I, Call_J, 1);
    GrB_Matrix Rs_raw = make_bool_matrix(Q, Q, Ret_I, Ret_J, 1);
    GrB_Matrix Ns_raw = make_bool_matrix(Q, Q, Ns_I, Ns_J, 1);

    LAGraph_Graph Ga[2], Na[2], Cs[1], Rs[1], Ns[1];
    Ga[0] = make_graph(&Ga_raw);
    Ga[1] = make_graph(&Gb_raw);
    Na[0] = make_graph(&Na_raw);
    Na[1] = make_graph(&Nb_raw);
    Cs[0] = make_graph(&Cs_raw);
    Rs[0] = make_graph(&Rs_raw);
    Ns[0] = make_graph(&Ns_raw);

    GrB_Index QS[] = {0}, QF[] = {3}, Source[] = {0};
    GrB_Vector result = GrB_NULL;

    int retval =
        LAGraph_RSM_reachability(&result, 2, Ga, Na, 1, Ns, Cs, Rs, QS, 1, QF, 1, Source, 1, msg);

    TEST_CHECK(GrB_SUCCESS == retval);
    TEST_MSG("test_based: retval = %d (%s)", retval, msg);

    const bool expected[6] = {false, false, false, false, true, false};
    check_result(result, V, expected, "test_based");

    GrB_free(&result);
    free_graphs(Ga, 2);
    free_graphs(Na, 2);
    free_graphs(Cs, 1);
    free_graphs(Rs, 1);
    free_graphs(Ns, 1);

    teardown();
}

void test_balanced_paren(void)
{
    setup();

    const GrB_Index Q = 5, V = 6;

    // Graph: 0 -a-> 1 -a-> 2 -b-> 3 -b-> 4 -b-> 5
    GrB_Index Ga_I[] = {0, 1, 4}, Ga_J[] = {1, 2, 5};
    GrB_Index Gb_I[] = {2, 3}, Gb_J[] = {3, 4};

    // RSM: S: q0 -(-> q1 -S-> q2 -)-> q3 -S-> q4
    GrB_Index Na_I[] = {0}, Na_J[] = {1};
    GrB_Index Nb_I[] = {1, 1, 2, 2}, Nb_J[] = {3, 4, 3, 4};

    GrB_Index Call_I[] = {1, 3}, Call_J[] = {0, 0};
    GrB_Index Ret_I[] = {4, 4}, Ret_J[] = {2, 4};
    GrB_Index Ns_I[] = {1, 3}, Ns_J[] = {2, 4};

    GrB_Matrix Ga_raw = make_bool_matrix(V, V, Ga_I, Ga_J, 3);
    GrB_Matrix Gb_raw = make_bool_matrix(V, V, Gb_I, Gb_J, 2);
    GrB_Matrix Na_raw = make_bool_matrix(Q, Q, Na_I, Na_J, 1);
    GrB_Matrix Nb_raw = make_bool_matrix(Q, Q, Nb_I, Nb_J, 4);
    GrB_Matrix Cs_raw = make_bool_matrix(Q, Q, Call_I, Call_J, 2);
    GrB_Matrix Rs_raw = make_bool_matrix(Q, Q, Ret_I, Ret_J, 2);
    GrB_Matrix Ns_raw = make_bool_matrix(Q, Q, Ns_I, Ns_J, 2);

    LAGraph_Graph Ga[2], Na[2], Cs[1], Rs[1], Ns[1];
    Ga[0] = make_graph(&Ga_raw);
    Ga[1] = make_graph(&Gb_raw);
    Na[0] = make_graph(&Na_raw);
    Na[1] = make_graph(&Nb_raw);
    Cs[0] = make_graph(&Cs_raw);
    Rs[0] = make_graph(&Rs_raw);
    Ns[0] = make_graph(&Ns_raw);

    GrB_Index QS[] = {0}, QF[] = {4}, Source[] = {0};
    GrB_Vector result = GrB_NULL;

    int retval =
        LAGraph_RSM_reachability(&result, 2, Ga, Na, 1, Ns, Cs, Rs, QS, 1, QF, 1, Source, 1, msg);

    TEST_CHECK(GrB_SUCCESS == retval);
    TEST_MSG("test_based: retval = %d (%s)", retval, msg);

    const bool expected[6] = {false, false, false, false, true, false};
    check_result(result, V, expected, "test_based");

    GrB_free(&result);
    free_graphs(Ga, 2);
    free_graphs(Na, 2);
    free_graphs(Cs, 1);
    free_graphs(Rs, 1);
    free_graphs(Ns, 1);

    teardown();
}

void parens_lrlrlr(void)
{
    setup();

    const GrB_Index Q = 5, V = 7;

    // Graph 0 -(-> 1 -)-> 2 -(-> 3 -)-> 4 -(-> 5 -)-> 6
    GrB_Index Ga_I[] = {0, 2, 4}, Ga_J[] = {1, 3, 5};
    GrB_Index Gb_I[] = {1, 3, 5}, Gb_J[] = {2, 4, 6};

    // RSM: S: q0 -(-> q1 -S-> q2 -)-> q3 -S-> q4
    GrB_Index Na_I[] = {0}, Na_J[] = {1};
    GrB_Index Nb_I[] = {1, 1, 2, 2}, Nb_J[] = {3, 4, 3, 4};

    GrB_Index Call_I[] = {1, 3}, Call_J[] = {0, 0};
    GrB_Index Ret_I[] = {4, 4}, Ret_J[] = {2, 4};
    GrB_Index Ns_I[] = {1, 3}, Ns_J[] = {2, 4};

    GrB_Matrix Ga_raw = make_bool_matrix(V, V, Ga_I, Ga_J, 3);
    GrB_Matrix Gb_raw = make_bool_matrix(V, V, Gb_I, Gb_J, 3);
    GrB_Matrix Na_raw = make_bool_matrix(Q, Q, Na_I, Na_J, 1);
    GrB_Matrix Nb_raw = make_bool_matrix(Q, Q, Nb_I, Nb_J, 4);
    GrB_Matrix Cs_raw = make_bool_matrix(Q, Q, Call_I, Call_J, 2);
    GrB_Matrix Rs_raw = make_bool_matrix(Q, Q, Ret_I, Ret_J, 2);
    GrB_Matrix Ns_raw = make_bool_matrix(Q, Q, Ns_I, Ns_J, 2);

    LAGraph_Graph Ga[2], Na[2], Cs[1], Rs[1], Ns[1];
    Ga[0] = make_graph(&Ga_raw);
    Ga[1] = make_graph(&Gb_raw);
    Na[0] = make_graph(&Na_raw);
    Na[1] = make_graph(&Nb_raw);
    Cs[0] = make_graph(&Cs_raw);
    Rs[0] = make_graph(&Rs_raw);
    Ns[0] = make_graph(&Ns_raw);

    GrB_Index QS[] = {0}, QF[] = {4}, Source[] = {0};
    GrB_Vector result = GrB_NULL;

    int retval =
        LAGraph_RSM_reachability(&result, 2, Ga, Na, 1, Ns, Cs, Rs, QS, 1, QF, 1, Source, 1, msg);

    TEST_CHECK(GrB_SUCCESS == retval);
    TEST_MSG("test_based: retval = %d (%s)", retval, msg);

    const bool expected[7] = {0, 0, 1, 0, 1, 0, 1};
    check_result(result, V, expected, "test_based");

    GrB_free(&result);
    free_graphs(Ga, 2);
    free_graphs(Na, 2);
    free_graphs(Cs, 1);
    free_graphs(Rs, 1);
    free_graphs(Ns, 1);

    teardown();
}

void very_big_parens(void)
{
    setup();

    const GrB_Index Q = 5, V = 13;

    // Graph (()(()()))()
    GrB_Index Ga_I[] = {0, 1, 3, 4, 6, 10}, Ga_J[] = {1, 2, 4, 5, 7, 11};
    GrB_Index Gb_I[] = {2, 5, 7, 8, 9, 11}, Gb_J[] = {3, 6, 8, 9, 10, 12};

    // RSM: S: q0 -(-> q1 -S-> q2 -)-> q3 -S-> q4
    GrB_Index Na_I[] = {0}, Na_J[] = {1};
    GrB_Index Nb_I[] = {1, 1, 2, 2}, Nb_J[] = {3, 4, 3, 4};

    GrB_Index Call_I[] = {1, 3}, Call_J[] = {0, 0};
    GrB_Index Ret_I[] = {4, 4}, Ret_J[] = {2, 4};
    GrB_Index Ns_I[] = {1, 3}, Ns_J[] = {2, 4};

    GrB_Matrix Ga_raw = make_bool_matrix(V, V, Ga_I, Ga_J, 6);
    GrB_Matrix Gb_raw = make_bool_matrix(V, V, Gb_I, Gb_J, 6);
    GrB_Matrix Na_raw = make_bool_matrix(Q, Q, Na_I, Na_J, 1);
    GrB_Matrix Nb_raw = make_bool_matrix(Q, Q, Nb_I, Nb_J, 4);
    GrB_Matrix Cs_raw = make_bool_matrix(Q, Q, Call_I, Call_J, 2);
    GrB_Matrix Rs_raw = make_bool_matrix(Q, Q, Ret_I, Ret_J, 2);
    GrB_Matrix Ns_raw = make_bool_matrix(Q, Q, Ns_I, Ns_J, 2);

    LAGraph_Graph Ga[2], Na[2], Cs[1], Rs[1], Ns[1];
    Ga[0] = make_graph(&Ga_raw);
    Ga[1] = make_graph(&Gb_raw);
    Na[0] = make_graph(&Na_raw);
    Na[1] = make_graph(&Nb_raw);
    Cs[0] = make_graph(&Cs_raw);
    Rs[0] = make_graph(&Rs_raw);
    Ns[0] = make_graph(&Ns_raw);

    GrB_Index QS[] = {0}, QF[] = {4}, Source[] = {0};
    GrB_Vector result = GrB_NULL;

    int retval =
        LAGraph_RSM_reachability(&result, 2, Ga, Na, 1, Ns, Cs, Rs, QS, 1, QF, 1, Source, 1, msg);

    TEST_CHECK(GrB_SUCCESS == retval);
    TEST_MSG("test_based: retval = %d (%s)", retval, msg);

    const bool expected[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1};
    check_result(result, V, expected, "test_based");

    GrB_free(&result);
    free_graphs(Ga, 2);
    free_graphs(Na, 2);
    free_graphs(Cs, 1);
    free_graphs(Rs, 1);
    free_graphs(Ns, 1);

    teardown();
}

TEST_LIST = {{"simple", test_based},
             {"test_balanced_paren", test_balanced_paren},
             {"parens_lrlrlr", parens_lrlrlr},
             {"very_big_parens", very_big_parens},
             {NULL, NULL}};
