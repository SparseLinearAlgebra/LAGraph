#include "LAGraph.h"
#include "LAGraphX.h"
#include "LAGraph_test.h"
#include <GraphBLAS.h>

char msg[LAGRAPH_MSG_LEN];

//==============================================================================
// test setup / teardown
//==============================================================================

static void setup(void)
{
    OK(LAGraph_Init(msg));
}

static void teardown(void)
{
    OK(LAGraph_Finalize(msg));
}

static GrB_Info build_call_matrix(GrB_Matrix *call, GrB_Matrix rsm_nt, GrB_Vector rsm_start,
                                  GrB_Index Q)
{
    GrB_Vector q_call_mask = GrB_NULL;

    GRB_TRY(GrB_Matrix_new(call, GrB_BOOL, Q, Q));
    GRB_TRY(GrB_Vector_new(&q_call_mask, GrB_BOOL, Q));

    // [q_call, q_ret] -> [q_call]
    GRB_TRY(GrB_reduce(q_call_mask, GrB_NULL, GrB_NULL, GrB_LOR_MONOID_BOOL, rsm_nt, GrB_NULL));
    // |Q|*1 * 1*|Q| -> |Q|*|Q|
    GRB_TRY(GrB_mxm(*call, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, (GrB_Matrix)q_call_mask,
                    (GrB_Matrix)rsm_start, GrB_DESC_T1));

    GRB_TRY(GrB_Vector_free(&q_call_mask));

    return GrB_SUCCESS;
}

static void check_reachable(GrB_Vector result, GrB_Index *expected, GrB_Index n_expected,
                            GrB_Index V, const char *test_name)
{
    GrB_Index nvals = 0;
    OK(GrB_Vector_nvals(&nvals, result));

    bool *got = (bool *)calloc(V, sizeof(bool));
    TEST_CHECK(got != NULL);

    GrB_Index *idx = GrB_NULL;
    bool *val = GrB_NULL;

    if (nvals > 0)
    {
        idx = (GrB_Index *)malloc(nvals * sizeof(GrB_Index));
        val = (bool *)malloc(nvals * sizeof(bool));
        if (!TEST_CHECK(idx != NULL && val != NULL))
        {
            free(got);
            free(idx);
            free(val);
            return;
        }
        OK(GrB_Vector_extractTuples_BOOL(idx, val, &nvals, result));
        for (GrB_Index k = 0; k < nvals; ++k)
            if (val[k])
                got[idx[k]] = true;
    }

    // every expected vertex must be reachable
    for (GrB_Index k = 0; k < n_expected; ++k)
    {
        TEST_CHECK(got[expected[k]]);
        TEST_MSG("%s: vertex %llu should be reachable", test_name, (unsigned long long)expected[k]);
    }

    GrB_Index actual_count = 0;
    for (GrB_Index v = 0; v < V; ++v)
        if (got[v])
            actual_count++;

    TEST_CHECK(actual_count == n_expected);
    TEST_MSG("%s: expected %llu reachable vertices, got %llu", test_name,
             (unsigned long long)n_expected, (unsigned long long)actual_count);

    free(got);
    free(idx);
    free(val);
}

//==============================================================================
// RSM builders
//==============================================================================
// Grammar 1:  S -> a S b | a b
static void rsm_aSb_ab(GrB_Matrix *term, GrB_Matrix *nt, GrB_Matrix *call, GrB_Vector *start,
                       GrB_Vector *final)
{
    GrB_Index Q = 4;
    OK(GrB_Matrix_new(&term[0], GrB_BOOL, Q, Q));
    OK(GrB_Matrix_new(&term[1], GrB_BOOL, Q, Q));
    OK(GrB_Matrix_new(&nt[0], GrB_BOOL, Q, Q));
    OK(GrB_Vector_new(&start[0], GrB_BOOL, Q));
    OK(GrB_Vector_new(&final[0], GrB_BOOL, Q));

    OK(GrB_Matrix_setElement_BOOL(term[0], true, 0, 1)); // 0 -a-> 1
    OK(GrB_Matrix_setElement_BOOL(term[1], true, 1, 3)); // 1 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(term[1], true, 2, 3)); // 2 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(nt[0], true, 1, 2));   // 1 -S-> 2
    OK(GrB_Vector_setElement_BOOL(start[0], true, 0));
    OK(GrB_Vector_setElement_BOOL(final[0], true, 3));

    OK(build_call_matrix(call, *nt, *start, Q));
}

// Grammar 2:  S -> a S b | c
static void rsm_aSb_c(GrB_Matrix *term, GrB_Matrix *nt, GrB_Matrix *call, GrB_Vector *start,
                      GrB_Vector *final)
{
    GrB_Index Q = 4;
    OK(GrB_Matrix_new(&term[0], GrB_BOOL, Q, Q));
    OK(GrB_Matrix_new(&term[1], GrB_BOOL, Q, Q));
    OK(GrB_Matrix_new(&term[2], GrB_BOOL, Q, Q));
    OK(GrB_Matrix_new(&nt[0], GrB_BOOL, Q, Q));
    OK(GrB_Vector_new(&start[0], GrB_BOOL, Q));
    OK(GrB_Vector_new(&final[0], GrB_BOOL, Q));

    OK(GrB_Matrix_setElement_BOOL(term[0], true, 0, 1)); // 0 -a-> 1
    OK(GrB_Matrix_setElement_BOOL(term[2], true, 0, 3)); // 1 -c-> 3
    OK(GrB_Matrix_setElement_BOOL(term[1], true, 2, 3)); // 2 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(nt[0], true, 1, 2));   // 1 -S-> 2
    OK(GrB_Vector_setElement_BOOL(start[0], true, 0));
    OK(GrB_Vector_setElement_BOOL(final[0], true, 3));

    OK(build_call_matrix(call, *nt, *start, Q));
}

// Grammar 3: S -> X | Y | S X | S Y
//            X -> a S b | a b
//            Y -> c S d | c d
static void rsm_complex(GrB_Matrix *term, GrB_Matrix *nt, GrB_Matrix *call, GrB_Vector *start,
                        GrB_Vector *final)
{
    GrB_Index Q = 11;
    for (int t = 0; t < 4; ++t)
        OK(GrB_Matrix_new(&term[t], GrB_BOOL, Q, Q));
    for (int n = 0; n < 3; ++n)
    {
        OK(GrB_Matrix_new(&nt[n], GrB_BOOL, Q, Q));
        OK(GrB_Vector_new(&start[n], GrB_BOOL, Q));
        OK(GrB_Vector_new(&final[n], GrB_BOOL, Q));
    }

    // --- S Component (States 0-2) ---
    // S -> X | Y | S X | S Y
    OK(GrB_Matrix_setElement_BOOL(nt[1], true, 0, 2)); // 0 -X-> 2
    OK(GrB_Matrix_setElement_BOOL(nt[2], true, 0, 2)); // 0 -Y-> 2
    OK(GrB_Matrix_setElement_BOOL(nt[0], true, 0, 1)); // 0 -S-> 1
    OK(GrB_Matrix_setElement_BOOL(nt[1], true, 1, 2)); // 1 -X-> 2
    OK(GrB_Matrix_setElement_BOOL(nt[2], true, 1, 2)); // 1 -Y-> 2
    OK(GrB_Vector_setElement_BOOL(start[0], true, 0));
    OK(GrB_Vector_setElement_BOOL(final[0], true, 2));

    // --- X Component (States 3-6) ---
    OK(GrB_Matrix_setElement_BOOL(term[0], true, 3, 4)); // 3 -a-> 4
    OK(GrB_Matrix_setElement_BOOL(term[1], true, 4, 6)); // 4 -b-> 6
    OK(GrB_Matrix_setElement_BOOL(term[1], true, 5, 6)); // 5 -b-> 6
    OK(GrB_Matrix_setElement_BOOL(nt[0], true, 4, 5));   // 4 -S-> 5
    OK(GrB_Vector_setElement_BOOL(start[1], true, 3));
    OK(GrB_Vector_setElement_BOOL(final[1], true, 6));

    // --- Y Component (States 7-10)
    OK(GrB_Matrix_setElement_BOOL(term[2], true, 7, 8));  // 7 -c-> 8
    OK(GrB_Matrix_setElement_BOOL(term[3], true, 8, 10)); // 8 -d-> 10
    OK(GrB_Matrix_setElement_BOOL(term[3], true, 9, 10)); // 9 -d-> 10
    OK(GrB_Matrix_setElement_BOOL(nt[0], true, 8, 9));    // 8 -S-> 9
    OK(GrB_Vector_setElement_BOOL(start[2], true, 7));
    OK(GrB_Vector_setElement_BOOL(final[2], true, 10));

    for (int n = 0; n < 3; n++)
        OK(build_call_matrix(&call[n], nt[n], start[n], Q));
}

//==============================================================================
// test cases
//==============================================================================

//------------------------------------------------------------------------------
// TC1: S -> a S b | a b
// Graph: 0-a->1; 1-a->2; 2-a->0; 0-b->3; 3-b->0
// Source: 1   Expected: {0, 3}
//------------------------------------------------------------------------------
void test_TC1_cyclic_ab(void)
{
    setup();

    GrB_Matrix term[2], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_ab(term, nt, call, start, final);

    GrB_Index V = 4;
    GrB_Matrix gterm[2];
    gterm[0] = gterm[1] = GrB_NULL;
    OK(GrB_Matrix_new(&gterm[0], GrB_BOOL, V, V));
    OK(GrB_Matrix_new(&gterm[1], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 1, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 2, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 0, 3));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 3, 0));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 2, term, gterm, 1, nt, call, start, final, 0, 1, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {0, 3};
        check_reachable(result, expected, 2, V, "TC1");
    }

    GrB_free(&result);
    for (int i = 0; i < 2; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

//------------------------------------------------------------------------------
// TC2: S -> a S b | a b
// Graph: 0-a->0; 0-a->1; 1-b->1
// Source: 0   Expected: {1}
//------------------------------------------------------------------------------
void test_TC2_self_loops(void)
{
    setup();

    GrB_Matrix term[2], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_ab(term, nt, call, start, final);

    GrB_Index V = 2;
    GrB_Matrix gterm[2];
    gterm[0] = gterm[1] = GrB_NULL;
    OK(GrB_Matrix_new(&gterm[0], GrB_BOOL, V, V));
    OK(GrB_Matrix_new(&gterm[1], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 1, 1));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 2, term, gterm, 1, nt, call, start, final, 0, 0, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {1};
        check_reachable(result, expected, 1, V, "TC2");
    }

    GrB_free(&result);
    for (int i = 0; i < 2; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

//------------------------------------------------------------------------------
// TC3: S -> a S b | a b
// Graph: 0-a->0; 0-b->0
// Source: 0   Expected: {0}
//------------------------------------------------------------------------------
void test_TC3_single_node(void)
{
    setup();

    GrB_Matrix term[2], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_ab(term, nt, call, start, final);

    GrB_Index V = 1;
    GrB_Matrix gterm[2];
    gterm[0] = gterm[1] = GrB_NULL;
    OK(GrB_Matrix_new(&gterm[0], GrB_BOOL, V, V));
    OK(GrB_Matrix_new(&gterm[1], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 0, 0));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 2, term, gterm, 1, nt, call, start, final, 0, 0, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {0};
        check_reachable(result, expected, 1, V, "TC3");
    }

    GrB_free(&result);
    for (int i = 0; i < 2; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

//------------------------------------------------------------------------------
// TC4: S -> a S b | c
// Graph: 0-a->0; 0-c->1; 1-b->1
// Source: 0   Expected: {1}
//------------------------------------------------------------------------------
void test_TC4_aSb_c_loops(void)
{
    setup();

    GrB_Matrix term[3], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = term[2] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_c(term, nt, call, start, final);

    GrB_Index V = 2;
    GrB_Matrix gterm[3];
    gterm[0] = gterm[1] = gterm[2] = GrB_NULL;
    for (int i = 0; i < 3; ++i)
        OK(GrB_Matrix_new(&gterm[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 1, 1));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 3, term, gterm, 1, nt, call, start, final, 0, 0, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {1};
        check_reachable(result, expected, 1, V, "TC4");
    }

    GrB_free(&result);
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

//------------------------------------------------------------------------------
// TC5: S -> a S b | c
// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->2
// Source: 1   Expected: {3}
//------------------------------------------------------------------------------
void test_TC5_aSb_c_cycle1(void)
{
    setup();

    GrB_Matrix term[3], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = term[2] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_c(term, nt, call, start, final);

    GrB_Index V = 4;
    GrB_Matrix gterm[3];
    gterm[0] = gterm[1] = gterm[2] = GrB_NULL;
    for (int i = 0; i < 3; ++i)
        OK(GrB_Matrix_new(&gterm[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 3, 2));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 3, term, gterm, 1, nt, call, start, final, 0, 1, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {3};
        check_reachable(result, expected, 1, V, "TC5");
    }

    GrB_free(&result);
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}
//------------------------------------------------------------------------------
// TC6: S -> a S b | c
// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->2
// Source: 1   Expected: {2}
//------------------------------------------------------------------------------
void test_TC6_aSb_c_cycle2(void)
{
    setup();

    GrB_Matrix term[3], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = term[2] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_c(term, nt, call, start, final);

    GrB_Index V = 3;
    GrB_Matrix gterm[3];
    gterm[0] = gterm[1] = gterm[2] = GrB_NULL;
    for (int i = 0; i < 3; ++i)
        OK(GrB_Matrix_new(&gterm[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 2, 2));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 3, term, gterm, 1, nt, call, start, final, 0, 1, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {2};
        check_reachable(result, expected, 1, V, "TC6");
    }

    GrB_free(&result);
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->4; 4-b->2
static void graph_tc7_8(GrB_Matrix gterm[3])
{
    GrB_Index V = 5;
    for (int i = 0; i < 3; ++i)
        OK(GrB_Matrix_new(&gterm[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 3, 4));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 4, 2));
}

//------------------------------------------------------------------------------
// TC7: S -> a S b | c
// same graph as TC7/8 above
// Source: 0   Expected: {2, 3, 4}
//------------------------------------------------------------------------------
void test_TC7_aSb_c_cycle3(void)
{
    setup();

    GrB_Matrix term[3], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = term[2] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_c(term, nt, call, start, final);

    GrB_Index V = 5;
    GrB_Matrix gterm[3];
    gterm[0] = gterm[1] = gterm[2] = GrB_NULL;
    graph_tc7_8(gterm);

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 3, term, gterm, 1, nt, call, start, final, 0, 0, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {2, 3, 4};
        check_reachable(result, expected, 3, V, "TC7");
    }

    GrB_free(&result);
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

//------------------------------------------------------------------------------
// TC8: S -> a S b | c
// same graph, different source vertex
// Source: 1   Expected: {2, 3, 4}
//------------------------------------------------------------------------------
void test_TC8_aSb_c_cycle4(void)
{
    setup();

    GrB_Matrix term[3], nt[1], call[1];
    GrB_Vector start[1], final[1];

    term[0] = term[1] = term[2] = GrB_NULL;
    nt[0] = GrB_NULL;
    call[0] = GrB_NULL;
    start[0] = final[0] = GrB_NULL;

    rsm_aSb_c(term, nt, call, start, final);

    GrB_Index V = 5;
    GrB_Matrix gterm[3];
    gterm[0] = gterm[1] = gterm[2] = GrB_NULL;
    graph_tc7_8(gterm);

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 3, term, gterm, 1, nt, call, start, final, 0, 1, 4, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {2, 3, 4};
        check_reachable(result, expected, 3, V, "TC8");
    }

    GrB_free(&result);
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    GrB_free(&nt[0]);
    GrB_free(&call[0]);
    GrB_free(&start[0]);
    GrB_free(&final[0]);

    OK(info);

    teardown();
}

void test_TC9_mult_recursive_cycles(void)
{
    setup();

    GrB_Matrix term[4], nt[3], call[3];
    GrB_Vector start[3], final[3];

    term[0] = term[1] = term[2] = term[3] = GrB_NULL;
    nt[0] = nt[1] = nt[2] = GrB_NULL;
    call[0] = call[1] = call[2] = GrB_NULL;
    start[0] = start[1] = start[2] = GrB_NULL;
    final[0] = final[1] = final[2] = GrB_NULL;

    rsm_complex(term, nt, call, start, final);

    GrB_Index V = 6;
    GrB_Matrix gterm[4];
    gterm[0] = gterm[1] = gterm[2] = gterm[3] = GrB_NULL;
    for (int i = 0; i < 4; ++i)
        OK(GrB_Matrix_new(&gterm[i], GrB_BOOL, V, V));

    // a-cycle (0,1,2)
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 1, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[0], true, 2, 0));
    // b-cycle (2,3)
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(gterm[1], true, 3, 2));
    // c-cycle (3,4,5)
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 3, 4));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 4, 5));
    OK(GrB_Matrix_setElement_BOOL(gterm[2], true, 5, 3));
    // d-cycle
    OK(GrB_Matrix_setElement_BOOL(gterm[3], true, 5, 2));
    OK(GrB_Matrix_setElement_BOOL(gterm[3], true, 2, 5));

    GrB_Vector result = GrB_NULL;

    GrB_Info info =
        LAGraph_CFPQ_RSM(&result, 4, term, gterm, 3, nt, call, start, final, 0, 0, 11, V, msg);
    if (info == GrB_SUCCESS)
    {
        GrB_Index expected[] = {2, 3, 5};
        check_reachable(result, expected, 3, V, "TC9");
    }

    GrB_free(&result);
    for (int i = 0; i < 4; ++i)
    {
        GrB_free(&term[i]);
        GrB_free(&gterm[i]);
    }
    for (int i = 0; i < 3; ++i)
    {
        GrB_free(&nt[i]);
        GrB_free(&call[i]);
        GrB_free(&start[i]);
        GrB_free(&final[i]);
    }

    OK(info);

    teardown();
}

TEST_LIST = {{"TC1_cyclic_ab", test_TC1_cyclic_ab},
             {"TC2_self_loops", test_TC2_self_loops},
             {"TC3_single_node", test_TC3_single_node},
             {"TC4_aSb_c_loops", test_TC4_aSb_c_loops},
             {"TC5_aSb_c_cycle1", test_TC5_aSb_c_cycle1},
             {"TC4_aSb_c_loops", test_TC4_aSb_c_loops},
             {"TC7_aSb_c_cycle3", test_TC7_aSb_c_cycle3},
             {"TC8_aSb_c_cycle4", test_TC8_aSb_c_cycle4},
             {"TC9_mult_recursive_cycles", test_TC9_mult_recursive_cycles},
             {NULL, NULL}};
