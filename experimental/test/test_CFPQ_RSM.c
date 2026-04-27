#include "LAGraph.h"
#include "LAGraphX.h"
#include "LAGraph_test.h"
#include <GraphBLAS.h>
#include <stddef.h>

char msg[LAGRAPH_MSG_LEN];

static RSM *rsm = NULL;
static size_t V = 0;
static size_t n_graph_term = 0;
static GrB_Matrix *graph_term = NULL;
static GrB_Vector result = GrB_NULL;

//==============================================================================
// test setup / teardown
//==============================================================================

static void setup(void) { LAGraph_Init(msg); }

static void teardown(void) { LAGraph_Finalize(msg); }

static RSM *alloc_rsm(GrB_Index state_count, GrB_Index term_count,
                      GrB_Index nonterm_count, GrB_Index start_nonterm) {
    RSM *r = NULL;
    LAGraph_Malloc((void **)&r, 1, sizeof(RSM), msg);
    if (r == NULL)
        return NULL;

    r->state_count = state_count;
    r->terminal_count = term_count;
    r->nonterminal_count = nonterm_count;
    r->start_nonterminal = start_nonterm;

    LAGraph_Calloc((void **)&r->terminal_matrices, term_count, sizeof(GrB_Matrix), msg);
    LAGraph_Calloc((void **)&r->nonterminal_matrices, nonterm_count, sizeof(GrB_Matrix),
                   msg);
    LAGraph_Calloc((void **)&r->start_states, nonterm_count, sizeof(GrB_Index), msg);
    LAGraph_Calloc((void **)&r->final_states, nonterm_count, sizeof(GrB_Vector), msg);

    for (size_t i = 0; i < term_count; ++i)
        GrB_Matrix_new(&r->terminal_matrices[i], GrB_BOOL, state_count, state_count);

    for (size_t i = 0; i < nonterm_count; ++i) {
        GrB_Matrix_new(&r->nonterminal_matrices[i], GrB_BOOL, state_count, state_count);
        GrB_Vector_new(&r->final_states[i], GrB_BOOL, state_count);
    }

    return r;
}

static void free_rsm(void) {
    if (rsm == NULL)
        return;

    if (rsm->terminal_matrices != NULL) {
        for (size_t i = 0; i < rsm->terminal_count; ++i)
            GrB_free(&rsm->terminal_matrices[i]);
        LAGraph_Free((void **)&rsm->terminal_matrices, msg);
    }
    if (rsm->nonterminal_matrices != NULL) {
        for (size_t i = 0; i < rsm->nonterminal_count; ++i)
            GrB_free(&rsm->nonterminal_matrices[i]);
        LAGraph_Free((void **)&rsm->nonterminal_matrices, msg);
    }

    if (rsm->final_states != NULL) {
        for (size_t i = 0; i < rsm->nonterminal_count; ++i) {
            GrB_free(&rsm->final_states[i]);
        }
        LAGraph_Free((void **)&rsm->final_states, msg);
    }

    LAGraph_Free((void **)&rsm->start_states, msg);
    LAGraph_Free((void **)&rsm, msg);
}

static void free_graph(void) {
    if (graph_term == NULL)
        return;
    for (size_t i = 0; i < n_graph_term; ++i)
        GrB_free(&graph_term[i]);
    LAGraph_Free((void **)&graph_term, msg);
    n_graph_term = 0;
    V = 0;
}

static void free_workspace(void) {
    GrB_free(&result);
    free_rsm();
    free_graph();
}

static void check_reachable(GrB_Index *expected, GrB_Index n_expected,
                            const char *test_name) {
    GrB_Index nvals = 0;
    OK(GrB_Vector_nvals(&nvals, result));

    bool *got = (bool *)calloc(V, sizeof(bool));
    TEST_CHECK(got != NULL);

    GrB_Index *idx = GrB_NULL;
    bool *val = GrB_NULL;

    if (nvals > 0) {
        idx = (GrB_Index *)malloc(nvals * sizeof(GrB_Index));
        val = (bool *)malloc(nvals * sizeof(bool));
        if (!TEST_CHECK(idx != NULL && val != NULL)) {
            free(got);
            free(idx);
            free(val);
            return;
        }
        OK(GrB_Vector_extractTuples_BOOL(idx, val, &nvals, result));
        for (GrB_Index k = 0; k < nvals; ++k)
            if (val[k])
                got[idx[k]] = true;
        free(idx);
        free(val);
    }

    // every expected vertex must be reachable
    for (GrB_Index k = 0; k < n_expected; ++k) {
        TEST_CHECK(got[expected[k]]);
        TEST_MSG("%s: vertex %llu should be reachable", test_name,
                 (unsigned long long)expected[k]);
    }

    GrB_Index actual_count = 0;
    for (GrB_Index v = 0; v < V; ++v)
        if (got[v])
            actual_count++;

    TEST_CHECK(actual_count == n_expected);
    TEST_MSG("%s: expected %llu reachable vertices, got %llu", test_name,
             (unsigned long long)n_expected, (unsigned long long)actual_count);

    free(got);
}

//==============================================================================
// RSM builders
//==============================================================================
// Grammar 1:  S -> a S b | a b
static void init_rsm_aSb_ab(void) {
    GrB_Index Q = 4;
    rsm = alloc_rsm(Q, /*term_count=*/2, /*nonterm_count=*/1, /*start_nonterm=*/0);

    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[0], true, 0, 1));    // 0 -a-> 1
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[1], true, 1, 3));    // 1 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[1], true, 2, 3));    // 2 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[0], true, 1, 2)); // 1 -S-> 2
    rsm->start_states[0] = 0;
    OK(GrB_Vector_setElement_BOOL(rsm->final_states[0], true, 3));
}

// Grammar 2:  S -> a S b | c
static void init_rsm_aSb_c(void) {
    GrB_Index Q = 4;
    rsm = alloc_rsm(Q, /*term_count=*/3, /*nonterm_count=*/1, /*start_nonterm=*/0);

    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[0], true, 0, 1));    // 0 -a-> 1
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[2], true, 0, 3));    // 0 -c-> 3
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[1], true, 2, 3));    // 2 -b-> 3
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[0], true, 1, 2)); // 1 -S-> 2
    rsm->start_states[0] = 0;
    OK(GrB_Vector_setElement_BOOL(rsm->final_states[0], true, 3));
}

// Grammar 3: S -> X | Y | S X | S Y
//            X -> a S b | a b
//            Y -> c S d | c d
static void init_rsm_complex(void) {
    GrB_Index Q = 11;
    rsm = alloc_rsm(Q, /*term_count=*/4, /*nonterm_count=*/3, /*start_nonterm=*/0);

    // --- S Component (States 0-2) ---
    // S -> X | Y | S X | S Y
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[1], true, 0, 2)); // 0 -X-> 2
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[2], true, 0, 2)); // 0 -Y-> 2
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[0], true, 0, 1)); // 0 -S-> 1
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[1], true, 1, 2)); // 1 -X-> 2
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[2], true, 1, 2)); // 1 -Y-> 2
    rsm->start_states[0] = 0;
    OK(GrB_Vector_setElement_BOOL(rsm->final_states[0], true, 2));

    // --- X Component (States 3-6) ---
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[0], true, 3, 4));    // 3 -a-> 4
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[1], true, 4, 6));    // 4 -b-> 6
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[1], true, 5, 6));    // 5 -b-> 6
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[0], true, 4, 5)); // 4 -S-> 5
    rsm->start_states[1] = 3;
    OK(GrB_Vector_setElement_BOOL(rsm->final_states[1], true, 6));

    // --- Y Component (States 7-10)
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[2], true, 7, 8));    // 7 -c-> 8
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[3], true, 8, 10));   // 8 -d-> 10
    OK(GrB_Matrix_setElement_BOOL(rsm->terminal_matrices[3], true, 9, 10));   // 9 -d-> 10
    OK(GrB_Matrix_setElement_BOOL(rsm->nonterminal_matrices[0], true, 8, 9)); // 8 -S-> 9
    rsm->start_states[2] = 7;
    OK(GrB_Vector_setElement_BOOL(rsm->final_states[2], true, 10));
}

//==============================================================================
// Graph builders
//==============================================================================
// 0-a->1; 1-a->2; 2-a->0; 0-b->3; 3-b->0 (V=4, 2 labels: a=0, b=1)
static void init_graph_double_cycle_ab(void) {
    V = 4;
    n_graph_term = 2;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 1, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 2, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 0, 3));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 3, 0));
}

// 0-a->0; 0-a->1; 1-b->1 (V=2, 2 labels: a=0, b=1)
static void init_graph_self_loops(void) {
    V = 2;
    n_graph_term = 2;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 1, 1));
}

// 0-a->0; 0-b->0 (V=2, 2 labels: a=0, b=1)
static void init_graph_single_node(void) {
    V = 2;
    n_graph_term = 2;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 0, 0));
}

// 0-a->0; 0-c->1; 1-b->1   (V=2, 3 labels: a=0, b=1, c=2)
static void init_graph_loops_abc(void) {
    V = 2;
    n_graph_term = 3;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 0, 0));
}

// 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->2   (V=4, 3 labels: a=0, b=1, c=2)
static void init_graph_loops_abc_1(void) {
    V = 4;
    n_graph_term = 3;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 3, 2));
}

// 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->4; 4-b->2   (V=5, 3 labels: a=0, b=1, c=2)
static void init_graph_loops_abc_2(void) {
    V = 5;
    n_graph_term = 3;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 3, 4));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 4, 2));
}

// a-cycle(0,1,2), b-cycle(2,3), c-cycle(3,4,5), d-edges(2<->5)   (V=6, 4 labels)
static void init_graph_multi_cycle(void) {
    V = 6;
    n_graph_term = 4;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 1, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 2, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 2, 3));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 3, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 3, 4));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 4, 5));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 5, 3));
    OK(GrB_Matrix_setElement_BOOL(graph_term[3], true, 5, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[3], true, 2, 5));
}

#define run_algorithm(sources, num_sources)                                              \
    LAGraph_CFPQ_RSM(&result, rsm, graph_term, (sources), (num_sources), V, msg);

//==============================================================================
// test cases
//==============================================================================

//------------------------------------------------------------------------------
// TC1: S -> a S b | a b
// Graph: 0-a->1; 1-a->2; 2-a->0; 0-b->3; 3-b->0
// Source: 1   Expected: {0, 3}
//------------------------------------------------------------------------------
void test_TC1_cyclic_ab(void) {
    setup();

    init_rsm_aSb_ab();
    init_graph_double_cycle_ab();

    GrB_Index sources = {1};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {0, 3};
        check_reachable(expected, /*n_expected=*/2, "TC1");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC2: S -> a S b | a b
// Graph: 0-a->0; 0-a->1; 1-b->1
// Source: 0   Expected: {1}
//------------------------------------------------------------------------------
void test_TC2_self_loops(void) {
    setup();
    init_rsm_aSb_ab();
    init_graph_self_loops();

    GrB_Index sources = {0};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {1};
        check_reachable(expected, /*n_expected=*/1, "TC2");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC3: S -> a S b | a b
// Graph: 0-a->0; 0-b->0
// Source: 0   Expected: {0}
//------------------------------------------------------------------------------
void test_TC3_single_node(void) {
    setup();
    init_rsm_aSb_ab();
    init_graph_single_node();

    GrB_Index sources = {0};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {0};
        check_reachable(expected, /*n_expected=*/1, "TC3");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC4: S -> a S b | c
// Graph: 0-a->0; 0-c->1; 1-b->1
// Source: 0   Expected: {1}
//------------------------------------------------------------------------------
void test_TC4_aSb_c_loops(void) {
    setup();
    init_rsm_aSb_c();
    init_graph_loops_abc();

    GrB_Index sources = {0};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {1};
        check_reachable(expected, 1, "TC4");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC5: S -> a S b | c
// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->2
// Source: 1   Expected: {3}
//------------------------------------------------------------------------------
void test_TC5_aSb_c_cycle1(void) {
    setup();
    init_rsm_aSb_c();
    init_graph_loops_abc_1();

    GrB_Index sources = {1};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {3};
        check_reachable(expected, 1, "TC5");
    }

    free_workspace();
    OK(info);
    teardown();
}
//------------------------------------------------------------------------------
// TC6: S -> a S b | c
// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->2
// Source: 1   Expected: {2}
//------------------------------------------------------------------------------
void test_TC6_aSb_c_cycle2(void) {
    setup();
    init_rsm_aSb_c();

    V = 3;
    n_graph_term = 3;
    LAGraph_Calloc((void **)&graph_term, n_graph_term, sizeof(GrB_Matrix), msg);
    for (size_t i = 0; i < n_graph_term; ++i)
        OK(GrB_Matrix_new(&graph_term[i], GrB_BOOL, V, V));

    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 0, 1));
    OK(GrB_Matrix_setElement_BOOL(graph_term[0], true, 1, 0));
    OK(GrB_Matrix_setElement_BOOL(graph_term[2], true, 0, 2));
    OK(GrB_Matrix_setElement_BOOL(graph_term[1], true, 2, 2));

    GrB_Index sources = {1};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {2};
        check_reachable(expected, 1, "TC6");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC7: S -> a S b | c
// Graph: 0-a->1; 1-a->0; 0-c->2; 2-b->3; 3-b->4; 4-b->2
// Source: 0   Expected: {2, 3, 4}
//------------------------------------------------------------------------------
void test_TC7_aSb_c_cycle3(void) {
    setup();
    init_rsm_aSb_c();
    init_graph_loops_abc_2();

    GrB_Index sources = {0};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {2, 3, 4};
        check_reachable(expected, 3, "TC7");
    }

    free_workspace();
    OK(info);
    teardown();
}

//------------------------------------------------------------------------------
// TC8: S -> a S b | c
// same graph, different source vertex
// Source: 1   Expected: {2, 3, 4}
//------------------------------------------------------------------------------
void test_TC8_aSb_c_cycle4(void) {
    setup();
    init_rsm_aSb_c();
    init_graph_loops_abc_2();

    GrB_Index sources = {1};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {2, 3, 4};
        check_reachable(expected, 3, "TC8");
    }

    free_workspace();
    OK(info);
    teardown();
}

void test_TC9_mult_recursive_cycles(void) {
    setup();

    init_rsm_complex();
    init_graph_multi_cycle();

    GrB_Index sources = {1};
    GrB_Info info = run_algorithm(&sources, 1);
    if (info == GrB_SUCCESS) {
        GrB_Index expected[] = {2, 3, 5};
        check_reachable(expected, 3, "TC9");
    }

    free_workspace();
    OK(info);
    teardown();
}

TEST_LIST = {{"TC1_cyclic_ab", test_TC1_cyclic_ab},
             {"TC2_self_loops", test_TC2_self_loops},
             {"TC3_single_node", test_TC3_single_node},
             {"TC4_aSb_c_loops", test_TC4_aSb_c_loops},
             {"TC5_aSb_c_cycle1", test_TC5_aSb_c_cycle1},
             {"TC6_aSb_c_cycle2", test_TC6_aSb_c_cycle2},
             {"TC7_aSb_c_cycle3", test_TC7_aSb_c_cycle3},
             {"TC8_aSb_c_cycle4", test_TC8_aSb_c_cycle4},
             {"TC9_mult_recursive_cycles", test_TC9_mult_recursive_cycles},
             {NULL, NULL}};