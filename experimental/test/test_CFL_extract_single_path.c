#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>

#define run_aux_algorithm()                                                             \
    LAGraph_CFL_single_path(outputs, adj_matrices, grammar.terms_count,                 \
                            grammar.nonterms_count, grammar.rules, grammar.rules_count, \
                            msg)

#define run_algorithm()                                                                               \
    LAGraph_CFL_extract_single_path(&path, start, end, 0, adj_matrices, outputs, grammar.terms_count, \
                                    grammar.nonterms_count, grammar.rules, grammar.rules_count,       \
                                    msg)

// Check the path through the string
#define check_result1(expected_ret, result)             \
    {                                                   \
        retval = run_algorithm();                       \
        TEST_CHECK(retval == expected_ret);             \
        TEST_MSG("retval = %d (%s)", retval, msg);      \
        char *expected = path_to_str();                 \
        TEST_CHECK(strcmp(result, expected) == 0);      \
        TEST_MSG("Wrong result. Actual: %s", expected); \
        LAGraph_Free((void **)&expected, msg);          \
    }

// Check the path according to its type
#define check_result2(expected_path)                   \
    {                                                  \
        retval = run_algorithm();                      \
        if (expected_path == non_exist)                \
        {                                              \
            TEST_CHECK(retval == GrB_NO_VALUE);        \
            TEST_MSG("retval = %d (%s)", retval, msg); \
            TEST_CHECK(check_empty_path());            \
            TEST_MSG("Wrong result");                  \
        }                                              \
        else if (expected_path == empty)               \
        {                                              \
            TEST_CHECK(retval == GrB_SUCCESS);         \
            TEST_MSG("retval = %d (%s)", retval, msg); \
            TEST_CHECK(check_empty_path());            \
            TEST_MSG("Wrong result");                  \
        }                                              \
        else                                           \
        {                                              \
            TEST_CHECK(retval == GrB_SUCCESS);         \
            TEST_MSG("retval = %d (%s)", retval, msg); \
            TEST_CHECK(check_non_empty_path());        \
            TEST_MSG("Wrong result");                  \
        }                                              \
    }

typedef struct
{
    size_t nonterms_count;
    size_t terms_count;
    size_t rules_count;
    LAGraph_rule_WCNF *rules;
} grammar_t;

GrB_Matrix *adj_matrices = NULL;
int n_adj_matrices = 0;
GrB_Matrix *outputs = NULL;
grammar_t grammar = {0, 0, 0, NULL};
Path path;
char msg[LAGRAPH_MSG_LEN];

typedef enum // Path type
{
    empty = 0,     // Empty path - through a rule with eps
    non_empty = 1, // Path with non-zero length
    non_exist = 2, // Path does not exist
} Type_of_path;

void setup() { LAGraph_Init(msg); }

void teardown(void) { LAGraph_Finalize(msg); }

void init_outputs()
{
    LAGraph_Calloc((void **)&outputs,
                   grammar.nonterms_count, sizeof(GrB_Matrix), msg);
}

bool check_empty_path()
{
    return path.len == 0 && path.path == NULL;
}

bool check_non_empty_path()
{
    if (path.len == 0)
    {
        return false;
    }
    for (size_t i = 0; i < path.len; i++)
    {
        Edge cur_edge = path.path[i];
        bool edge_exist;
        if (GrB_Matrix_extractElement_BOOL(&edge_exist, adj_matrices[cur_edge.label], cur_edge.start, cur_edge.end) != GrB_SUCCESS)
        {
            return false;
        }
    }
    return true;
}

char *path_to_str()
{
    char *result_str = NULL;
    // 15 - size of "%ld->(%d)->%ld "
    // 15 - size of "len: %zu path: "
    // we need 16 + len * 15
    LAGraph_Malloc((void **)&result_str, 15 + path.len * 15, sizeof(char), msg);
    result_str[0] = '\0';

    sprintf(result_str + strlen(result_str), "len: %zu path: ", path.len);

    if (path.len > 0)
    {
        for (size_t i = 0; i < path.len; i++)
        {
            sprintf(result_str + strlen(result_str), "%" PRIu64 "->(%" PRId32 ")->%" PRIu64 " ", path.path[i].start, path.path[i].label, path.path[i].end);
        }
    }

    return result_str;
}

void free_workspace()
{

    if (adj_matrices != NULL)
    {
        for (size_t i = 0; i < n_adj_matrices; i++)
        {
            GrB_free(&adj_matrices[i]);
        }
    }
    LAGraph_Free((void **)&adj_matrices, msg);

    if (outputs != NULL)
    {
        for (size_t i = 0; i < grammar.nonterms_count; i++)
        {
            GrB_free(&outputs[i]);
        }
    }
    LAGraph_Free((void **)&outputs, msg);

    LAGraph_Free((void **)&grammar.rules, msg);
    grammar = (grammar_t){0, 0, 0, NULL};
}

//====================
// Grammars
//====================

// S -> aSb | ab in WCNF
//
// Terms: [0 a] [1 b]
// Nonterms: [0 S] [1 A] [2 B] [3 C]
// S -> AB [0 1 2 0]
// S -> AC [0 1 3 0]
// C -> SB [3 0 2 0]
// A -> a  [1 0 -1 0]
// B -> b  [2 1 -1 0]
void init_grammar_aSb()
{
    LAGraph_rule_WCNF *rules = NULL;
    LAGraph_Calloc((void **)&rules, 5, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 1, 2, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 1, 3, 0};
    rules[2] = (LAGraph_rule_WCNF){3, 0, 2, 0};
    rules[3] = (LAGraph_rule_WCNF){1, 0, -1, 0};
    rules[4] = (LAGraph_rule_WCNF){2, 1, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 4, .terms_count = 2, .rules_count = 5, .rules = rules};
}

// S -> aS | a | eps in WCNF
//
// Terms: [0 a]
// Nonterms: [0 S]
// S -> SS [0 0 0 0]
// S -> a  [0 0 -1 0]
// S -> eps [0 -1 -1 0]
void init_grammar_aS()
{
    LAGraph_rule_WCNF *rules = NULL;
    LAGraph_Calloc((void **)&rules, 3, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 0, 0, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 0, -1, 0};
    rules[2] = (LAGraph_rule_WCNF){0, -1, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 1, .terms_count = 1, .rules_count = 3, .rules = rules};
}

// Complex grammar
// aaaabbbb or aaabbb
//
// Terms: [0 a] [1 b]
// Nonterms: [0 S] [n Sn]
// S -> S1 S2       [0 1 2 0]
// S -> S15 S16     [0 15 16 0]
// S1 -> S3 S4      [1 3 4 0]
// S2 -> S5 S6      [2 5 6 0]
// S3 -> S7 S8      [3 7 8 0]
// S4 -> S9 S10     [4 9 10 0]
// S5 -> S11 S12    [5 11 12 0]
// S6 -> S13 S14    [6 13 14 0]
// S16 -> S17 S18   [16 17 18 0]
// S17 -> S19 S20   [17 19 20 0]
// S18 -> S21 S22   [18 21 22 0]
// S22 -> S23 S24   [22 23 24 0]
// S7 -> a          [7 0 -1 0]
// S8 -> a          [8 0 -1 0]
// S9 -> a          [9 0 -1 0]
// S10 -> a         [10 0 -1 0]
// S11 -> b         [11 1 -1 0]
// S12 -> b         [12 1 -1 0]
// S13 -> b         [13 1 -1 0]
// S14 -> b         [14 1 -1 0]
// S15 -> a         [15 0 -1 0]
// S19 -> a         [19 0 -1 0]
// S20 -> a         [20 0 -1 0]
// S21 -> b         [21 1 -1 0]
// S23 -> b         [23 1 -1 0]
// S24 -> b         [24 1 -1 0]
void init_grammar_complex()
{
    LAGraph_rule_WCNF *rules = NULL;
    LAGraph_Calloc((void **)&rules, 26, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 1, 2, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 15, 16, 0};
    rules[2] = (LAGraph_rule_WCNF){1, 3, 4, 0};
    rules[3] = (LAGraph_rule_WCNF){2, 5, 6, 0};
    rules[4] = (LAGraph_rule_WCNF){3, 7, 8, 0};
    rules[5] = (LAGraph_rule_WCNF){4, 9, 10, 0};
    rules[6] = (LAGraph_rule_WCNF){5, 11, 12, 0};
    rules[7] = (LAGraph_rule_WCNF){6, 13, 14, 0};
    rules[8] = (LAGraph_rule_WCNF){16, 17, 18, 0};
    rules[9] = (LAGraph_rule_WCNF){17, 19, 20, 0};
    rules[10] = (LAGraph_rule_WCNF){18, 21, 22, 0};
    rules[11] = (LAGraph_rule_WCNF){22, 23, 24, 0};
    rules[12] = (LAGraph_rule_WCNF){7, 0, -1, 0};
    rules[13] = (LAGraph_rule_WCNF){8, 0, -1, 0};
    rules[14] = (LAGraph_rule_WCNF){9, 0, -1, 0};
    rules[15] = (LAGraph_rule_WCNF){10, 0, -1, 0};
    rules[16] = (LAGraph_rule_WCNF){11, 1, -1, 0};
    rules[17] = (LAGraph_rule_WCNF){12, 1, -1, 0};
    rules[18] = (LAGraph_rule_WCNF){13, 1, -1, 0};
    rules[19] = (LAGraph_rule_WCNF){14, 1, -1, 0};
    rules[20] = (LAGraph_rule_WCNF){15, 0, -1, 0};
    rules[21] = (LAGraph_rule_WCNF){19, 0, -1, 0};
    rules[22] = (LAGraph_rule_WCNF){20, 0, -1, 0};
    rules[23] = (LAGraph_rule_WCNF){21, 1, -1, 0};
    rules[24] = (LAGraph_rule_WCNF){23, 1, -1, 0};
    rules[25] = (LAGraph_rule_WCNF){24, 1, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 25, .terms_count = 2, .rules_count = 26, .rules = rules};
}

//====================
// Graphs
//====================

// Graph:
//
// 0 -a-> 1
// 1 -a-> 2
// 2 -a-> 0
// 0 -b-> 3
// 3 -b-> 0
void init_graph_double_cycle()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    OK(GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 4, 4));
    OK(GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 4, 4));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 3));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 3, 0));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

// Graph:
//
// 0 -a-> 1
// 1 -a-> 2
// 2 -a-> 0
void init_graph_one_cycle()
{
    LAGraph_Calloc((void **)&adj_matrices, 1, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 1;

    GrB_Matrix adj_matrix_a;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 3, 3);

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));

    adj_matrices[0] = adj_matrix_a;
}

// Graph:
//
// 0 -a-> 1
// 1 -a-> 2
// 2 -a-> 3
// 3 -a-> 4
// 3 -b-> 5
// 4 -b-> 3
// 5 -b-> 6
// 6 -b-> 7
void init_graph_1()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    OK(GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 8, 8));
    OK(GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 8, 8));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 3));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 3, 4));

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 3, 5));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 4, 3));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 5, 6));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 6, 7));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

// Graph:
//
// 0 -a-> 2
// 1 -a-> 2
// 3 -a-> 5
// 4 -a-> 5
// 2 -a-> 6
// 5 -a-> 6
// 2 -b-> 0
// 2 -b-> 1
// 5 -b-> 3
// 5 -b-> 4
// 6 -b-> 2
// 6 -b-> 5
void init_graph_tree()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    OK(GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 7, 7));
    OK(GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 7, 7));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 3, 5));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 4, 5));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 6));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 5, 6));

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 2, 0));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 2, 1));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 5, 3));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 5, 4));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 6, 2));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 6, 5));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

// Graph:
//
// 0 -a-> 1
// 1 -a-> 2
// 2 -b-> 3
// 3 -b-> 4
void init_graph_line()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 5, 5);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 5, 5);

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 2, 3));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 3, 4));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

// Graph:
//
// 0 -a-> 0
// 0 -b-> 1
// 1 -c-> 2
void init_graph_2()
{
    LAGraph_Calloc((void **)&adj_matrices, 3, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 3;

    GrB_Matrix adj_matrix_a, adj_matrix_b, adj_matrix_c;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 3, 3);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 3, 3);
    GrB_Matrix_new(&adj_matrix_c, GrB_BOOL, 3, 3);

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 0));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_c, true, 1, 2));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
    adj_matrices[2] = adj_matrix_c;
}

// Graph:
//
// 0 -a-> 1
// 1 -a-> 0
// 0 -b-> 0
void init_graph_3()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 2, 2);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 2, 2);

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 0));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 0));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

// Graph:
//
// 0 -b-> 1
// 1 -b-> 0
void init_graph_4()
{
    LAGraph_Calloc((void **)&adj_matrices, 2, sizeof(GrB_Matrix), msg);
    n_adj_matrices = 2;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 2, 2);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 2, 2);

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 1, 0));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

//=====================
// Tests full result
//=====================

void test_CFL_extract_single_path_two_cycle(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs();
    int expected_ret[16] = {GrB_SUCCESS, GrB_NO_VALUE, GrB_NO_VALUE, GrB_SUCCESS, GrB_SUCCESS, GrB_NO_VALUE, GrB_NO_VALUE, GrB_SUCCESS, GrB_SUCCESS, GrB_NO_VALUE, GrB_NO_VALUE, GrB_SUCCESS, GrB_NO_VALUE, GrB_NO_VALUE, GrB_NO_VALUE, GrB_NO_VALUE};
    char *expected_path[16] = {"len: 12 path: 0->(0)->1 1->(0)->2 2->(0)->0 0->(0)->1 1->(0)->2 2->(0)->0 0->(1)->3 3->(1)->0 0->(1)->3 3->(1)->0 0->(1)->3 3->(1)->0 ",
                               "len: 0 path: ",
                               "len: 0 path: ",
                               "len: 6 path: 0->(0)->1 1->(0)->2 2->(0)->0 0->(1)->3 3->(1)->0 0->(1)->3 ",
                               "len: 4 path: 1->(0)->2 2->(0)->0 0->(1)->3 3->(1)->0 ",
                               "len: 0 path: ",
                               "len: 0 path: ",
                               "len: 10 path: 1->(0)->2 2->(0)->0 0->(0)->1 1->(0)->2 2->(0)->0 0->(1)->3 3->(1)->0 0->(1)->3 3->(1)->0 0->(1)->3 ",
                               "len: 8 path: 2->(0)->0 0->(0)->1 1->(0)->2 2->(0)->0 0->(1)->3 3->(1)->0 0->(1)->3 3->(1)->0 ",
                               "len: 0 path: ",
                               "len: 0 path: ",
                               "len: 2 path: 2->(0)->0 0->(1)->3 ",
                               "len: 0 path: ",
                               "len: 0 path: ",
                               "len: 0 path: ",
                               "len: 0 path: "};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 4; start++)
    {
        for (GrB_Index end = 0; end < 4; end++)
        {
            check_result1(expected_ret[start * 4 + end], expected_path[start * 4 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

//==========================================
// Tests that the path exists in the graph
//==========================================

void test_CFL_extract_single_path_cycle(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aS();
    init_graph_one_cycle();
    init_outputs();
    int expected[9] = {empty, non_empty, non_empty,
                       non_empty, empty, non_empty,
                       non_empty, non_empty, empty};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 3; start++)
    {
        for (GrB_Index end = 0; end < 3; end++)
        {
            check_result2(expected[start * 3 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_labels_more_than_nonterms(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_2();
    init_outputs();

    int expected[9] = {non_exist, non_empty, non_exist,
                       non_exist, non_exist, non_exist,
                       non_exist, non_exist, non_exist};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 3; start++)
    {
        for (GrB_Index end = 0; end < 3; end++)
        {
            check_result2(expected[start * 3 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_complex_grammar(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_complex();
    init_graph_1();
    init_outputs();
    int expected[64] = {non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_empty,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_empty, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 8; start++)
    {
        for (GrB_Index end = 0; end < 8; end++)
        {
            check_result2(expected[start * 8 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_tree(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_tree();
    init_outputs();

    int expected[49] = {non_empty, non_empty, non_exist, non_empty, non_empty, non_exist, non_exist,
                        non_empty, non_empty, non_exist, non_empty, non_empty, non_exist, non_exist,
                        non_exist, non_exist, non_empty, non_exist, non_exist, non_empty, non_exist,
                        non_empty, non_empty, non_exist, non_empty, non_empty, non_exist, non_exist,
                        non_empty, non_empty, non_exist, non_empty, non_empty, non_exist, non_exist,
                        non_exist, non_exist, non_empty, non_exist, non_exist, non_empty, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist, non_exist, non_exist};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 7; start++)
    {
        for (GrB_Index end = 0; end < 7; end++)
        {
            check_result2(expected[start * 7 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_line(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_line();
    init_outputs();
    int expected[25] = {non_exist, non_exist, non_exist, non_exist, non_empty,
                        non_exist, non_exist, non_exist, non_empty, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist,
                        non_exist, non_exist, non_exist, non_exist, non_exist};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 5; start++)
    {
        for (GrB_Index end = 0; end < 5; end++)
        {
            check_result2(expected[start * 5 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_two_nodes_cycle(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_3();
    init_outputs();
    int expected[4] = {non_empty, non_exist,
                       non_empty, non_exist};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 2; start++)
    {
        for (GrB_Index end = 0; end < 2; end++)
        {
            check_result2(expected[start * 2 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

void test_CFL_extract_single_path_with_empty_adj_matrix(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aS();
    init_graph_4();
    init_outputs();

    int expected[4] = {empty, non_exist,
                       non_exist, empty};
    OK(run_aux_algorithm());
    for (GrB_Index start = 0; start < 2; start++)
    {
        for (GrB_Index end = 0; end < 2; end++)
        {
            check_result2(expected[start * 2 + end]);
            if (path.len > 0)
            {
                LAGraph_Free((void **)&path.path, msg);
            }
        }
    }
    free_workspace();
    teardown();
#endif
}

TEST_LIST = {
    {"CFL_extract_single_path_two_cycle", test_CFL_extract_single_path_two_cycle},
    {"CFL_extract_single_path_cycle", test_CFL_extract_single_path_cycle},
    {"CFL_extract_single_path_labels_more_than_nonterms", test_CFL_extract_single_path_labels_more_than_nonterms},
    {"CFL_extract_single_path_complex_grammar", test_CFL_extract_single_path_complex_grammar},
    {"CFL_extract_single_path_tree", test_CFL_extract_single_path_tree},
    {"CFL_extract_single_path_line", test_CFL_extract_single_path_line},
    {"CFL_extract_single_path_two_nodes_cycle", test_CFL_extract_single_path_two_nodes_cycle},
    {"CFL_extract_single_path_with_empty_adj_matrix", test_CFL_extract_single_path_with_empty_adj_matrix},
    {NULL, NULL}};
