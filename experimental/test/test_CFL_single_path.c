#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>

#define run_algorithm()                                                          \
    LAGraph_CFL_single_path(outputs, adj_matrices, grammar.terms_count,          \
                     grammar.nonterms_count, grammar.rules, grammar.rules_count, \
                     msg)

#define check_error(error)                         \
    {                                              \
        retval = run_algorithm();                  \
        TEST_CHECK(retval == error);               \
        TEST_MSG("retval = %d (%s)", retval, msg); \
    }

#define check_result(result)                            \
    {                                                   \
        char *expected = output_to_str(0);              \
        TEST_CHECK(strcmp(result, expected) == 0);      \
        TEST_MSG("Wrong result. Actual: %s", expected); \
        LAGraph_Free((void **)&expected, msg);          \
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
char msg[LAGRAPH_MSG_LEN];

void setup() { LAGraph_Init(msg); }

void teardown(void) { LAGraph_Finalize(msg); }

void init_outputs()
{
    LAGraph_Calloc((void **)&outputs,
                   grammar.nonterms_count, sizeof(GrB_Matrix), msg);
}

char *output_to_str(size_t nonterm)
{
    GrB_Index nnz = 0;
    GrB_Matrix_nvals(&nnz, outputs[nonterm]);
    GrB_Index *row = NULL;
    GrB_Index *col = NULL;
    PathIndex *val = NULL;
    LAGraph_Malloc((void **)&row, nnz, sizeof(GrB_Index), msg);
    LAGraph_Malloc((void **)&col, nnz, sizeof(GrB_Index), msg);
    LAGraph_Malloc((void **)&val, nnz, sizeof(PathIndex), msg);

    GrB_Matrix_extractTuples_UDT(row, col, val, &nnz, outputs[nonterm]);

    char *result_str = NULL;
    LAGraph_Malloc((void **)&result_str, nnz * 40, sizeof(char), msg);
    result_str[0] = '\0';

    for (GrB_Index i = 0; i < nnz; i++)
    {
        char buf[64];
        sprintf(buf,
                i == 0
                    ? "(%" PRIu64 ", %" PRIu64 "): m=%ld h=%d"
                    : " (%" PRIu64 ", %" PRIu64 "): m=%ld h=%d",
                row[i], col[i],
                val[i].middle,
                val[i].height);

        strcat(result_str, buf);
    }
    LAGraph_Free((void **)&row, msg);
    LAGraph_Free((void **)&col, msg);
    LAGraph_Free((void **)&val, msg);

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

void test_CFL_single_path_cycle(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aS();
    init_graph_one_cycle();
    init_outputs();

    OK(run_algorithm());
    check_result("(0, 0): m=0 h=1 (0, 1): m=0 h=1 (0, 2): m=1 h=2 (1, 0): m=2 h=2 (1, 1): m=1 h=1 (1, 2): m=1 h=1 (2, 0): m=2 h=1 (2, 1): m=0 h=2 (2, 2): m=2 h=1");

    free_workspace();
    teardown();
#endif
}

void test_CFL_single_path_two_cycle(void)
{
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs();

    OK(run_algorithm());
    check_result("(0, 0): m=1 h=12 (0, 3): m=1 h=6 (1, 0): m=2 h=4 (1, 3): m=2 h=10 (2, 0): m=0 h=8 (2, 3): m=0 h=2");

    free_workspace();
    teardown();
#endif
}

TEST_LIST = {
    {"CFL_reachability_cycle", test_CFL_single_path_cycle},
    {"CFL_path_two_cycle", test_CFL_single_path_two_cycle},
    {NULL, NULL}};
