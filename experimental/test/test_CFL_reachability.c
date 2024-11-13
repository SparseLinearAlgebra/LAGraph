#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>

#define check_invalid_value()                                                            \
    {                                                                                    \
        retval = LAGraph_CFL_reachability(outputs, adj_matrices, grammar.terms_count,    \
                                          grammar.nonterms_count, grammar.rules,         \
                                          grammar.rules_count, msg);                     \
        TEST_CHECK(retval == GrB_INVALID_VALUE);                                         \
        TEST_MSG("retval = %d (%s)", retval, msg);                                       \
    }

typedef struct {
    size_t nonterms_count;
    size_t terms_count;
    size_t rules_count;
    LAGraph_rule_WCNF *rules;
} grammar_t;

GrB_Matrix *adj_matrices = NULL;
GrB_Matrix *outputs = NULL;
grammar_t grammar = {0, 0, 0, NULL};
char msg[LAGRAPH_MSG_LEN];

void setup() { LAGraph_Init(msg); }

void teardown(void) { LAGraph_Finalize(msg); }

void init_outputs() { outputs = calloc(grammar.nonterms_count, sizeof(GrB_Matrix)); }

void free_workspace() {
    for (size_t i = 0; i < grammar.terms_count; i++) {
        GrB_free(&adj_matrices[i]);
    }
    free(adj_matrices);
    adj_matrices = NULL;

    for (size_t i = 0; i < grammar.nonterms_count; i++) {
        GrB_free(&outputs[i]);
    }
    free(outputs);
    outputs = NULL;

    free(grammar.rules);
    grammar = (grammar_t){0, 0, 0, NULL};
}

//====================
// Grammars
//====================

// Terms: [0 a] [1 b]
// Nonterms: [0 S] [1 A] [2 B] [3 C]
// S -> AB [0 1 2 0]
// S -> AC [0 1 3 0]
// C -> SB [3 0 2 0]
// A -> a  [1 0 -1 0]
// B -> b  [2 1 -1 0]
grammar_t init_grammar_1() {
    LAGraph_rule_WCNF *rules = calloc(5, sizeof(LAGraph_rule_WCNF));
    rules[0] = (LAGraph_rule_WCNF){0, 1, 2, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 1, 3, 0};
    rules[2] = (LAGraph_rule_WCNF){3, 0, 2, 0};
    rules[3] = (LAGraph_rule_WCNF){1, 0, -1, 0};
    rules[4] = (LAGraph_rule_WCNF){2, 1, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 4, .terms_count = 2, .rules_count = 5, .rules = rules};
}

//====================
// Graphs
//====================

// Graph:
// 0 -a-> 1
// 1 -a-> 2
// 2 -a-> 0
// 0 -b-> 3
// 3 -b-> 0
void init_graph_1() {
    adj_matrices = calloc(2, sizeof(GrB_Matrix));
    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 4, 4);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 4, 4);
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 3));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 3, 0));
    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

void test_CFL_reachability_invalid_rules(void) {
    setup();
    GrB_Info retval;

    init_grammar_1();
    init_graph_1();
    init_outputs();

    // Rule [Variable -> _ B]
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = -1, .prod_B = 1, .index = 0};
    check_invalid_value();

    // Rule [_ -> A B]
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = -1, .prod_A = 1, .prod_B = 2, .index = 0};
    check_invalid_value();

    // Rule [C -> A B], where C >= nonterms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 10, .prod_A = 1, .prod_B = 2, .index = 0};
    check_invalid_value();

    // Rule [C -> t], where t >= terms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = 10, .prod_B = -1, .index = 0};
    check_invalid_value();

    free_workspace();
    teardown();

    return;
}

TEST_LIST = {{"CFG_reach_basic_invalid_rules", test_CFL_reachability_invalid_rules},
             {NULL, NULL}};