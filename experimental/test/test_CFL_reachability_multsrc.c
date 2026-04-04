//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_CFL_reachability_multsrc.c: test cases for
// Multiple-Source Context-Free Language Reachability Matrix-Based Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vlasenco Daniel, Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State University.

//------------------------------------------------------------------------------

/*
 * Note:
 * Tests than end with *_allsrc are the same as tests for CFL_reachability, with
 * all vertices selected as sources. Their behaviour should be identical to CFL_reachability.
*/

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>

#define run_algorithm()                                                                  \
    LAGraph_CFL_reachability_multsrc(&output, adj_matrices, src, src_count,              \
        grammar.terms_count, grammar.nonterms_count, grammar.rules, grammar.rules_count, \
        msg)

#define check_error(error)                                                               \
    {                                                                                    \
        retval = run_algorithm();                                                        \
        TEST_CHECK(retval == error);                                                     \
        TEST_MSG("%d (retval) != %d (error) (%s)", retval, error, msg);                                       \
    }

#define check_result(result)                                                             \
    {                                                                                    \
        char *expected = output_to_str(0);                                               \
        TEST_CHECK(strcmp(result, expected) == 0);                                       \
        TEST_MSG("Wrong result. Actual: %s", expected);                                  \
    }

typedef struct {
    size_t nonterms_count;
    size_t terms_count;
    size_t rules_count;
    LAGraph_rule_WCNF *rules;
} grammar_t;

GrB_Matrix *adj_matrices = NULL;
GrB_Matrix output = NULL;
grammar_t grammar = {0, 0, 0, NULL};
char msg[LAGRAPH_MSG_LEN];

void setup() { LAGraph_Init(msg); }

void teardown(void) { LAGraph_Finalize(msg); }

char *output_to_str(size_t nonterm) {
    GrB_Index nnz = 0;
    OK(GrB_Matrix_nvals(&nnz, output));
    GrB_Index *row = NULL ;
    GrB_Index *col = NULL ;
    bool *val = NULL ;
    LAGraph_Malloc ((void **) &row, nnz, sizeof (GrB_Index), msg) ;
    LAGraph_Malloc ((void **) &col, nnz, sizeof (GrB_Index), msg) ;
    LAGraph_Malloc ((void **) &val, nnz, sizeof (GrB_Index), msg) ;

    OK(GrB_Matrix_extractTuples(row, col, val, &nnz, output));

    // 11 - size of " (%ld, %ld)"
//  char *result_str = malloc(11 * nnz * sizeof(char));
    char *result_str = NULL ;
    LAGraph_Malloc ((void **) &result_str, 11*nnz, sizeof (char), msg) ;

    result_str[0] = '\0';
    for (size_t i = 0; i < nnz; i++) {
        sprintf(result_str + strlen(result_str), i == 0 ?
            "(%" PRIu64 ", %" PRIu64 ")" : " (%" PRIu64 ", %" PRIu64 ")",
            row[i], col[i]);
    }

    LAGraph_Free ((void **) &row, msg);
    LAGraph_Free ((void **) &col, msg);
    LAGraph_Free ((void **) &val, msg);

    return result_str;
}

void free_workspace() {
    for (size_t i = 0; i < grammar.terms_count; i++) {
        GrB_free(&adj_matrices[i]);
    }
    LAGraph_Free ((void **) &adj_matrices, msg);

    GrB_free(&output);

    LAGraph_Free ((void **) &grammar.rules, msg);
    grammar = (grammar_t){0, 0, 0, NULL};
}

//================================
// Grammars
//================================

// S -> aSb | ab in WCNF
//
// Terms: [0 a] [1 b]
// Nonterms: [0 S] [1 A] [2 B] [3 C]
// S -> AB [0 1 2 0]
// S -> AC [0 1 3 0]
// C -> SB [3 0 2 0]
// A -> a  [1 0 -1 0]
// B -> b  [2 1 -1 0]
void init_grammar_aSb() {
//  LAGraph_rule_WCNF *rules = calloc(5, sizeof(LAGraph_rule_WCNF));
    LAGraph_rule_WCNF *rules = NULL ;
    LAGraph_Calloc ((void **) &rules, 5, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 1, 2, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 1, 3, 0};
    rules[2] = (LAGraph_rule_WCNF){3, 0, 2, 0};
    rules[3] = (LAGraph_rule_WCNF){1, 0, -1, 0};
    rules[4] = (LAGraph_rule_WCNF){2, 1, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 4, .terms_count = 2, .rules_count = 5, .rules = rules};
}

// S -> aSb | ab | eps in WCNF
//
// Terms: [0 a] [1 b]
// Nonterms: [0 S] [1 A] [2 B] [3 C]
// S -> AB  [0 1 2 0]
// S -> AC  [0 1 3 0]
// S -> eps [0 -1 -1 0]
// C -> SB  [3 0 2 0]
// A -> a   [1 0 -1 0]
// B -> b   [2 1 -1 0]
void init_grammar_aSb_eps() {
//  LAGraph_rule_WCNF *rules = calloc(5, sizeof(LAGraph_rule_WCNF));
    LAGraph_rule_WCNF *rules = NULL ;
    LAGraph_Calloc ((void **) &rules, 6, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 1, 2, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 1, 3, 0};
    rules[2] = (LAGraph_rule_WCNF){3, 0, 2, 0};
    rules[3] = (LAGraph_rule_WCNF){1, 0, -1, 0};
    rules[4] = (LAGraph_rule_WCNF){2, 1, -1, 0};
    rules[5] = (LAGraph_rule_WCNF){0, -1, -1, 0};


    grammar = (grammar_t){
        .nonterms_count = 4, .terms_count = 2, .rules_count = 6, .rules = rules};
}

// S -> aS | a in WCNF
//
// Terms: [0 a]
// Nonterms: [0 S]
// S -> SS [0 0 0 0]
// S -> a  [0 0 -1 0]
void init_grammar_aS() {
//  LAGraph_rule_WCNF *rules = calloc(2, sizeof(LAGraph_rule_WCNF));
    LAGraph_rule_WCNF *rules = NULL ;
    LAGraph_Calloc ((void **) &rules, 2, sizeof(LAGraph_rule_WCNF), msg);

    rules[0] = (LAGraph_rule_WCNF){0, 0, 0, 0};
    rules[1] = (LAGraph_rule_WCNF){0, 0, -1, 0};

    grammar = (grammar_t){
        .nonterms_count = 1, .terms_count = 1, .rules_count = 2, .rules = rules};
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
void init_grammar_complex() {
//  LAGraph_rule_WCNF *rules = calloc(26, sizeof(LAGraph_rule_WCNF));
    LAGraph_rule_WCNF *rules = NULL ;
    LAGraph_Calloc ((void **) &rules, 26, sizeof(LAGraph_rule_WCNF), msg);

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

//================================
// Graphs without vertex labels
//================================

// Graph:
//
// 0 -a-> 1
// 1 -a-> 2
// 2 -a-> 0
// 0 -b-> 3
// 3 -b-> 0
void init_graph_double_cycle() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;

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
// 2 -a-> 3
// 3 -a-> 4
// 3 -b-> 5
// 4 -b-> 3
// 5 -b-> 6
// 6 -b-> 7
void init_graph_1() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;

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
void init_graph_tree() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;

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
// 2 -a-> 0
void init_graph_one_cycle() {
//  adj_matrices = calloc(1, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 1, sizeof (GrB_Matrix), msg) ;

    GrB_Matrix adj_matrix_a;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 3, 3);

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));

    adj_matrices[0] = adj_matrix_a;
}

// Graph:

// 0 -a-> 1
// 1 -a-> 2
// 2 -b-> 3
// 3 -b-> 4
void init_graph_line() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;

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

// 0 -a-> 0
// 0 -b-> 1
// 1 -c-> 2
void init_graph_2() {
//  adj_matrices = calloc(3, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 3, sizeof (GrB_Matrix), msg) ;

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

// 0 -a-> 1
// 1 -a-> 0
// 0 -b-> 0
void init_graph_3() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;

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
// 1 -a-> 0
// 2 -a-> 0
// 3 -a-> 0
// 4 -a-> 0
// 5 -a-> 0
// 1 -a-> 2
// 2 -a-> 3
// 3 -a-> 4
// 4 -a-> 5
// 5 -a-> 1
void init_graph_whirlpool() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 1, sizeof (GrB_Matrix), msg) ;

    GrB_Matrix adj_matrix_a;
    OK(GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 6, 6));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 0));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 3, 0));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 4, 0));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 5, 0));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 2, 3));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 3, 4));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 4, 5));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 5, 1));

    adj_matrices[0] = adj_matrix_a;
}

// Graph:
//
// 0 -a-> 1
// 0 -a-> 2
// 0 -a-> 3
// 0 -a-> 4
// 0 -b-> 0
void init_graph_allout() {
//  adj_matrices = calloc(2, sizeof(GrB_Matrix));
    LAGraph_Calloc ((void **) &adj_matrices, 1, sizeof (GrB_Matrix), msg) ;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    OK(GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 5, 5));

    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 3));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 2));
    OK(GrB_Matrix_setElement(adj_matrix_a, true, 0, 4));

    adj_matrices[0] = adj_matrix_a;
}

//================================
// Tests with all vertices as sources
//================================

void test_CFL_reachability_cycle_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_one_cycle();

    OK(run_algorithm());
    check_result("(0, 0) (0, 1) (0, 2) (1, 0) (1, 1) (1, 2) (2, 0) (2, 1) (2, 2)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_two_cycle_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2, 3 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_double_cycle();

    OK(run_algorithm());
    check_result("(0, 0) (0, 3) (1, 0) (1, 3) (2, 0) (2, 3)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_labels_more_than_nonterms_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_2();

    OK(run_algorithm());
    check_result("(0, 1)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_complex_grammar_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_complex();
    init_graph_1();

    OK(run_algorithm());
    check_result("(0, 7) (1, 6)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_tree_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2, 3, 4, 5, 6 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_tree();

    OK(run_algorithm());
    check_result("(0, 0) (0, 1) (0, 3) (0, 4) (1, 0) (1, 1) (1, 3) (1, 4) (2, 2) (2, 5) "
                 "(3, 0) (3, 1) (3, 3) (3, 4) (4, 0) (4, 1) (4, 3) (4, 4) (5, 2) (5, 5)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_line_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1, 2, 3, 4 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_line();

    OK(run_algorithm());
    check_result("(0, 4) (1, 3)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_two_nodes_cycle_allsrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0, 1 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_3();

    OK(run_algorithm());
    check_result("(0, 0) (1, 0)");

    free_workspace();
    teardown();
}

//================================
// Tests with some vertices as sources
//================================

void test_CFL_reachability_tree_msrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 2, 3 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_tree();

    OK(run_algorithm());
    check_result("(2, 2) (2, 5) (3, 0) (3, 1) (3, 3) (3, 4)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_allin_1_4(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 1, 4 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_whirlpool();

    OK(run_algorithm());
    check_result(
        "(1, 0) (1, 1) (1, 2) (1, 3) (1, 4) (1, 5) (4, 0) (4, 1) (4, 2) (4, 3) (4, 4) (4, 5)");

    free_workspace();
    teardown();
}

//================================
// Tests with one source vertex
//================================

void test_CFL_reachability_cycle_onesrc(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_one_cycle();

    OK(run_algorithm());
    check_result("(0, 0) (0, 1) (0, 2)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_allin_1(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_whirlpool();

    OK(run_algorithm());
    check_result("");

    free_workspace();
    teardown();
}

void test_CFL_reachability_allout_0(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 0 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_allout();

    OK(run_algorithm());
    check_result("(0, 1) (0, 2) (0, 3) (0, 4)");

    free_workspace();
    teardown();
}

void test_CFL_reachability_allout_1(void) {
    setup();
    GrB_Info retval;
    GrB_Index src[] = { 1 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aS();
    init_graph_allout();

    OK(run_algorithm());
    check_result("");

    free_workspace();
    teardown();
}

//====================
// Tests with invalid parameters
//====================

void test_CFL_reachability_invalid_rules(void) {
    setup();
    GrB_Info retval;

    GrB_Index src[] = { 0 };
    int32_t src_count = sizeof(src) / sizeof(GrB_Index);

    init_grammar_aSb_eps();
    init_graph_double_cycle();

    // Rule [Variable -> _ B]
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = -1, .prod_B = 1, .index = 0};
    check_error(GrB_INVALID_VALUE);

    // Rule [_ -> A B]
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = -1, .prod_A = 1, .prod_B = 2, .index = 0};
    check_error(GrB_INVALID_VALUE);

    // Rule [C -> A B], where C >= nonterms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 10, .prod_A = 1, .prod_B = 2, .index = 0};
    check_error(GrB_INVALID_VALUE);

    // Rule [C -> t], where t >= terms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = 10, .prod_B = -1, .index = 0};
    check_error(GrB_INVALID_VALUE);

    free_workspace();
    teardown();
}

// At least one rule of form S -> _ and
// at least one source node should be present.
// nonterms_count, rules_count, src_count > 0
void test_CFL_reachability_invalid_params(void) {
    setup();
    GrB_Info retval;

    GrB_Index src[] = { 0, 1, 2 };

    init_grammar_aSb();
    init_graph_double_cycle();

    int32_t src_count = 0;
    check_error(GrB_INVALID_VALUE);
    src_count = sizeof(src) / sizeof(GrB_Index);

    int32_t temp = grammar.nonterms_count;
    grammar.nonterms_count = 0;
    check_error(GrB_INVALID_VALUE);
    grammar.nonterms_count = temp;

    temp = grammar.rules_count;
    grammar.rules_count = 0;
    check_error(GrB_INVALID_VALUE);
    grammar.rules_count = temp;

    free_workspace();
    teardown();
}

// output, rules, adj_matrices, src should not be NULL
void test_CFL_reachability_null_params(void) {
    setup();
    GrB_Info retval;

    GrB_Index _src[] = { 0, 1, 2 };
    GrB_Index* src = _src;
    int32_t src_count = sizeof(_src) / sizeof(GrB_Index);

    init_grammar_aSb();
    init_graph_double_cycle();

    GrB_Matrix *temp_m = adj_matrices;
    adj_matrices = NULL;
    check_error(GrB_NULL_POINTER);
    adj_matrices = temp_m;

    GrB_Index *temp_src = src;
    src = NULL;
    check_error(GrB_NULL_POINTER);
    src = temp_src;

    LAGraph_rule_WCNF *temp_rules = grammar.rules;
    grammar.rules = NULL;
    check_error(GrB_NULL_POINTER);
    grammar.rules = temp_rules;

    #define run_algorithm()                                                                  \
    LAGraph_CFL_reachability_multsrc(NULL, adj_matrices, src, src_count,              \
        grammar.terms_count, grammar.nonterms_count, grammar.rules, grammar.rules_count, \
        msg)
    check_error(GrB_NULL_POINTER);

    free_workspace();
    teardown();
}


TEST_LIST = {
             // Tests from LAGraph_CFL_reachability as special cases
             {"test_CFL_reachability_cycle_allsrc", test_CFL_reachability_cycle_allsrc},
             {"CFL_reachability_two_cycle_allsrc", test_CFL_reachability_two_cycle_allsrc},
             {"CFL_reachability_labels_more_than_nonterms_allsrc", test_CFL_reachability_labels_more_than_nonterms_allsrc},
             {"CFL_reachability_complex_grammar_allsrc", test_CFL_reachability_complex_grammar_allsrc},
             {"CFL_reachability_tree_allsrc", test_CFL_reachability_tree_allsrc},
             {"CFL_reachability_line_allsrc", test_CFL_reachability_line_allsrc},
             {"CFL_reachability_two_nodes_cycle_allsrc", test_CFL_reachability_two_nodes_cycle_allsrc},
             //  Tests for some source vertices.
             {"CFL_reachability_cycle_onesrc", test_CFL_reachability_cycle_onesrc},
             {"CFL_reachability_tree_msrc", test_CFL_reachability_tree_msrc},
             {"CFL_reachability_allin_1", test_CFL_reachability_allin_1},
             {"CFL_reachability_allin_1_4", test_CFL_reachability_allin_1_4},
             {"CFL_reachability_allout_0", test_CFL_reachability_allout_0},
             {"CFL_reachability_allout_1", test_CFL_reachability_allout_1},
            //  Tests with invalid parameters
             {"CFL_reachability_invalid_rules", test_CFL_reachability_invalid_rules},
             {"CFL_reachability_invalid_params", test_CFL_reachability_invalid_params},
             {"CFL_reachability_null_params", test_CFL_reachability_null_params},
             {NULL, NULL}
            };
