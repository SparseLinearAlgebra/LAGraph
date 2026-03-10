#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdarg.h>

#define run_algorithm()                                                                  \
    LAGraph_CFL_AllPaths(outputs, adj_matrices, grammar.terms_count,                 \
                             grammar.nonterms_count, grammar.rules, grammar.rules_count, \
                             msg)

#define check_error(error)                                                               \
    {                                                                                    \
        retval = run_algorithm();                                                        \
        TEST_CHECK(retval == error);                                                     \
        TEST_MSG("retval = %d (%s)", retval, msg);                                       \
    }

#define check_result(nonterm, result) \
    {                                                                                    \
        char *expected = output_to_str(nonterm);                                               \
        TEST_CHECK(strcmp(result, expected) == 0);                                       \
        TEST_MSG("Wrong result. Actual: %s", expected);                                  \
        LAGraph_Free ((void **) &expected, msg);                                         \
    }

typedef struct {
    size_t nonterms_count;
    size_t terms_count;
    size_t rules_count;
    LAGraph_rule_WCNF *rules;
} grammar_t;

GrB_Matrix *adj_matrices = NULL;
int n_adj_matrices = 0 ;
GrB_Matrix *outputs = NULL;
grammar_t grammar = {0, 0, 0, NULL};
char msg[LAGRAPH_MSG_LEN];

void setup() { LAGraph_Init(msg); }

void teardown(void) { LAGraph_Finalize(msg); }

void init_outputs()
{
    LAGraph_Calloc ((void **) &outputs,
        grammar.nonterms_count, sizeof(GrB_Matrix), msg) ;
}

char *output_to_str(size_t nonterm)
{
    GrB_Index nnz = 0;
    OK(GrB_Matrix_nvals(&nnz, outputs[nonterm]));

    if (nnz == 0)
    {
        char *empty = NULL;
        LAGraph_Malloc ((void **) &empty, 1, sizeof(char), msg);
        empty[0] = '\0';
        return empty;
    }

  GrB_Index *row = NULL;
   GrB_Index *col = NULL;
   void *val_void = NULL;

   LAGraph_Malloc ((void **) &row, nnz, sizeof(GrB_Index), msg);
   LAGraph_Malloc ((void **) &col, nnz, sizeof(GrB_Index), msg);
   LAGraph_Malloc ((void **) &val_void, nnz, sizeof(AllPathsElem), msg);

   OK(GrB_Matrix_extractTuples(row, col, val_void, &nnz, outputs[nonterm]));

   AllPathsElem *val = (AllPathsElem *) val_void;

    size_t cap = 256;
    size_t len = 0;
    char *result_str = NULL;
    LAGraph_Malloc ((void **) &result_str, cap, sizeof(char), msg);
    result_str[0] = '\0';

    for (size_t i = 0; i < nnz; i++)
    {
        size_t needed = 64 + (val[i].n) * 20;
        if (len + needed + 1 > cap)
        {
            while (len + needed + 1 > cap) cap *= 2;
            result_str = (char *) realloc(result_str, cap);
            if (result_str == NULL)
            {
                LAGraph_Free ((void **) &row, msg);
                LAGraph_Free ((void **) &col, msg);
                LAGraph_Free ((void **) &val, msg);
                LAGraph_Free ((void **) &result_str, msg);
                OK(GrB_OUT_OF_MEMORY);
            }
        }

        int wrote = 0;
        if (val[i].n == 0 || val[i].middle == NULL)
        {
            wrote = sprintf(result_str + len,
                    (i == 0 ? "(%" PRIu64 ", %" PRIu64 "):[]" : " (%" PRIu64 ", %" PRIu64 "):[]"),
                    (uint64_t) row[i], (uint64_t) col[i]);
            len += (size_t) wrote;
        }
        else
        {
            wrote = sprintf(result_str + len,
                    (i == 0 ? "(%" PRIu64 ", %" PRIu64 "):[" : " (%" PRIu64 ", %" PRIu64 "):["),
                    (uint64_t) row[i], (uint64_t) col[i]);
            len += (size_t) wrote;

            for (GrB_Index k = 0; k < val[i].n; k++)
            {
                if (k == 0)
                {
                  if (val[i].middle[k] == GrB_INDEX_MAX) {
                      wrote = sprintf(result_str + len, "INDEX_MAX");
                  } else {
                      wrote = sprintf(result_str + len, "%" PRIu64, (uint64_t)val[i].middle[k]);
                  }
                }
                else
                {
                  if (val[i].middle[k] == GrB_INDEX_MAX) {
                      wrote = sprintf(result_str + len, ",INDEX_MAX");
                  } else {
                      wrote = sprintf(result_str + len, ",%" PRIu64, (uint64_t)val[i].middle[k]);
                  }
                }
                len += (size_t) wrote;
            }

          result_str[len++] = ']';
            result_str[len] = '\0';
        }
    }

    LAGraph_Free ((void **) &row, msg);
    LAGraph_Free ((void **) &col, msg);
    LAGraph_Free ((void **) &val, msg);

    return result_str;
}

void print_outputs(void)
{
    for (size_t i = 0; i < grammar.nonterms_count; i++)
    {
        char *s = output_to_str(i);
        printf("Output[%zu]: %s\n", i, s[0] ? s : "(empty)");
        LAGraph_Free ((void **) &s, msg);
    }
}

void free_workspace() {

    if (adj_matrices != NULL)
    {
        for (size_t i = 0; i < n_adj_matrices ; i++)
        {
            GrB_free(&adj_matrices[i]);
        }
    }
    LAGraph_Free ((void **) &adj_matrices, msg);

    if (outputs != NULL)
    {
        for (size_t i = 0; i < grammar.nonterms_count; i++)
        {
          void* val_void = NULL;
          GrB_Index nnz = 0;
          GrB_Matrix_nvals(&nnz, outputs[i]);
          LAGraph_Malloc ((void **) &val_void, nnz, sizeof (AllPathsElem), msg) ;
          GrB_Matrix_extractTuples(NULL, NULL, val_void, &nnz, outputs[i]);
          AllPathsElem* val = (AllPathsElem*) val_void;
          for(size_t j = 0; j<nnz; ++j) if(val[j].middle) free(val[j].middle);
          LAGraph_Free (&val_void, msg);
            GrB_free(&outputs[i]);
        }
    }
    LAGraph_Free ((void **) &outputs, msg);

    LAGraph_Free ((void **) &grammar.rules, msg);
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
void init_grammar_aSb() {
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

// S -> aS | a | eps in WCNF
//
// Terms: [0 a]
// Nonterms: [0 S]
// S -> SS [0 0 0 0]
// S -> a  [0 0 -1 0]
// S -> eps [0 -1 -1 0]
void init_grammar_aS() {
    LAGraph_rule_WCNF *rules = NULL ;
    LAGraph_Calloc ((void **) &rules, 3, sizeof(LAGraph_rule_WCNF), msg);

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
void init_grammar_complex() {
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
void init_graph_double_cycle() {
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 1, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 1 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 3, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 3 ;

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
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

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

// 0 -b-> 1
// 1 -b-> 0
void init_graph_4() {
    LAGraph_Calloc ((void **) &adj_matrices, 2, sizeof (GrB_Matrix), msg) ;
    n_adj_matrices = 2 ;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 2, 2);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 2, 2);

    OK(GrB_Matrix_setElement(adj_matrix_b, true, 0, 1));
    OK(GrB_Matrix_setElement(adj_matrix_b, true, 1, 0));

    adj_matrices[0] = adj_matrix_a;
    adj_matrices[1] = adj_matrix_b;
}

//====================
// Tests with valid result
//====================

void test_CFL_AllPaths_cycle(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aS();
    init_graph_one_cycle();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0, "(0, 0):[0,1,2,INDEX_MAX] (0, 1):[0,1,2,INDEX_MAX] (0, 2):[0,1,2] (1, 0):[0,1,2] (1, 1):[0,1,2,INDEX_MAX] (1, 2):[0,1,2,INDEX_MAX] (2, 0):[0,1,2,INDEX_MAX] (2, 1):[0,1,2] (2, 2):[0,1,2,INDEX_MAX]");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_two_cycle(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0, "(0, 0):[1] (0, 3):[1] (1, 0):[2] (1, 3):[2] (2, 0):[0] (2, 3):[0]");
    check_result(1,"(0, 1):[INDEX_MAX] (1, 2):[INDEX_MAX] (2, 0):[INDEX_MAX]");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_labels_more_than_nonterms(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_2();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0, "(0, 1):[0]");
    check_result(1,"(0, 0):[INDEX_MAX]");
    check_result(2,"(0, 1):[INDEX_MAX]");
    check_result(3,"");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_complex_grammar(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_complex();
    init_graph_1();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0, "(0, 7):[1,4] (1, 6):[2]");
    check_result(1, "(0, 4):[2]");
    check_result(6,"(3, 6):[5] (4, 5):[3] (5, 7):[6]");
    check_result(7,"(0, 1):[INDEX_MAX] (1, 2):[INDEX_MAX] (2, 3):[INDEX_MAX] (3, 4):[INDEX_MAX]");
    
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_tree(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_tree();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0,"(0, 0):[2] (0, 1):[2] (0, 3):[2] (0, 4):[2] (1, 0):[2] (1, 1):[2] (1, 3):[2] (1, 4):[2] (2, 2):[6] (2, 5):[6] (3, 0):[5] (3, 1):[5] (3, 3):[5] (3, 4):[5] (4, 0):[5] (4, 1):[5] (4, 3):[5] (4, 4):[5] (5, 2):[6] (5, 5):[6]");
    check_result(1, "(0, 2):[INDEX_MAX] (1, 2):[INDEX_MAX] (2, 6):[INDEX_MAX] (3, 5):[INDEX_MAX] (4, 5):[INDEX_MAX] (5, 6):[INDEX_MAX]");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_line(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_line();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0,"(0, 4):[1] (1, 3):[2]");
    check_result(1,"(0, 1):[INDEX_MAX] (1, 2):[INDEX_MAX]");
    check_result(2,"(2, 3):[INDEX_MAX] (3, 4):[INDEX_MAX]");
    check_result(3,"(1, 4):[3]");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_two_nodes_cycle(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_3();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0,"(0, 0):[1] (1, 0):[0]");
    check_result(1,"(0, 1):[INDEX_MAX] (1, 0):[INDEX_MAX]");
    check_result(2, "(0, 0):[INDEX_MAX]");
    check_result(3, "(0, 0):[0] (1, 0):[0]");
  
    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_with_empty_adj_matrix(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aS();
    init_graph_4();
    init_outputs() ;

    OK(run_algorithm());
    check_result(0, "(0, 0):[0,INDEX_MAX] (1, 1):[1,INDEX_MAX]");

    free_workspace();
    teardown();
#endif
}

//====================
// Tests with invalid result
//====================

void test_CFL_AllPaths_invalid_rules(void) {
#if LAGRAPH_SUITESPARSE
    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs() ;

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

    // Rule [S -> A B], where A >= nonterms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = 10, .prod_B = 2, .index = 0};
    check_error(GrB_INVALID_VALUE);

    // Rule [C -> t], where t >= terms_count
    grammar.rules[0] =
        (LAGraph_rule_WCNF){.nonterm = 0, .prod_A = 10, .prod_B = -1, .index = 0};
    check_error(GrB_INVALID_VALUE);

    free_workspace();
    teardown();
#endif
}

void test_CFL_AllPaths_null_pointers(void) {
#if LAGRAPH_SUITESPARSE

    setup();
    GrB_Info retval;

    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs() ;

//  adj_matrices[0] = NULL;
//  adj_matrices[1] = NULL;
    GrB_free(&adj_matrices[0]);
    GrB_free(&adj_matrices[1]);

    check_error(GrB_NULL_POINTER);

//  adj_matrices = NULL;
    LAGraph_Free ((void **) &adj_matrices, msg);
    check_error(GrB_NULL_POINTER);

    free_workspace();
    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs() ;

//  outputs = NULL;
    LAGraph_Free ((void **) &outputs, msg);
    check_error(GrB_NULL_POINTER);

    free_workspace();
    init_grammar_aSb();
    init_graph_double_cycle();
    init_outputs() ;

//  grammar.rules = NULL;
    LAGraph_Free ((void **) &grammar.rules, msg);
    check_error(GrB_NULL_POINTER);

    free_workspace();
    teardown();
#endif
}

TEST_LIST = {{"CFL_AllPaths_complex_grammar", test_CFL_AllPaths_complex_grammar},
             {"CFL_AllPaths_cycle", test_CFL_AllPaths_cycle},
             {"CFL_AllPaths_two_cycle", test_CFL_AllPaths_two_cycle},
             {"CFL_AllPaths_labels_more_than_nonterms",
              test_CFL_AllPaths_labels_more_than_nonterms},
             {"CFL_AllPaths_tree", test_CFL_AllPaths_tree},
             {"CFL_AllPaths_line", test_CFL_AllPaths_line},
             {"CFL_AllPaths_two_nodes_cycle", test_CFL_AllPaths_two_nodes_cycle},
             {"CFG_reach_basic_invalid_rules", test_CFL_AllPaths_invalid_rules},
             {"test_CFL_AllPaths_with_empty_adj_matrix", test_CFL_AllPaths_with_empty_adj_matrix},
             #if !defined ( GRAPHBLAS_HAS_CUDA )
             {"CFL_AllPaths_null_pointers", test_CFL_AllPaths_null_pointers},
             #endif
             {NULL, NULL}};

