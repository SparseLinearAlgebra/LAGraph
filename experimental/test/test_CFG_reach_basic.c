#include "LG_internal.h"
#include <LAGraphX.h>

int main() {

    char *msg = malloc(256);
    LAGRAPH_TRY(LAGraph_Init(msg));

    GrB_Matrix output;

    GrB_Matrix adj_matrix_a, adj_matrix_b;
    GrB_Matrix_new(&adj_matrix_a, GrB_BOOL, 4, 4);
    GrB_Matrix_new(&adj_matrix_b, GrB_BOOL, 4, 4);
    GRB_TRY(GrB_Matrix_setElement(adj_matrix_a, true, 0, 1));
    GRB_TRY(GrB_Matrix_setElement(adj_matrix_a, true, 1, 2));
    GRB_TRY(GrB_Matrix_setElement(adj_matrix_a, true, 2, 0));
    GRB_TRY(GrB_Matrix_setElement(adj_matrix_b, true, 0, 3));
    GRB_TRY(GrB_Matrix_setElement(adj_matrix_b, true, 3, 0));
    GrB_Matrix adj_matrices[2] = {adj_matrix_a, adj_matrix_b};

    // Rule defined by [NONTERM, PROD_A, PROD_B, INDEX] in Chomsky Normal Form
    // Variable -> EPS: [NONTERM, -1, -1, INDEX]
    // Variable -> term: [NONTERM, TERM, -1, INDEX]
    // Variable -> AB: [NONTERM, TERM1, TERM2, INDEX]
    // Terms: [0 a] [1 b]
    // Nonterms: [0 S] [1 A] [2 B] [3 C]
    // S -> AB [0 1 2 0]
    // S -> AC [0 1 3 0]
    // C -> SB [3 0 2 0]
    // A -> a  [1 0 -1 0]
    // B -> b  [2 1 -1 0]
    LAGraph_rule_WCNF rules[5] = {
        {0, 1, 2, 0}, {0, 1, 3, 0}, {3, 0, 2, 0}, {1, 0, -1, 0}, {2, 1, -1, 0}};

    LAGraph_CFG_reach_basic(&output, adj_matrices, 2, 6, rules,
                            sizeof(rules) / (sizeof(LAGraph_rule_WCNF)), msg);

    GrB_Index nnz = 0;
    GRB_TRY(GrB_Matrix_nvals(&nnz, output));
    GrB_Index *row = malloc(nnz * sizeof(GrB_Index));
    GrB_Index *col = malloc(nnz * sizeof(GrB_Index));
    bool *val = malloc(nnz * sizeof(GrB_Index));
    GRB_TRY(GrB_Matrix_extractTuples(row, col, val, &nnz, output));

    printf("RESULT:\n");
    for (size_t i = 0; i < nnz; i++) {
        printf(i == 0 ? "(%ld, %ld)" : " (%ld, %ld)", row[i], col[i]);
    }
    printf("\n");

    return 0;
}