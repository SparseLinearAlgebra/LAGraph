#define LG_FREE_WORK                                                                     \
    do {                                                                                 \
        for (size_t i = 0; i < T_size; i++) {                                            \
            GrB_free(&T[i]);                                                             \
        }                                                                                \
        free(nnz);                                                                       \
    } while (0)

#include "LG_internal.h"
#include <LAGraphX.h>

GrB_Info LAGraph_MTX_reach_basic(GrB_Matrix *output, GrB_Matrix *adj_matrices,
                                 size_t terms_count, size_t nonterms_count,
                                 LAGraph_rule_WCNF *rules, size_t rules_count,
                                 char *msg) {
    // Declare workspace and clear the msg string, if not NULL
    GrB_Matrix T[nonterms_count];
    size_t T_size = 0; // Variable for correct free
    uint64_t *nnz = NULL;
    LG_CLEAR_MSG;

    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "Count of terms must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "Count of nonterms must be greater than zero.")
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "Count of rules must be greater than zero.");

    GrB_Index n;
    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    // Create nonterms matrices
    for (size_t i = 0; i < nonterms_count; i++) {
        GRB_TRY(GrB_Matrix_new(&T[i], GrB_BOOL, n, n));
    }

    // Rule [Variable -> term]
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_WCNF rule = rules[i];
        // Only [Variable -> term] type of rule accepted
        if (rule.prod_A == -1 || rule.prod_B != -1)
            continue;

        GrB_Index nnz;
        GrB_Matrix_nvals(&nnz, adj_matrices[rule.prod_A]);

        GrB_Index row[nnz], col[nnz];
        bool values[nnz];

        GrB_Matrix_extractTuples(row, col, values, &nnz, adj_matrices[rule.prod_A]);

        for (size_t j = 0; j < nnz; j++) {
#ifdef DEBUG
            printf("[TERM] SET ELEMENT [TRUE], NONTERM: %d, ROW: %ld, COL: %ld\n",
                   nonterm, row[j], col[j]);
#endif
            GrB_Matrix_setElement(T[rule.nonterm], 1, row[j], col[j]);
        }
    }

    // Rule [Variable -> eps]
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_WCNF rule = rules[i];

        // Only [Variable -> eps] type of rule accepted
        if (rule.prod_A != -1) {
            continue;
        }

        for (size_t j = 0; j < n; j++) {
#ifdef DEBUG
            printf("[EPS] SET ELEMENT [TRUE], NONTERM: %d, INDEX: %ld\n", nonterm, j);
#endif
            GrB_Matrix_setElement(T[rule.nonterm], true, j, j);
        }
    }

    // Rule [Variable -> Variable1 Variable2]
    nnz = calloc(nonterms_count, sizeof(uint64_t));
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < rules_count; i++) {
            LAGraph_rule_WCNF rule = rules[i];

            // Only [Variable -> Variable1 Variable2] type of rule accepted
            if (rule.prod_A == -1 || rule.prod_B == -1) {
                continue;
            }

            GRB_TRY(GrB_mxm(T[rule.nonterm], GrB_NULL, GxB_LOR_BOOL, GxB_ANY_PAIR_BOOL,
                            T[rule.prod_A], T[rule.prod_B], GrB_NULL));

            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, T[rule.nonterm]));

            changed = changed | (nnz[rule.nonterm] != new_nnz);
#ifdef DEBUG
            printf("[TERM1 TERM2] MULTIPLY, S: %d, A: %d, B: %d, "
                   "I: %ld\n",
                   nonterm, nonterm_A, nonterm_B, i);
#endif
            nnz[rule.nonterm] = new_nnz;
        }
    }

#ifdef DEBUG
    for (size_t i = 0; i < nonterms_count; i++) {
        printf("MATRIX WITH INDEX %ld:\n", i);
        GxB_print(T[i], 5);
    }
#endif

    GrB_Matrix_dup(output, T[0]);
    return GrB_SUCCESS;
}