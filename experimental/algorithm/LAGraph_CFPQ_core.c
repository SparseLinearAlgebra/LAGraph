#define LG_FREE_WORK                                 \
    {                                                \
        LAGraph_Free((void **)&nnzs, NULL);          \
        GrB_free(&false_scalar);                     \
        GrB_free(&identity_matrix);                  \
        LAGraph_Free((void **)&T, NULL);             \
        LAGraph_Free((void **)&indexes, NULL);       \
        LAGraph_Free((void **)&t_empty_flags, NULL); \
        LAGraph_Free((void **)&eps_rules, NULL);     \
        LAGraph_Free((void **)&term_rules, NULL);    \
        LAGraph_Free((void **)&bin_rules, NULL);     \
    }

#define LG_FREE_ALL                                  \
    {                                                \
        for (int64_t i = 0; i < nonterms_count; i++) \
        {                                            \
            GrB_free(&T[i]);                         \
        }                                            \
                                                     \
        LG_FREE_WORK;                                \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

#define ADD_TO_MSG(...)                                                   \
    {                                                                     \
        if (msg_len == 0)                                                 \
        {                                                                 \
            msg_len +=                                                    \
                snprintf(msg, LAGRAPH_MSG_LEN,                            \
                         "LAGraph failure (file %s, line %d): ",          \
                         __FILE__, __LINE__);                             \
        }                                                                 \
        if (msg_len < LAGRAPH_MSG_LEN)                                    \
        {                                                                 \
            msg_len += snprintf(msg + msg_len, LAGRAPH_MSG_LEN - msg_len, \
                                __VA_ARGS__);                             \
        }                                                                 \
    }

// LAGraph_CFPQ_core: Context-Free Path Querying Matrix-Based Algorithm
//
// Internal core function for context-free language path finding.
// Computes path information for all non-terminals using matrix operations defined by the input semiring.
// The semiring determines the specific problem variant being solved.
GrB_Info LAGraph_CFPQ_core(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) contains a corresponding semiring value if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the non-terminal 'k' of the specified CFG.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, terms_count - 1]).
    int64_t terms_count,            // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    const CFL_Semiring *semiring,   // The algebraic structure that defines operations on matrices for a specific problem
    char *msg                       // Message string for error reporting.
)
{

#if LAGRAPH_SUITESPARSE
    // Declare workspace and clear the msg string, if not NULL
    GrB_Matrix *T;
    bool *t_empty_flags = NULL; // t_empty_flags[i] == true <=> T[i] is empty
    GrB_Matrix identity_matrix = NULL;
    uint64_t *nnzs = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    bool iso_flag = false;
    GrB_Index *indexes = NULL;

    // Arrays for processing rules
    size_t *eps_rules = NULL, eps_rules_count = 0;   // [Variable -> eps]
    size_t *term_rules = NULL, term_rules_count = 0; // [Variable -> term]
    size_t *bin_rules = NULL, bin_rules_count = 0;   // [Variable -> AB]

    GrB_Scalar false_scalar;
    GRB_TRY(GrB_Scalar_new(&false_scalar, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(false_scalar, false));

    LG_TRY(LAGraph_Calloc((void **)&T, nonterms_count, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&t_empty_flags, nonterms_count, sizeof(bool), msg));

    LG_TRY(LAGraph_CFL_check_base_inputs(adj_matrices, terms_count, nonterms_count, rules_count, rules, msg, &msg_len));
    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "The outputs array cannot be null.");
    LG_ASSERT_MSG(semiring != NULL, GrB_NULL_POINTER,
                  "The semiring cannot be null.");

    GrB_Index n;
    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    // Create nonterms matrices
    for (int64_t i = 0; i < nonterms_count; i++)
    {
        GRB_TRY(GrB_Matrix_new(&T[i], semiring->type, n, n));
        t_empty_flags[i] = true;
    }

    LG_TRY(LAGraph_Calloc((void **)&eps_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&term_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&bin_rules, rules_count, sizeof(size_t), msg));

    LAGraph_CFL_classify_rules(eps_rules, &eps_rules_count, term_rules, &term_rules_count, bin_rules, &bin_rules_count, rules, rules_count);

    // Rule [Variable -> term]
    for (int64_t i = 0; i < term_rules_count; i++)
    {
        LAGraph_rule_WCNF term_rule = rules[term_rules[i]];
        GrB_Index adj_matrix_nnz = 0;
        GRB_TRY(GrB_Matrix_nvals(&adj_matrix_nnz, adj_matrices[term_rule.prod_A]));

        if (adj_matrix_nnz == 0)
        {
            continue;
        }
        GxB_eWiseUnion(
            T[term_rule.nonterm], GrB_NULL, GrB_NULL, semiring->init_path,
            T[term_rule.nonterm], semiring->bottom_scalar, adj_matrices[term_rule.prod_A], false_scalar, GrB_NULL);

        t_empty_flags[term_rule.nonterm] = false;
    }

    GrB_Vector v_diag;
    GRB_TRY(GrB_Vector_new(&v_diag, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v_diag, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&identity_matrix, v_diag, 0));
    GRB_TRY(GrB_free(&v_diag));

    // Rule [Variable -> eps]
    for (int64_t i = 0; i < eps_rules_count; i++)
    {
        LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];
        GrB_BinaryOp acc_op = t_empty_flags[eps_rule.nonterm] ? GrB_NULL : semiring->add;
        GxB_eWiseUnion(
            T[eps_rule.nonterm], GrB_NULL, acc_op, semiring->init_path,
            T[eps_rule.nonterm], semiring->bottom_scalar, identity_matrix, false_scalar, GrB_NULL);

        t_empty_flags[eps_rule.nonterm] = false;
    }

    // Rule [Variable -> Variable1 Variable2]
    LG_TRY(LAGraph_Calloc((void **)&nnzs, nonterms_count, sizeof(uint64_t), msg));
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (int64_t i = 0; i < bin_rules_count; i++)
        {
            LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];

            // If one of matrices is empty then their product will be empty
            if (t_empty_flags[bin_rule.prod_A] || t_empty_flags[bin_rule.prod_B])
            {
                continue;
            }

            GrB_BinaryOp acc_op = t_empty_flags[bin_rule.nonterm] ? GrB_NULL : semiring->add;
            GRB_TRY(GrB_mxm(T[bin_rule.nonterm], GrB_NULL, acc_op,
                            semiring->semiring, T[bin_rule.prod_A], T[bin_rule.prod_B],
                            GrB_NULL))

            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, T[bin_rule.nonterm]));
            if (new_nnz != 0)
                t_empty_flags[bin_rule.nonterm] = false;

            changed = changed || (nnzs[bin_rule.nonterm] != new_nnz);
            nnzs[bin_rule.nonterm] = new_nnz;
        }
    }

    for (int64_t i = 0; i < nonterms_count; i++)
    {
        outputs[i] = T[i];
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
