#include "LG_internal.h"
#include <LAGraphX.h>

#define LG_FREE_WORK                                 \
    {                                                \
        LAGraph_Free((void **)&T, msg);              \
        GrB_free(&true_scalar);                      \
        GrB_free(&bottom_scalar);                    \
        GrB_free(&IPI_set);                          \
        GrB_free(&IPI_mult);                         \
        GrB_free(&PI_set);                           \
        GrB_free(&PI_mult);                          \
        GrB_free(&PI_add);                           \
        GrB_free(&PI_semiring);                      \
        GrB_free(&Theta);                            \
        GrB_free(&PI_type);                          \
        GrB_free(&PI_monoid);                        \
        GrB_free(&identity_matrix);                  \
        LAGraph_Free((void **)&t_empty_flags, NULL); \
        LAGraph_Free((void **)&nnzs, msg);           \
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

#define ERROR_RULE(msg, i)                                                  \
    {                                                                       \
        LG_ASSERT_MSGF(false, GrB_INVALID_VALUE,                            \
                       "Rule with index %" PRId64 " is invalid. ", msg, i); \
    }

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

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                     \
    {                                                        \
        rule.len_indexes_str += snprintf(                    \
            rule.indexes_str + rule.len_indexes_str,         \
            LAGRAPH_MSG_LEN - rule.len_indexes_str,          \
            rule.count == 0 ? "%" PRId64 : ", %" PRId64, i); \
        rule.count++;                                        \
    }

void add_path_index(PathIndex *z, const PathIndex *x, const PathIndex *y)
{
    // If one is ⊥, then we take another.
    if (x->height == 0 && x->middle == 0)
    {
        *z = *y;
    }
    else if (y->height == 0 && y->middle == 0)
    {
        *z = *x;
    }
    // Take the path with the minimum height
    else if (x->height < y->height || (x->height == y->height && x->middle <= y->middle))
    {
        *z = *x;
    }
    else
    {
        *z = *y;
    }
}

void mult_path_index(PathIndex *z,
                     const PathIndex *x, GrB_Index ix, GrB_Index jx,
                     const PathIndex
                         *y,
                     GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
    // Сoncatenation with ⊥ gives ⊥
    if ((x->height == 0 && x->middle == 0) || (y->height == 0 && y->middle == 0))
    {
        z->height = 0;
        z->middle = 0;
    }
    else
    {
        z->middle = jx;
        z->height = (x->height > y->height ? x->height : y->height) + 1;
    }
}

void set_path_index(PathIndex *z,
                    const PathIndex *x, GrB_Index ix, GrB_Index jx,
                    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,
                    const void *theta)
{
    if (*edge_exist)
    {
        z->middle = ix;
        z->height = 1;
    }
    else
    {
        z->middle = 0;
        z->height = 0;
    }
}

#define SET_PATH_INDEX_DEFN                                                     \
    "void set_path_index(                                                   \n" \
    "    PathIndex *z,                                                      \n" \
    "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n" \
    "    const bool *edge_exist, GrB_Index i_edge, GrB_Index j_edge,        \n" \
    "    const void *theta)                                                 \n" \
    "{                                                                      \n" \
    "    if (*edge_exist)                                                   \n" \
    "    {                                                                  \n" \
    "        z->middle = ix;                                                \n" \
    "        z->height = 1;                                                 \n" \
    "    }                                                                  \n" \
    "    else                                                               \n" \
    "    {                                                                  \n" \
    "        z->middle = 0;                                                 \n" \
    "        z->height = 0;                                                 \n" \
    "    }                                                                  \n" \
    "}"

#define MULT_PATH_INDEX_DEFN                                                     \
    "void mult_path_index(                                                  \n"  \
    "    PathIndex *z,                                                      \n"  \
    "    const PathIndex *x, GrB_Index ix, GrB_Index jx,                    \n"  \
    "    const PathIndex *y, GrB_Index iy, GrB_Index jy,                    \n"  \
    "    const void *theta)                                                 \n"  \
    "{                                                                      \n"  \
    "    if ((x->height == 0 && x->middle == 0) ||                           \n" \
    "        (y->height == 0 && y->middle == 0))                             \n" \
    "    {                                                                  \n"  \
    "        z->height = 0;                                                 \n"  \
    "        z->middle = 0;                                                 \n"  \
    "    }                                                                  \n"  \
    "    else                                                               \n"  \
    "    {                                                                  \n"  \
    "        z->middle = jx;                                                \n"  \
    "        z->height = (x->height > y->height ?                           \n"  \
    "                     x->height : y->height) + 1;                       \n"  \
    "    }                                                                  \n"  \
    "}"

GrB_Info LAGraph_CFL_single_path(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) contains a PathIndex structure if and only if there is a path
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
    char *msg                       // Message string for error reporting.
)
{
    // Declare workspace
    GrB_Matrix *T;
    bool *t_empty_flags = NULL; // t_empty_flags[i] == true <=> T[i] is empty
    uint64_t *nnzs = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    GrB_Type PI_type = NULL;
    GrB_Type_new(&PI_type, sizeof(PathIndex));

    // Create semiring
    GrB_BinaryOp PI_add;
    GrB_BinaryOp_new(
        &PI_add,
        (void *)add_path_index,
        PI_type,
        PI_type,
        PI_type);

    GrB_Monoid PI_monoid;
    PathIndex identity = {0, 0}; // ⊥ - neutral element for the addition operation
    GrB_Monoid_new(
        &PI_monoid,
        PI_add,
        (void *)(&identity));

    GrB_Scalar Theta;
    // Theta cannot be NULL
    GrB_Scalar_new(&Theta, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(Theta, false);

    GxB_IndexBinaryOp IPI_mult;
    GxB_IndexBinaryOp_new(
        &IPI_mult,
        (void *)mult_path_index,
        PI_type,
        PI_type,
        PI_type,
        GrB_BOOL,
        "mult_path_index",
        MULT_PATH_INDEX_DEFN);

    GrB_BinaryOp PI_mult;
    GxB_BinaryOp_new_IndexOp(
        &PI_mult,
        IPI_mult,
        Theta);

    GrB_Semiring PI_semiring;
    GrB_Semiring_new(
        &PI_semiring,
        PI_monoid,
        PI_mult);

    GxB_IndexBinaryOp IPI_set;
    GxB_IndexBinaryOp_new(
        &IPI_set,
        (void *)set_path_index,
        PI_type,
        PI_type,
        GrB_BOOL,
        GrB_BOOL,
        "set_path_index",
        SET_PATH_INDEX_DEFN);

    GrB_BinaryOp PI_set;
    GxB_BinaryOp_new_IndexOp(
        &PI_set,
        IPI_set,
        Theta);

    GrB_Matrix identity_matrix = NULL;
    GrB_Scalar true_scalar, bottom_scalar;
    PathIndex bottom = {0, 0};
    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);
    GrB_Scalar_new(&bottom_scalar, PI_type);
    GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom));

    // Arrays for processing rules
    size_t *eps_rules = NULL, eps_rules_count = 0;   // [Variable -> eps]
    size_t *term_rules = NULL, term_rules_count = 0; // [Variable -> term]
    size_t *bin_rules = NULL, bin_rules_count = 0;   // [Variable -> AB]

    LG_TRY(LAGraph_Calloc((void **)&T, nonterms_count, sizeof(GrB_Matrix), msg));
    LG_TRY(LAGraph_Calloc((void **)&t_empty_flags, nonterms_count, sizeof(bool), msg));

    LG_ASSERT_MSG(terms_count > 0, GrB_INVALID_VALUE,
                  "The number of terminals must be greater than zero.");
    LG_ASSERT_MSG(nonterms_count > 0, GrB_INVALID_VALUE,
                  "The number of non-terminals must be greater than zero.");
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "The number of rules must be greater than zero.");
    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "The outputs array cannot be null.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "The rules array cannot be null.");
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "The adjacency matrices array cannot be null.");

    // Find null adjacency matrices
    bool found_null = false;
    for (int64_t i = 0; i < terms_count; i++)
    {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null)
        {
            ADD_TO_MSG("Adjacency matrices with these indexes are null: ");
            ADD_TO_MSG("%" PRId64, i);
        }
        else
        {
            ADD_TO_MSG("%" PRId64, i);
        }

        found_null = true;
    }

    if (found_null)
    {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }
    GrB_Index n;
    GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    // Create nonterms matrices
    for (int64_t i = 0; i < nonterms_count; i++)
    {
        GRB_TRY(GrB_Matrix_new(&T[i], PI_type, n, n));
        t_empty_flags[i] = true;
    }
    LG_TRY(LAGraph_Calloc((void **)&eps_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&term_rules, rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&bin_rules, rules_count, sizeof(size_t), msg));

    // Process rules
    typedef struct
    {
        size_t count;
        size_t len_indexes_str;
        char indexes_str[LAGRAPH_MSG_LEN];
    } rule_error_s;
    rule_error_s term_err = {0};
    rule_error_s nonterm_err = {0};
    rule_error_s invalid_err = {0};
    for (int64_t i = 0; i < rules_count; i++)
    {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // Check that all rules are well-formed
        if (rule.nonterm < 0 || rule.nonterm >= nonterms_count)
        {
            ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
        }

        // [Variable -> eps]
        if (is_rule_eps)
        {
            eps_rules[eps_rules_count++] = i;

            continue;
        }

        // [Variable -> term]
        if (is_rule_term)
        {
            term_rules[term_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= terms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin)
        {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || rule.prod_A >= nonterms_count || rule.prod_B < -1 ||
                rule.prod_B >= nonterms_count)
            {
                ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
            }

            continue;
        }

        // [Variable -> _ B]
        ADD_INDEX_TO_ERROR_RULE(invalid_err, i);
    }
    if (term_err.count + nonterm_err.count + invalid_err.count > 0)
    {
        ADD_TO_MSG("Count of invalid rules: %" PRId64 ".\n",
                   (int64_t)(term_err.count + nonterm_err.count + invalid_err.count));

        if (nonterm_err.count > 0)
        {
            ADD_TO_MSG("Non-terminals must be in range [0, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indexes_str);
        }
        if (term_err.count > 0)
        {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indexes_str);
        }
        if (invalid_err.count > 0)
        {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indexes_str);
        }

        LG_FREE_ALL;
        return GrB_INVALID_VALUE;
    }

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
        GrB_BinaryOp acc_op = t_empty_flags[term_rule.nonterm] ? GrB_NULL : PI_add;
        GxB_eWiseUnion(
            T[term_rule.nonterm], GrB_NULL, acc_op, PI_set,
            T[term_rule.nonterm], bottom_scalar, adj_matrices[term_rule.prod_A], true_scalar, GrB_NULL);

        t_empty_flags[term_rule.nonterm] = false;
    }

    // Rule [Variable -> eps]
    GrB_Vector v_diag;
    GRB_TRY(GrB_Vector_new(&v_diag, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v_diag, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&identity_matrix, v_diag, 0));
    GRB_TRY(GrB_free(&v_diag));

    for (int64_t i = 0; i < eps_rules_count; i++)
    {
        LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];
        GrB_BinaryOp acc_op = t_empty_flags[eps_rule.nonterm] ? GrB_NULL : PI_add;
        GxB_eWiseUnion(
            T[eps_rule.nonterm], GrB_NULL, acc_op, PI_set,
            T[eps_rule.nonterm], bottom_scalar, identity_matrix, true_scalar, GrB_NULL);

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

            GrB_BinaryOp acc_op = t_empty_flags[bin_rule.nonterm] ? GrB_NULL : PI_add;
            GRB_TRY(GrB_mxm(T[bin_rule.nonterm], GrB_NULL, acc_op,
                            PI_semiring, T[bin_rule.prod_A], T[bin_rule.prod_B],
                            GrB_NULL))

            GrB_Index new_nnz;
            GRB_TRY(GrB_Matrix_nvals(&new_nnz, T[bin_rule.nonterm]));
            if (new_nnz != 0)
                t_empty_flags[bin_rule.nonterm] = false;

            changed = changed || (nnzs[bin_rule.nonterm] != new_nnz);
            nnzs[bin_rule.nonterm] = new_nnz;
        }
    }
    // finish
    for (int64_t i = 0; i < nonterms_count; i++)
    {
        outputs[i] = T[i];
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
