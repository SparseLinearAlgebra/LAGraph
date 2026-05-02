#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        if (outputs_reachability != NULL) {                                              \
            for (size_t i = 0; i < symbols_amount; i++) {                                \
                GrB_free(&outputs_reachability[i]);                                      \
            }                                                                            \
        }                                                                                \
        LAGraph_Free((void **)&outputs_reachability, NULL);                              \
        LAGraph_Free((void **)&T, NULL);                                                 \
        GrB_free(&false_scalar);                                                         \
        GrB_free(&identity_matrix);                                                      \
        GrB_free(&v_diag);                                                               \
        LAGraph_Free((void **)&t_empty_flags, NULL);                                     \
        LAGraph_Free((void **)&eps_rules, NULL);                                         \
        LAGraph_Free((void **)&term_rules, NULL);                                        \
        LAGraph_Free((void **)&bin_rules, NULL);                                         \
        LAGraph_Free((void **)&new_rules, NULL);                                         \
        GrB_free(&AllPaths_semiring);                                                    \
        GrB_free(&AllPaths_monoid);                                                      \
        GrB_free(&bottom_scalar);                                                        \
        GrB_free(&IAllPaths_mult);                                                       \
        GrB_free(&AllPaths_set);                                                         \
        GrB_free(&AllPaths_mult);                                                        \
        GrB_free(&AllPaths_add);                                                         \
        GrB_free(&Theta);                                                                \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

// Merging two ordered arrays of internal vertices in the add function
static GrB_Index *merge_all_paths(size_t *n, const GrB_Index *a, const size_t na,
                                  const GrB_Index *b, const size_t nb) {
    GrB_Index *tmp = malloc((na + nb) * sizeof(GrB_Index));

    size_t ia = 0, ib = 0, outn = 0;
    // Handle the first elements of both arrays to initialize tmp and avoid checking if
    // tmp is empty in the loop
    if (na > 0 && nb > 0) {
        if (a[0] < b[0]) {
            tmp[outn++] = a[ia++];
        } else if (b[0] < a[0]) {
            tmp[outn++] = b[ib++];
        } else {
            tmp[outn++] = a[ia++];
            ib++;
        }
    } else if (na > 0) {
        tmp[outn++] = a[ia++];
    } else {
        tmp[outn++] = b[ib++];
    }

    while (ia < na && ib < nb) {
        GrB_Index va = a[ia];
        GrB_Index vb = b[ib];
        if (va < vb) {
            if (tmp[outn - 1] != va)
                tmp[outn++] = va;
            ia++;
        } else if (vb < va) {
            if (tmp[outn - 1] != vb)
                tmp[outn++] = vb;
            ib++;
        } else {
            if (tmp[outn - 1] != va)
                tmp[outn++] = va;
            ia++;
            ib++;
        }
    }
    while (ia < na) {
        GrB_Index va = a[ia++];
        if (tmp[outn - 1] != va)
            tmp[outn++] = va;
    }
    while (ib < nb) {
        GrB_Index vb = b[ib++];
        if (tmp[outn - 1] != vb)
            tmp[outn++] = vb;
    }

    tmp = realloc(tmp, outn * sizeof(GrB_Index));
    *n = outn;

    return tmp;
}

static GrB_Index *insert_all_paths(size_t *n, const GrB_Index *arr, size_t len,
                                   GrB_Index value) {
    // Array is ordered, so we can use binary search to find the position to insert value
    size_t l = 0, r = len, m = 0;
    while (l < r) {
        m = l + (r - l) / 2;
        if (arr[m] < value)
            l = m + 1;
        else
            r = m;
    }

    // If value is already in the array, return the original array
    if (l < len && arr[l] == value) {
        GrB_Index *tmp = malloc(len * sizeof(GrB_Index));
        memcpy(tmp, arr, len * sizeof(GrB_Index));
        *n = len;
        return tmp;
    }

    GrB_Index *tmp = malloc((len + 1) * sizeof(GrB_Index));
    memcpy(tmp, arr, l * sizeof(GrB_Index));
    tmp[l] = value;
    memcpy(tmp + l + 1, arr + l, (len - l) * sizeof(GrB_Index));

    *n = len + 1;
    return tmp;
}

static void add_all_paths(AllPathsElem *z, AllPathsElem *x, AllPathsElem *y) {
    // temp is needed to avoid freeing the memory of z in case z == x or z == y
    AllPathsElem temp;

    // x->n and y->n are not equal to zero
    if (x->n == 1 && y->n == 1) {
        if (x->data.single_elem == y->data.single_elem) {
            temp.n = 1;
            temp.data.single_elem = x->data.single_elem;
        } else {
            temp.n = 2;
            temp.data.middle = malloc(2 * sizeof(GrB_Index));
            if (x->data.single_elem < y->data.single_elem) {
                temp.data.middle[0] = x->data.single_elem;
                temp.data.middle[1] = y->data.single_elem;
            } else {
                temp.data.middle[0] = y->data.single_elem;
                temp.data.middle[1] = x->data.single_elem;
            }
        }
    } else if (x->n == 1) {
        temp.data.middle =
            insert_all_paths(&temp.n, y->data.middle, y->n, x->data.single_elem);
    } else if (y->n == 1) {
        temp.data.middle =
            insert_all_paths(&temp.n, x->data.middle, x->n, y->data.single_elem);
    } else {
        temp.data.middle =
            merge_all_paths(&temp.n, x->data.middle, x->n, y->data.middle, y->n);
    }

    if (x->n > 1) {
        free(x->data.middle);
    }
    if (y->n > 1) {
        free(y->data.middle);
    }

    *z = temp;
}

static void mult_all_paths_post(AllPathsElem *z, const void *x, GrB_Index ix,
                                GrB_Index jx, const void *y, GrB_Index iy, GrB_Index jy,
                                const void *theta) {
    z->data.single_elem = jx;
    z->n = 1;
}

static void set_all_paths(AllPathsElem *z, const AllPathsElem *x,
                          const bool *edge_exist) {
    z->data.single_elem =
        GrB_INDEX_MAX; // A special value to indicate that this path corresponds to a
                       // terminal rule (A->t) or an epsilon rule (A->eps)
    z->n = 1;
}

#define MULT_PATH_POST_INDEX_DEFN                                                        \
    "static void mult_all_paths_post(AllPathsElem *z, \n"                                \
    "                     const void *x, GrB_Index ix, GrB_Index jx, \n"                 \
    "                     const void *y, GrB_Index iy, GrB_Index jy, \n"                 \
    "                     const void *theta) \n"                                         \
    "{ \n"                                                                               \
    "  z->data.single_elem = jx; \n"                                                     \
    "  z->n = 1; \n"                                                                     \
    "}"

GrB_Info LAGraph_CFL_AllPaths_adv(GrB_Matrix *outputs, GrB_Type *all_paths_ptr_t,
                                  const GrB_Matrix *adj_matrices, size_t symbols_amount,
                                  const LAGraph_rule_EWCNF *rules, size_t rules_count,
                                  char *msg, int8_t optimizations) {
    LG_CLEAR_MSG;
    size_t msg_len = 0;

    GrB_Type AllPaths_type = NULL;

    GrB_Matrix *T = NULL;
    GrB_Scalar false_scalar = NULL;
    GrB_Matrix *outputs_reachability = NULL;
    bool *t_empty_flags = NULL;
    LAGraph_rule_EWCNF *new_rules = NULL;
    size_t new_rules_count = 0;
    size_t *eps_rules = NULL, eps_rules_count = 0;
    size_t *term_rules = NULL, term_rules_count = 0;
    size_t *bin_rules = NULL, bin_rules_count = 0;
    GrB_Matrix identity_matrix = NULL;
    GrB_Vector v_diag = NULL;

#if GxB_IMPLEMENTATION < GxB_VERSION(9, 4, 5)
    return (GrB_NOT_IMPLEMENTED);
#else
    GrB_BinaryOp AllPaths_add = NULL;
    GrB_Monoid AllPaths_monoid = NULL;
    GxB_IndexBinaryOp IAllPaths_mult = NULL;
    GrB_BinaryOp AllPaths_mult = NULL;
    GrB_Semiring AllPaths_semiring = NULL;
    GrB_BinaryOp AllPaths_set = NULL;
    GrB_Scalar Theta = NULL;
    GrB_Scalar bottom_scalar = NULL;

    GrB_free(all_paths_ptr_t);
    GRB_TRY(GxB_Type_new(all_paths_ptr_t, sizeof(AllPathsElem), "AllPathsElem",
                         "typedef struct{size_t n;union{GrB_Index single_elem;GrB_Index* "
                         "middle;}data;}AllPathsElem;"));
    AllPaths_type = *all_paths_ptr_t;

    GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, false));

    AllPathsElem bottom = {0};
    GRB_TRY(GrB_Scalar_new(&bottom_scalar, AllPaths_type));
    GRB_TRY(GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom)));

    GRB_TRY(GrB_BinaryOp_new(&AllPaths_add, (void *)add_all_paths, AllPaths_type,
                             AllPaths_type, AllPaths_type));
    GRB_TRY(GrB_Monoid_new(&AllPaths_monoid, AllPaths_add, (void *)(&bottom)));

    GRB_TRY(GxB_IndexBinaryOp_new(&IAllPaths_mult, (void *)mult_all_paths_post,
                                  AllPaths_type, GrB_BOOL, GrB_BOOL, GrB_BOOL,
                                  "mult_all_paths_post", MULT_PATH_POST_INDEX_DEFN));

    GRB_TRY(GxB_BinaryOp_new_IndexOp(&AllPaths_mult, IAllPaths_mult, Theta));
    GRB_TRY(GrB_Semiring_new(&AllPaths_semiring, AllPaths_monoid, AllPaths_mult));
    GRB_TRY(GrB_BinaryOp_new(&AllPaths_set, (void *)set_all_paths, AllPaths_type,
                             AllPaths_type, GrB_BOOL));

    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "The outputs array cannot be null.");
    LG_CLEAR_MSG;
    msg_len = 0;

    LAGraph_Calloc((void **)&outputs_reachability, symbols_amount, sizeof(GrB_Matrix),
                   msg);
    LG_TRY(LAGraph_CFL_reachability_adv(outputs_reachability, adj_matrices,
                                        symbols_amount, rules, rules_count, msg,
                                        optimizations));

    for (size_t i = 0; i < symbols_amount; i++) {
        GrB_set(outputs_reachability[i], GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT);
    }

    for (size_t i = 0; i < rules_count; i++) {
        new_rules_count += (rules[i].indexed_count == 0) ? 1 : rules[i].indexed_count;
    }

    LG_TRY(LAGraph_Calloc((void **)&new_rules, new_rules_count,
                          sizeof(LAGraph_rule_EWCNF), msg));
    size_t nr_idx = 0;
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_EWCNF rule = rules[i];
        if (rule.indexed_count == 0) {
            new_rules[nr_idx++] = rule;
        } else {
            for (size_t rule_index = 0; rule_index < rule.indexed_count; rule_index++) {
                LAGraph_rule_EWCNF exploded_rule = rule;
                exploded_rule.indexed_count = 0;
                exploded_rule.indexed = 0;
                if (rule.nonterm != -1 && rule.indexed & LAGraph_EWNCF_INDEX_NONTERM) {
                    exploded_rule.nonterm = rule.nonterm + rule_index;
                }
                if (rule.prod_A != -1 && rule.indexed & LAGraph_EWNCF_INDEX_PROD_A) {
                    exploded_rule.prod_A = rule.prod_A + rule_index;
                }
                if (rule.prod_B != -1 && rule.indexed & LAGraph_EWNCF_INDEX_PROD_B) {
                    exploded_rule.prod_B = rule.prod_B + rule_index;
                }
                new_rules[nr_idx++] = exploded_rule;
            }
        }
    }

    LG_TRY(LAGraph_Calloc((void **)&T, symbols_amount, sizeof(GrB_Matrix), msg));
    GRB_TRY(GrB_Scalar_new(&false_scalar, GrB_BOOL));
    GRB_TRY(GrB_Scalar_setElement_BOOL(false_scalar, false));
    LG_TRY(LAGraph_Calloc((void **)&t_empty_flags, symbols_amount, sizeof(bool), msg));

    GrB_Index n;
    bool found_n = false;
    for (size_t i = 0; i < symbols_amount; i++) {
        if (adj_matrices[i] != NULL) {
            GRB_TRY(GrB_Matrix_ncols(&n, adj_matrices[i]));
            found_n = true;
            break;
        }
    }
    if (!found_n) {
        n = 0;
    }

    for (size_t i = 0; i < symbols_amount; i++) {
        GRB_TRY(GrB_Matrix_new(&T[i], AllPaths_type, n, n));
        t_empty_flags[i] = true;
    }

    LG_TRY(LAGraph_Calloc((void **)&eps_rules, new_rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&term_rules, new_rules_count, sizeof(size_t), msg));
    LG_TRY(LAGraph_Calloc((void **)&bin_rules, new_rules_count, sizeof(size_t), msg));

    for (size_t i = 0; i < new_rules_count; i++) {
        LAGraph_rule_EWCNF rule = new_rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        if (is_rule_eps) {
            eps_rules[eps_rules_count++] = i;
            continue;
        }
        if (is_rule_term) {
            term_rules[term_rules_count++] = i;
            continue;
        }
        if (is_rule_bin) {
            bin_rules[bin_rules_count++] = i;
            continue;
        }
    }

    for (size_t i = 0; i < term_rules_count; i++) {
        LAGraph_rule_EWCNF term_rule = new_rules[term_rules[i]];
        if (adj_matrices[term_rule.prod_A] == NULL)
            continue;

        GrB_Index adj_matrix_nnz = 0;
        GRB_TRY(GrB_Matrix_nvals(&adj_matrix_nnz, adj_matrices[term_rule.prod_A]));
        if (adj_matrix_nnz == 0)
            continue;

        GxB_eWiseUnion(T[term_rule.nonterm], GrB_NULL, GrB_NULL, AllPaths_set,
                       T[term_rule.nonterm], bottom_scalar,
                       adj_matrices[term_rule.prod_A], false_scalar, GrB_NULL);
        t_empty_flags[term_rule.nonterm] = false;
    }

    if (eps_rules_count > 0) {
        GRB_TRY(GrB_Vector_new(&v_diag, GrB_BOOL, n));
        GRB_TRY(
            GrB_Vector_assign_BOOL(v_diag, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
        GRB_TRY(GrB_Matrix_diag(&identity_matrix, v_diag, 0));
        GRB_TRY(GrB_free(&v_diag));

        for (size_t i = 0; i < eps_rules_count; i++) {
            LAGraph_rule_EWCNF eps_rule = new_rules[eps_rules[i]];
            GrB_BinaryOp acc_op =
                t_empty_flags[eps_rule.nonterm] ? GrB_NULL : AllPaths_add;
            GxB_eWiseUnion(T[eps_rule.nonterm], GrB_NULL, acc_op, AllPaths_set,
                           T[eps_rule.nonterm], bottom_scalar, identity_matrix,
                           false_scalar, GrB_NULL);
            t_empty_flags[eps_rule.nonterm] = false;
        }
        GrB_free(&identity_matrix);
    }

    for (size_t i = 0; i < symbols_amount; i++) {
        GrB_Index temp_nvals = 0;
        GrB_Matrix_nvals(&temp_nvals, outputs_reachability[i]);
        if (temp_nvals != 0) {
            t_empty_flags[i] = false;
        }
    }

    for (size_t i = 0; i < bin_rules_count; i++) {
        LAGraph_rule_EWCNF bin_rule = new_rules[bin_rules[i]];

        if (t_empty_flags[bin_rule.prod_A] || t_empty_flags[bin_rule.prod_B])
            continue;

        GrB_BinaryOp acc_op = t_empty_flags[bin_rule.nonterm] ? GrB_NULL : AllPaths_add;
        GRB_TRY(GrB_mxm(T[bin_rule.nonterm], GrB_NULL, acc_op, AllPaths_semiring,
                        outputs_reachability[bin_rule.prod_A],
                        outputs_reachability[bin_rule.prod_B], GrB_NULL));
    }

    for (size_t i = 0; i < symbols_amount; i++) {
        outputs[i] = T[i];
        GrB_Matrix_wait(outputs[i], GrB_MATERIALIZE);
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
#endif
}

// Helper function to free the output matrix of LAGraph_CFL_AllPaths, which contains
// elements of type AllPathsElem with dynamically allocated arrays of intermediate
// vertices.
static void free_AllPaths_matrix(GrB_Matrix *ptr_output) {
    GxB_Iterator iterator;
    GxB_Iterator_new(&iterator);
    GrB_Info info = GxB_Matrix_Iterator_attach(iterator, *ptr_output, NULL);
    info = GxB_Matrix_Iterator_seek(iterator, 0);
    AllPathsElem val;

    while (info != GxB_EXHAUSTED) {
        GxB_Iterator_get_UDT(iterator, (void *)&val);
        if (val.n > 1 && val.data.middle != NULL) {
            free(val.data.middle);
        }
        info = GxB_Matrix_Iterator_next(iterator);
    }

    GrB_free(&iterator);
    GrB_free(ptr_output);
}

// Free outputs and all_paths_ptr_t after you have finished working with the output
// matrices from LAGraph_CFL_AllPaths. do outputs = NULL, all_paths_ptr_t = NULL after
// LAGraph_CFL_AllPaths_free_outputs
GrB_Info LAGraph_CFL_AllPaths_adv_free_outputs(GrB_Matrix *outputs,
                                               int64_t nonterms_count,
                                               GrB_Type *all_paths_ptr_t) {
#if GxB_IMPLEMENTATION < GxB_VERSION(9, 4, 5)
    return (GrB_NOT_IMPLEMENTED);
#else
    if (outputs) {
        for (size_t i = 0; i < nonterms_count; i++) {
            if (outputs[i] == NULL)
                continue;
            free_AllPaths_matrix(&outputs[i]);
            outputs[i] = NULL;
        }
        free(outputs);
    }
    GrB_free(all_paths_ptr_t);
    return GrB_SUCCESS;
#endif
}
