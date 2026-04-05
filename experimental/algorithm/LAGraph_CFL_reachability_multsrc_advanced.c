//------------------------------------------------------------------------------
// LAGraph_CFL_reachability_multsrc.c: Optimized Multiple-Source Context-Free
// Language Reachability Matrix-Based Algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Ilhom Kombaev, Vlasenco Daniel, Semyon Grigoriev, St. Petersburg State University.

//------------------------------------------------------------------------------

// Code is based on the "multiple-source CFPQ", described in the following paper:
//  * Arseny Terekhov et al., "Multiple-Source Context-Free Path Querying in Terms of Linear Algebra"
//  * URL: https://openproceedings.org/2021/conf/edbt/p48.pdf

// #define DEBUG_CFL_REACHABILITY

#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        for (size_t i = 0; i < new_symbols_amount; i++) {                                    \
            CFL_matrix_free(&T[i]);                                                      \
            CFL_matrix_free(&dT[i]);                                                     \
            CFL_matrix_free(&TSrc[i]);                                                   \
        }                                                                                \
        CFL_matrix_free(&iden);                                                          \
        LAGraph_Free((void **) &nnzs_T, msg);                                            \
        LAGraph_Free((void **) &nnzs_TSrc_B, msg);                                       \
        LAGraph_Free((void **) &nnzs_TSrc_C, msg);                                       \
        LAGraph_Free((void **) &ones_vec, msg);                                          \
        LAGraph_Free((void **) &T, msg);                                                 \
        LAGraph_Free((void **) &TSrc, msg);                                              \
        LAGraph_Free((void **) &MSrc, msg);                                              \
        LAGraph_Free ((void **) &A, msg);                                                \
        LAGraph_Free ((void **) &a, msg);                                                \
        GrB_free(&true_scalar);                                                          \
    }

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        LG_FREE_WORK;                                                                    \
    }

#include "LG_internal.h"
#include "LAGraph_CFL_optimized_matrix_opt.h"
#include <LAGraphX.h>

#define ERROR_RULE(msg)                                                                  \
    {                                                                                    \
        LG_ASSERT_MSGF(false, GrB_INVALID_VALUE, "Rule with index %ld is invalid. " msg, \
                       i);                                                               \
    }

#define ADD_TO_MSG(...)                                                                  \
    {                                                                                    \
        if (msg_len == 0) {                                                              \
            msg_len +=                                                                   \
                snprintf(msg, LAGRAPH_MSG_LEN,                                           \
                         "LAGraph failure (file %s, line %d): ", __FILE__, __LINE__);    \
        }                                                                                \
        if (msg_len < LAGRAPH_MSG_LEN) {                                                 \
            msg_len += snprintf(msg + msg_len, LAGRAPH_MSG_LEN - msg_len, __VA_ARGS__);  \
        }                                                                                \
    }

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                                                \
    {                                                                                   \
        rule.len_indices_str += snprintf(rule.indices_str + rule.len_indices_str,       \
        LAGRAPH_MSG_LEN - rule.len_indices_str,                                         \
                                         rule.count == 0 ? "%ld" : ", %ld", i);         \
        rule.count++;                                                                   \
    }

#define PRINT_MATRIX(_m) {                                                              \
    for (size_t _i = 0; _i < n; _i++) {                                                 \
        for (size_t _j = 0; _j < n; _j++) {                                             \
            if (GxB_Matrix_isStoredElement(_m, _i, _j) == GrB_SUCCESS) {                \
                printf("1 ");                                                           \
            }                                                                           \
            else {                                                                      \
                printf("0 ");                                                           \
            }                                                                           \
        }                                                                               \
        printf("\n");                                                                   \
    }                                                                                   \
}

#define PRINT_VECTOR(_v) {                                                              \
    for (size_t _i = 0; _i < n; _i++) {                                                 \
        if (GxB_Vector_isStoredElement(_v, _i) == GrB_SUCCESS) {                        \
            printf("1 ");                                                               \
        }                                                                               \
        else {                                                                          \
            printf("0 ");                                                               \
        }                                                                               \
    }                                                                                   \
    printf("\n");                                                                       \
}

#define PRINT_RULE(_r) {                                                                \
    printf("%c -> ", _r.nonterm + 'A' - 1);                                             \
    if (_r.prod_A == -1 && _r.prod_B == -1) {                                           \
        printf("eps\n");                                                                \
    } else if (_r.prod_A != -1 && _r.prod_B == -1) {                                    \
        printf("%c\n", _r.prod_A + 'a');                                                \
    } else {                                                                            \
        printf("%c %c\n", _r.prod_A + 'A'- 1 ,_r.prod_B + 'A' - 1);                     \
    }                                                                                   \
}

#define TRY(GrB_method)                                                                  \
    {                                                                                    \
        GrB_Info LG_GrB_Info = GrB_method;                                               \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                 \
            fprintf(stderr, "LAGraph failure (file %s, line %d): \n", __FILE__,          \
                    __LINE__);                                                           \
            LG_FREE_ALL;                                                                 \
            return (LG_GrB_Info);                                                        \
        }                                                                                \
    }

// Checks the return value of a GraphBLAS/LAGraph call inside a helper function.
// On failure, logs the error location to stderr, invokes FREE_INNER() to release
// any resources allocated within the current function, and returns the error code.
//
// FREE_INNER() must be defined by the caller before using this macro:
//
//   #define FREE_INNER()       \
//       {                      \
//           LAGraph_Free(&a)   \
//           GrB_free(&M)       \
//       }
//
// Note: use TRY() instead when the function is not a helper (i.e., it has
// its own top-level FREE_ALL cleanup).
#define TRY_INNER(GrB_method)                                                            \
    {                                                                                    \
        GrB_Info LG_GrB_Info = GrB_method;                                               \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                 \
            fprintf(stderr, "LAGraph failure (file %s, line %d): \n", __FILE__,          \
                    __LINE__);                                                           \
            FREE_INNER();                                                                \
            return (LG_GrB_Info);                                                        \
        }                                                                                \
    }


#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

typedef struct {
    int32_t index;
    int32_t base_index;
    int32_t count;
} CFL_Symbol;

// When using the OPT_BLOCK optimization, indexed symbols must be grouped together.
// This produces a mapping: [old_index -> (new_index, base_index, indexed_count)]
//   - new_index:     the index of the symbol in the new numeration
//   - base_index:    the original index of the first symbol in the indexed group
//   - indexed_count: the number of indexed symbols in the group (0 if not indexed)
//
// Example:
//   (0) S   -> (0, 0, 0)  - non-indexed symbol
//   (1) A_0 -> (1, 1, 3)  - indexed group of 3, starting at old index 1
//   (2) A_1 ->              }
//   (3) A_2 ->              } members of the A group
//   (4) B_0 -> (2, 4, 2)  - indexed group of 2, starting at old index 4
//   (5) B_1 ->             } member of the B group
//   (6) C   -> (3, 6, 0)  - non-indexed symbol
//   (7) a   -> (4, 7, 0)  - non-indexed symbol
//
// This mapping is used to build a compact matrix array and to expand
// production rules, both with and without the OPT_BLOCK optimization.
//
// Output: CFL_Symbol **symbols and size_t *size
static GrB_Info get_new_symbols(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                                size_t symbols_amount, CFL_Symbol **symbols, size_t *size,
                                char *msg) {
    *symbols = NULL;
    bool *checked = NULL;

#undef FREE_INNER_WORK
#undef FREE_INNER

#define FREE_INNER_WORK()                                                                \
    {                                                                                    \
        LAGraph_Free((void **)&checked, msg);                                            \
    };

#define FREE_INNER()                                                                     \
    {                                                                                    \
        FREE_INNER_WORK();                                                               \
        LAGraph_Free((void **)symbols, msg);                                             \
    }

    TRY_INNER(LAGraph_Calloc((void **)&checked, symbols_amount, sizeof(bool), msg));
    for (size_t i = 0; i < symbols_amount; i++) {
        checked[i] = false;
    }

    size_t capacity = 1;
    *size = 0;
    TRY_INNER(LAGraph_Calloc((void **)symbols, capacity, sizeof(CFL_Symbol), msg));

    for (size_t i = 0; i < capacity; i++) {
        CFL_Symbol sym = {0};
        (*symbols)[i] = sym;
    }

    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_EWCNF rule = rules[i];

        int32_t prods[3] = {rule.nonterm, rule.prod_A, rule.prod_B};
        int bitmasks[3] = {LAGraph_EWNCF_INDEX_NONTERM, LAGraph_EWNCF_INDEX_PROD_A,
                           LAGraph_EWNCF_INDEX_PROD_B};

        for (size_t j = 0; j < 3; j++) {
            if (prods[j] == -1) {
                continue;
            }

            if (checked[prods[j]])
                continue;

            checked[prods[j]] = true;

            CFL_Symbol sym;
            sym.base_index = prods[j];
            sym.index = *size;
            sym.count = rule.indexed & bitmasks[j] ? rule.indexed_count : 0;

            for (size_t k = sym.base_index; k < sym.base_index + sym.count; k++) {
                checked[k] = true;
            }

            if (*size == capacity) {
                capacity *= 2;
                TRY_INNER(LAGraph_Realloc((void **)symbols, capacity, capacity / 2,
                                          sizeof(CFL_Symbol), msg));
            }

            (*symbols)[(*size)++] = sym;
        }
    }

    for (size_t i = 0; i < symbols_amount; i++) {
        if (checked[i])
            continue;

        if (*size == capacity) {
            capacity *= 2;
            TRY_INNER(LAGraph_Realloc((void **)symbols, capacity, capacity / 2,
                                      sizeof(CFL_Symbol), msg));
        }

        CFL_Symbol sym;
        sym.base_index = i;
        sym.index = *size;
        sym.count = 0;

        (*symbols)[(*size)++] = sym;
        checked[i] = true;
        // printf("Inserted (%ld, %ld, %ld)\n", sym.index, sym.base_index, sym.count);
    }

    FREE_INNER_WORK();

    return GrB_SUCCESS;
}

// Expands indexed grammar rules into a set of concrete rules.
//
//   Before: A_i -> B_i C     (indexed_count = 3)
//           D   -> E F
//   After:  A_0 -> B_0 C
//           A_1 -> B_1 C
//           A_2 -> B_2 C
//           D   -> E F
//
// Output: LAGraph_rule_EWCNF **new_rules and size_t *new_rules_count
static GrB_Info explode_rules(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                              LAGraph_rule_EWCNF **new_rules, size_t *new_rules_count,
                              char *msg) {
    *new_rules = NULL;

#undef FREE_INNER

#define FREE_INNER()                                                                     \
    {                                                                                    \
        LAGraph_Free((void **)new_rules, msg);                                           \
    }

    size_t new_rules_size = 0;
    size_t new_rules_capacity = 1;

    TRY_INNER(LAGraph_Calloc((void **)new_rules, new_rules_capacity,
                             sizeof(LAGraph_rule_EWCNF), msg));

    // explode rules
    for (size_t i_rule = 0; i_rule < rules_count; i_rule++) {
        LAGraph_rule_EWCNF rule = rules[i_rule];

        if (new_rules_size == new_rules_capacity) {
            TRY_INNER(LAGraph_Realloc((void **)new_rules, new_rules_capacity * 2,
                                      new_rules_capacity, sizeof(LAGraph_rule_EWCNF),
                                      msg));
            new_rules_capacity *= 2;
        }

        if (rule.indexed_count == 0) {
            (*new_rules)[new_rules_size++] = rule;
        } else {
            while (new_rules_size + rule.indexed_count >= new_rules_capacity) {
                TRY_INNER(LAGraph_Realloc((void **)new_rules, new_rules_capacity * 2,
                                          new_rules_capacity, sizeof(LAGraph_rule_EWCNF),
                                          msg));
                new_rules_capacity *= 2;
            }

            for (size_t rule_index = 0; rule_index < rule.indexed_count; rule_index++) {
                LAGraph_rule_EWCNF new_rule = rule;
                new_rule.indexed_count = 0;
                new_rule.indexed = 0;
                if (rule.nonterm != -1 && rule.indexed & LAGraph_EWNCF_INDEX_NONTERM) {
                    new_rule.nonterm = rule.nonterm + rule_index;
                }
                if (rule.prod_A != -1 && rule.indexed & LAGraph_EWNCF_INDEX_PROD_A) {
                    new_rule.prod_A = rule.prod_A + rule_index;
                }
                if (rule.prod_B != -1 && rule.indexed & LAGraph_EWNCF_INDEX_PROD_B) {
                    new_rule.prod_B = rule.prod_B + rule_index;
                }

                (*new_rules)[new_rules_size++] = new_rule;
            }
        }
    }

    *new_rules_count = new_rules_size;

    return GrB_SUCCESS;
}

// Splits a CFL_Matrix into an array of GrB_Matrix matrices

// If the matrix is not a horizontal or vertical block vector (i.e., its
// block_type is CELL), the underlying GrB_Matrix is extracted and copied
// into outputs[0]
//
// Otherwise, the matrix split into square sub-matrices of size (graph_size x graph_size)
// Parameters:
//   outputs       - [out] Caller-allocated array of GrB_Matrix to write results into
//                         Must have enough space for all sub-matrices
//   matrix        - [in]  Source CFL_Matrix to split
static GrB_Info split_CFL_matrix(GrB_Matrix *outputs, CFL_Matrix *matrix,
                                 int8_t optimizations) {
    char msg[LAGRAPH_MSG_LEN];
    CFL_Matrix *base_matrix;
    GrB_Index *nrows = NULL, *ncols = NULL;

#undef FREE_INNER

#define FREE_INNER()                                                                     \
    {                                                                                    \
        CFL_matrix_free(&base_matrix);                                                   \
        LAGraph_Free((void **)&nrows, msg);                                              \
        LAGraph_Free((void **)&ncols, msg);                                              \
    }

    if (matrix->block_type == CELL) {
        CFL_Matrix *result;
        TRY_INNER(CFL_matrix_to_base(&result, matrix, optimizations));
        TRY_INNER(GrB_Matrix_dup(outputs, result->base));
        TRY_INNER(CFL_matrix_free(&result));
        return GrB_SUCCESS;
    }

    TRY_INNER(CFL_matrix_to_base(&base_matrix, matrix, optimizations));

    // we can create nrows and ncols array with the same size, it will be mush easier than
    // calculate size of each array :)
    GrB_Index matrices_count = matrix->nrows > matrix->ncols
                                   ? matrix->nrows / matrix->ncols
                                   : matrix->ncols / matrix->nrows;

    TRY_INNER(LAGraph_Calloc((void **)&nrows, matrices_count, sizeof(GrB_Index), msg));
    TRY_INNER(LAGraph_Calloc((void **)&ncols, matrices_count, sizeof(GrB_Index), msg));

    GrB_Index graph_size = matrix->nrows < matrix->ncols ? matrix->nrows : matrix->ncols;
    for (size_t i = 0; i < matrices_count; i++) {
        nrows[i] = graph_size;
        ncols[i] = graph_size;
    }

    GrB_Index m = 0, n = 0;
    if (matrix->block_type == VEC_VERT) {
        m = matrices_count;
        n = 1;
    } else {
        m = 1;
        n = matrices_count;
    }

    TRY_INNER(GxB_Matrix_split(outputs, m, n, nrows, ncols, base_matrix->base, GrB_NULL));

    TRY_INNER(LAGraph_Free((void **)&nrows, msg));
    TRY_INNER(LAGraph_Free((void **)&ncols, msg));
    TRY_INNER(CFL_matrix_free(&base_matrix));

    return GrB_SUCCESS;
}

// Builds the mapping [old_index -> (new_index, base_index, indexed_count)]
// See get_new_symbols() for a detailed description of the mapping format
//
// Parameters:
//   map           - [out] Allocated array of CFL_Symbol. Caller must free
//   size          - [out] Number of entries in map
static GrB_Info get_new_symbols_map(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                                    size_t symbols_amount, CFL_Symbol **map, size_t *size,
                                    char *msg, int8_t optimizations) {
#undef FREE_INNER

#define FREE_INNER()                                                                     \
    {                                                                                    \
        LAGraph_Free((void **)map, msg);                                                 \
    }

    if (optimizations & OPT_BLOCK) {
        TRY_INNER(get_new_symbols(rules, rules_count, symbols_amount, map, size, msg));
    } else {
        *size = symbols_amount;
        TRY_INNER(LAGraph_Calloc((void **)map, symbols_amount, sizeof(CFL_Symbol), msg));
        for (size_t i = 0; i < symbols_amount; i++) {
            CFL_Symbol *sym = &(*map)[i];
            sym->index = i;
            sym->base_index = i;
            sym->count = 0;
        }
    }

    return GrB_SUCCESS;
}

// Builds a new array of adjacency matrices according to the symbol mapping
//
// Parameters:
//   new_adj_matrices_p - [out] Resulting matrix array. Caller must free
//                              (only if OPT_BLOCK is set).
static GrB_Info get_new_adj_matrices(const GrB_Matrix *adj_matrices, CFL_Symbol *map,
                                     GrB_Index map_size, GrB_Matrix **new_adj_matrices_p,
                                     char *msg, int8_t optimizations) {
#undef FREE_INNER

#define FREE_INNER()                                                                     \
    {                                                                                    \
        LAGraph_Free((void **)new_adj_matrices_p, msg);                                  \
    }

    GrB_Index n;
    TRY_INNER(GrB_Matrix_ncols(&n, adj_matrices[0]));

    if (!(optimizations & OPT_BLOCK)) {
        *new_adj_matrices_p = (GrB_Matrix *)adj_matrices;
        return GrB_SUCCESS;
    }

    TRY_INNER(
        LAGraph_Calloc((void **)new_adj_matrices_p, map_size, sizeof(GrB_Matrix), msg));

    for (size_t i = 0; i < map_size; i++) {
        CFL_Symbol sym = map[i];

        if (sym.count == 0) {
            TRY_INNER(
                GrB_Matrix_dup(&(*new_adj_matrices_p)[i], adj_matrices[sym.base_index]));
            continue;
        }

        GrB_Matrix new_col_matrix;
        TRY_INNER(GrB_Matrix_new(&new_col_matrix, GrB_BOOL, n * sym.count, n));
        GrB_Matrix *Tiles = (GrB_Matrix *)adj_matrices + sym.base_index;
        TRY_INNER(GxB_Matrix_concat(new_col_matrix, Tiles, sym.count, 1, GrB_NULL));
        (*new_adj_matrices_p)[i] = new_col_matrix;
    }

    return GrB_SUCCESS;
}

// Remaps rule symbol indices according to the symbol mapping
//
// Parameters:
//   new_rules       - [out] Allocated output rule array. Caller must free
//   new_rules_count - [out] Number of rules written to new_rules
static GrB_Info get_new_rules(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                              CFL_Symbol *map, size_t map_size,
                              LAGraph_rule_EWCNF **new_rules, size_t *new_rules_count,
                              char *msg, int8_t optimizations) {
    *new_rules_count = 0;
#undef FREE_INNER

#define FREE_INNER()                                                                     \
    {                                                                                    \
        LAGraph_Free((void **)new_rules, msg);                                           \
    }

    if (!(optimizations & OPT_BLOCK)) {
        TRY_INNER(explode_rules(rules, rules_count, new_rules, new_rules_count, msg));
        return GrB_SUCCESS;
    }

    TRY_INNER(
        LAGraph_Calloc((void **)new_rules, rules_count, sizeof(LAGraph_rule_EWCNF), msg));
    for (size_t i = 0; i < rules_count; i++) {
        LAGraph_rule_EWCNF rule = rules[i];
        LAGraph_rule_EWCNF new_rule = rule;

        for (size_t i_sym = 0; i_sym < map_size; i_sym++) {
            CFL_Symbol sym = map[i_sym];

            if (rule.nonterm != -1 && rule.nonterm == sym.base_index) {
                new_rule.nonterm = sym.index;
            }
            if (rule.prod_A != -1 && rule.prod_A == sym.base_index) {
                new_rule.prod_A = sym.index;
            }
            if (rule.prod_B != -1 && rule.prod_B == sym.base_index) {
                new_rule.prod_B = sym.index;
            }
        }

        (*new_rules)[(*new_rules_count)++] = new_rule;
    }

    return GrB_SUCCESS;
}

// LAGraph_CFL_reachability_multsrc_adv: Optimized Multiple-Source Context-Free
// Language Reachability Matrix-Based Algorithm
//
// This function determines the set of vertex pairs (u, v) in a graph (represented by
// adjacency matrices) such that u is in the given set of source vertices `src`, and
// there is a path from u to v, where the edge labels form a word from the language
// generated by the context-free grammar (represented by `rules`).

GrB_Info LAGraph_CFL_reachability_multsrc_adv
(
    // Output
    GrB_Matrix *output, // A handle for the matrix containing results.
                        //
                        // output: (i, j) = true if and only if there is a path
                        // from node i to node j whose edge labels form a word
                        // derivable from the specified CFG.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, symbols_amount - 1]).
    GrB_Index *src,                  // Array of source vertices
    size_t src_count,                // The total number of source vertices
    size_t symbols_amount,           // The total number of terminals and nonterminals
    const LAGraph_rule_EWCNF *rules, // The rules of the CFG.
    size_t rules_count,              // The total number of rules in the CFG.
    char *msg,                       // Message string for error reporting.
    int8_t opt_mask                  // Optimizations mask
)
{
    CFL_Matrix **T;
    CFL_Matrix **dT;
    CFL_Matrix **TSrc;
    CFL_Matrix **Adj;
    CFL_Matrix *MSrc;
    CFL_Matrix *M1;
    CFL_Matrix *M2;
    CFL_Matrix *Temp1;
    CFL_Matrix *Temp2;
    CFL_Matrix *iden = NULL;
    CFL_Matrix *A;
    GrB_Vector a;
    GrB_Matrix identity_matrix = NULL;
    GrB_Index *nnzs_T = NULL;
    GrB_Index *nnzs_TSrc_B = NULL;
    GrB_Index *nnzs_TSrc_C = NULL;
    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting
    bool iso_flag = false;
    GrB_Scalar true_scalar;
    GrB_Vector ones_vec;

    // for OPT_BLOCK optimization
    size_t new_symbols_amount = 0;
    GrB_Matrix *new_adj_matrices = NULL;
    LAGraph_rule_EWCNF *new_rules = NULL;
    size_t new_rules_count = 0;
    CFL_Symbol *to_new_symbols_map = NULL;

    // TODO: LG_ASSERT_MSG
    if (symbols_amount <= 0) {
        return GrB_INVALID_VALUE;
    }

    if (rules_count <= 0) {
        return GrB_INVALID_VALUE;
    }

    if (!output || !rules || !adj_matrices || !src) {
        return GrB_NULL_POINTER;
    }

    if (src_count <= 0) {
        return GrB_INVALID_VALUE;
    }

    // Find null adjacency matrices
    bool found_null = false;
    for (int32_t i = 0; i < symbols_amount; i++) {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null) {
            ADD_TO_MSG("Adjacency matrices with these indices are null: ");
            ADD_TO_MSG("%d", i);
        } else {
            ADD_TO_MSG(", %d", i);
        }

        found_null = true;
    }
    if (found_null) {
        return GrB_NULL_POINTER;
    }

    GrB_Index n;
    TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));
    if (n < src_count) {
        return GrB_INVALID_VALUE;
    }

    TRY(get_new_symbols_map(rules, rules_count, symbols_amount, &to_new_symbols_map,
                            &new_symbols_amount, msg, opt_mask));
    TRY(get_new_adj_matrices(adj_matrices, to_new_symbols_map, new_symbols_amount,
                             &new_adj_matrices, msg, opt_mask));
    TRY(get_new_rules(rules, rules_count, to_new_symbols_map, new_symbols_amount,
                      &new_rules, &new_rules_count, msg, opt_mask));


    // Arrays for processing rules
    size_t eps_rules[new_rules_count], eps_rules_count = 0;   // [Variable -> eps]
    size_t term_rules[new_rules_count], term_rules_count = 0; // [Variable -> term]
    size_t bin_rules[new_rules_count], bin_rules_count = 0;   // [Variable -> AB]

    // Process rules
    typedef struct {
        size_t count;
        size_t len_indices_str;
        char indices_str[LAGRAPH_MSG_LEN];
    } rule_error_s;

    rule_error_s term_err = {0};
    rule_error_s nonterm_err = {0};
    rule_error_s invalid_err = {0};

    for (size_t i = 0; i < rules_count; i++) {
        // LAGraph_rule_WCNF rule = rules[i];
        LAGraph_rule_EWCNF rule = new_rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // Check that all rules are well-formed
        if (rule.nonterm < 0 || (size_t)rule.nonterm >= new_symbols_amount) {
            ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
        }

        // [Variable -> eps]
        if (is_rule_eps) {
            eps_rules[eps_rules_count++] = i;

            continue;
        }

        // [Variable -> term]
        if (is_rule_term) {
            term_rules[term_rules_count++] = i;

            if (rule.prod_A < -1 || (size_t)rule.prod_A >= new_symbols_amount) {
                ADD_INDEX_TO_ERROR_RULE(term_err, i);
            }

            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin) {
            bin_rules[bin_rules_count++] = i;

            if (rule.prod_A < -1 || (size_t)rule.prod_A >= new_symbols_amount ||
                rule.prod_B < -1 || (size_t)rule.prod_B >= new_symbols_amount) {
                ADD_INDEX_TO_ERROR_RULE(nonterm_err, i);
            }

            continue;
        }

        // [Variable -> _ B]
        ADD_INDEX_TO_ERROR_RULE(invalid_err, i);
    }

    if (term_err.count + nonterm_err.count + invalid_err.count > 0) {
        ADD_TO_MSG("Count of invalid rules: %ld.\n",
                   term_err.count + nonterm_err.count + invalid_err.count);

        if (nonterm_err.count > 0) {
            ADD_TO_MSG("Non-terminals must be in range [0, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indices_str)
        }
        if (term_err.count > 0) {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indices_str)
        }
        if (invalid_err.count > 0) {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indices_str)
        }

        return GrB_INVALID_VALUE;
    }

    GrB_Scalar_new(&true_scalar, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(true_scalar, true);

    TRY(LAGraph_Calloc((void **) &T, symbols_amount, sizeof(CFL_Matrix), msg));
    TRY(LAGraph_Calloc((void **) &dT, symbols_amount, sizeof(CFL_Matrix), msg));
    TRY(LAGraph_Calloc((void **) &TSrc, symbols_amount, sizeof(CFL_Matrix), msg));

    TRY(LAGraph_Calloc((void **) &nnzs_T, symbols_amount, sizeof(GrB_Index), msg));
    TRY(LAGraph_Calloc((void **) &nnzs_TSrc_B, symbols_amount, sizeof(GrB_Index), msg));
    TRY(LAGraph_Calloc((void **) &nnzs_TSrc_C, symbols_amount, sizeof(GrB_Index), msg));

    TRY(LAGraph_Calloc((void **) &Adj, symbols_amount, sizeof(CFL_Matrix), msg));

    TRY(GrB_Vector_new(&ones_vec, GrB_BOOL, n));
    TRY(GrB_Vector_assign_BOOL(ones_vec, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));

    TRY(GrB_Matrix_diag(&identity_matrix, ones_vec, 0));
    TRY(CFL_matrix_from_base(&iden, identity_matrix));

    for (int32_t i = 0; i < new_symbols_amount; i++) {
        TRY(CFL_matrix_from_base(&Adj[i], new_adj_matrices[i]));
    }

    // Create nonterms matrices
    for (int32_t i = 0; i < new_symbols_amount; i++) {
        if (opt_mask & OPT_LAZY) {
            TRY(CFL_matrix_create_lazy(&T[i], n, n));
        } else {
            TRY(CFL_matrix_create(&T[i], n, n));
        }

        TRY(CFL_matrix_create(&dT[i], n, n));

        TRY(CFL_matrix_create(&TSrc[i], n, n));
    }

    for (int32_t i = 0; i < src_count; i++) {
        GrB_Matrix_setElement(TSrc[0]->base, true, src[i], src[i]);
    }
    TRY(CFL_matrix_update(TSrc[0]));
    // TRY(CFL_matrix_from_base(&TSrc[0], TSrc[0]->base));

    TRY(CFL_matrix_create(&MSrc, n, n));
    TRY(CFL_matrix_create(&M1, n, n));
    TRY(CFL_matrix_create(&M2, n, n));
    TRY(CFL_matrix_create(&Temp1, n, n));
    TRY(CFL_matrix_create(&Temp2, n, n));
    TRY(CFL_matrix_create(&A, n, n));
    TRY(CFL_dup(MSrc, TSrc[0], opt_mask));
    TRY(GrB_Vector_new(&a, GrB_BOOL, n));

    // Rule [Variable -> term]
    for (size_t i = 0; i < term_rules_count; i++) {
        // LAGraph_rule_WCNF term_rule = rules[term_rules[i]];
        LAGraph_rule_EWCNF term_rule = new_rules[term_rules[i]];

        // I'd like to get rid of double operation
        GRB_TRY(CFL_wise(T[term_rule.nonterm], T[term_rule.nonterm], Adj[term_rule.prod_A], true, opt_mask));
        GRB_TRY(CFL_wise(dT[term_rule.nonterm], dT[term_rule.nonterm], Adj[term_rule.prod_A], true, opt_mask));
    }

    // Rule [Variable -> eps]
    for (size_t i = 0; i < eps_rules_count; i++) {
        // LAGraph_rule_WCNF eps_rule = rules[eps_rules[i]];
        LAGraph_rule_EWCNF eps_rule = new_rules[eps_rules[i]];

        GRB_TRY(CFL_wise(T[eps_rule.nonterm], T[eps_rule.nonterm], iden, true, opt_mask));
        GRB_TRY(CFL_wise(dT[eps_rule.nonterm], dT[eps_rule.nonterm], iden, true, opt_mask));
    }

    // Rule [Variable -> Variable1 Variable2]
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < bin_rules_count; i++) {
            // LAGraph_rule_WCNF bin_rule = rules[bin_rules[i]];
            LAGraph_rule_EWCNF bin_rule = new_rules[bin_rules[i]];

            TRY(CFL_mxm(M1, TSrc[bin_rule.nonterm], T[bin_rule.prod_A], false, false, opt_mask));

            TRY(CFL_wise(T[bin_rule.nonterm], T[bin_rule.nonterm], dT[bin_rule.nonterm], false, opt_mask));

            // ? <- this means that without this step, unit tests pass,
            // although algorithm requires it.

            TRY(CFL_mxm(M2, TSrc[bin_rule.nonterm], dT[bin_rule.prod_A], false, false, opt_mask));

            TRY(CFL_mxm(Temp1, M1, dT[bin_rule.prod_B], false, false, opt_mask));

            // ?
            TRY(CFL_mxm(Temp2, M2, T[bin_rule.prod_B], false, false, opt_mask));

            TRY(CFL_dup(dT[bin_rule.nonterm], Temp1, opt_mask));

            TRY(CFL_wise(dT[bin_rule.nonterm], dT[bin_rule.nonterm], Temp2, false, opt_mask));

            // ?
            TRY(CFL_rsub(dT[bin_rule.nonterm], T[bin_rule.nonterm], opt_mask));

            // Update source vertices matrix to find appropriate paths only
            // M1[i, j] == 1 => A[j, j] == 1
            TRY(GrB_vxm(a, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, ones_vec, M1->base, GrB_NULL));
            TRY(GrB_Matrix_free(&A->base));
            TRY(GrB_Matrix_diag(&A->base, a, 0));
            TRY(CFL_matrix_from_base(&A, A->base));

            TRY(CFL_wise(TSrc[bin_rule.prod_A], TSrc[bin_rule.prod_A], TSrc[bin_rule.nonterm], false, opt_mask));

            TRY(CFL_wise(TSrc[bin_rule.prod_B], TSrc[bin_rule.prod_B], A, false, opt_mask));

            // Check if any of the matrices changed. If not, job is done.
            GrB_Index nnz_T, nnz_TSrc_B, nnz_TSrc_C;

            CFL_matrix_update(T[bin_rule.nonterm]);
            CFL_matrix_update(TSrc[bin_rule.prod_A]);
            CFL_matrix_update(TSrc[bin_rule.prod_B]);

            nnz_T = T[bin_rule.nonterm]->nvals;
            nnz_TSrc_B = TSrc[bin_rule.prod_A]->nvals;
            nnz_TSrc_C = TSrc[bin_rule.prod_B]->nvals;

            changed = changed || (nnzs_T[bin_rule.nonterm] != nnz_T);
            changed = changed || (nnzs_TSrc_B[bin_rule.prod_A] != nnz_TSrc_B);
            changed = changed || (nnzs_TSrc_C[bin_rule.prod_B] != nnz_TSrc_C);

            nnzs_T[bin_rule.nonterm] = nnz_T;
            nnzs_TSrc_B[bin_rule.prod_A] = nnz_TSrc_B;
            nnzs_TSrc_C[bin_rule.prod_B] = nnz_TSrc_C;
        }
    }

    GRB_TRY(CFL_mxm(MSrc, MSrc, T[0], false, false, opt_mask));

    // get outputs matrix
    if (MSrc->base_matrices_count == 0) {
        TRY(GrB_Matrix_dup(output, MSrc->base));
    } else {
        CFL_Matrix *res;
        CFL_matrix_to_base(&res, MSrc, opt_mask);
        TRY(GrB_Matrix_dup(output, res->base));
        CFL_matrix_free(&res);
    }

    // TRY(split_CFL_matrix(output, MSrc, opt_mask));

    LG_FREE_ALL;
    return GrB_SUCCESS;
}
