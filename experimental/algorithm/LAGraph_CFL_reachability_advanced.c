//------------------------------------------------------------------------------
// LAGraph_CFL_reachability.c: Context-Free Language Reachability Matrix-Based
// Algorithm
// ------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Ilhom Kombaev, Semyon Grigoriev, St. Petersburg State University.

//------------------------------------------------------------------------------

// Code is based on the "A matrix-based CFPQ algorithm" described in the
// following paper: * Rustam Azimov, Semyon Grigorev, "Context-Free Path
// Querying Using Linear Algebra", URL:
// https://disser.spbu.ru/files/2022/disser_azimov.pdf

#define LG_FREE_WORK                                                                     \
    {                                                                                    \
        TRY_INNER(CFL_matrix_free(&iden));                                               \
        TRY_INNER(LAGraph_Free((void **)&to_new_symbols_map, msg));                      \
        TRY_INNER(LAGraph_Free((void **)&new_rules, msg));                               \
        for (size_t i = 0; i < new_symbols_amount; i++) {                                \
            TRY_INNER(CFL_matrix_free(&temp_matrices[i]));                               \
            TRY_INNER(CFL_matrix_free(&delta_matrices[i]));                              \
            TRY_INNER(CFL_matrix_free(&matrices[i]));                                    \
            if (new_adj_matrices != adj_matrices) {                                      \
                TRY_INNER(GrB_free(&new_adj_matrices[i]));                               \
            }                                                                            \
        }                                                                                \
        if (new_adj_matrices != adj_matrices) {                                          \
            TRY_INNER(LAGraph_Free((void **)&new_adj_matrices, msg));                    \
        }                                                                                \
        TRY_INNER(LAGraph_Free((void **)&delta_matrices, msg));                          \
        TRY_INNER(LAGraph_Free((void **)&matrices, msg));                                \
        TRY_INNER(LAGraph_Free((void **)&temp_matrices, msg));                           \
    }

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        LG_FREE_WORK;                                                                    \
    }

#include "LAGraph_CFL_optimized_matrix_opt.h"
#include "LG_internal.h"
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

#define ADD_INDEX_TO_ERROR_RULE(rule, i)                                                 \
    {                                                                                    \
        rule.len_indexes_str += snprintf(rule.indexes_str + rule.len_indexes_str,        \
                                         LAGRAPH_MSG_LEN - rule.len_indexes_str,         \
                                         rule.count == 0 ? "%ld" : ", %ld", i);          \
        rule.count++;                                                                    \
    }

#define BENCH_CFL_REACHBILITY false

#if BENCH_CFL_REACHBILITY
    #define IS_ISO(matrix, str)                                                          \
        {                                                                                \
            bool iso_flag;                                                               \
            GrB_Index nnz;                                                               \
            GxB_Matrix_iso(&iso_flag, matrix);                                           \
            GrB_Matrix_nvals(&nnz, matrix);                                              \
            if (!iso_flag && nnz) {                                                      \
                printf("-----ISO ALERT----- (%s)\n", str);                               \
                GxB_print(matrix, 1);                                                    \
                printf("-------------------\n");                                         \
            }                                                                            \
        }

    #define TIMER_START()                                                                \
        {                                                                                \
            start_time = LAGraph_WallClockTime();                                        \
        }

    #define TIMER_STOP(label, accumulator)                                               \
        {                                                                                \
            end_time = LAGraph_WallClockTime();                                          \
            printf("%s %.3fs\n", label, end_time - start_time);                          \
            if (accumulator != NULL) {                                                   \
                *(accumulator) += (end_time - start_time);                               \
            }                                                                            \
        }

    #define IS_ROW(matrix, str)                                                          \
        {                                                                                \
            int32_t orientation;                                                         \
            GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT);                 \
            if (orientation != GrB_ROWMAJOR) {                                           \
                printf("-----NOT A ROW----- (%s)\n", str);                               \
                GxB_print(matrix, 1);                                                    \
                printf("-------------------\n");                                         \
            }                                                                            \
        }

    #define IS_COL(matrix, str)                                                          \
        {                                                                                \
            int32_t orientation;                                                         \
            GrB_get(matrix, &orientation, GrB_STORAGE_ORIENTATION_HINT);                 \
            if (orientation != GrB_COLMAJOR) {                                           \
                printf("-----NOT A COL----- (%s)\n", str);                               \
                GxB_print(matrix, 1);                                                    \
                printf("-------------------\n");                                         \
            }                                                                            \
        }
#else
    #define IS_ISO(matrix, str)
    #define TIMER_START()
    #define TIMER_STOP(label, accumulator)
    #define IS_ROW(matrix, str)
    #define IS_COL(matrix, str)
#endif
// clang-format on

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

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

#define TRY_I(GrB_method)                                                                \
    {                                                                                    \
        GrB_Info LG_GrB_Info = GrB_method;                                               \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                 \
            fprintf(stderr,                                                              \
                    "LAGraph failure (file %s, line %d) (Iteration: %d, i: %d): \n",     \
                    __FILE__, __LINE__, iteration, i);                                   \
            return (LG_GrB_Info);                                                        \
        }                                                                                \
    }

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

// LAGraph_CFL_reachability_adv: Context-Free Language Reachability Matrix-Based
// Algorithm
//
// This function determines the set of vertex pairs (u, v) in a graph (represented by
// adjacency matrices) such that there is a path from u to v, where the edge labels
// form a word from the language generated by the context-free grammar (represented by
// `rules`).
//
// Terminals and non-terminals are enumerated by integers starting from zero.
// The start non-terminal is the non-terminal with index 0.
//
// Example:
//
// Graph:
// ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
// │ 0 ├───► 1 ├───► 2 ├───► 3 ├───► 4 │
// └───┘ a └─┬─┘ a └─▲─┘ b └───┘ b └───┘
//           │       │
//           │ ┌───┐ │
//          a└─► 5 ├─┘b
//             └───┘
//
// Grammar: S -> aSb | ab
//
// There are paths from node [1] to node [3] and from node [1] to node [2] that form
// the word "ab" ([1]-a->[2]-b->[3] and [1]-a->[5]-b->[2]). The word "ab" is in the
// language generated by our context-free grammar, so the pairs (1, 3) and (1, 2) will
// be included in the result.
//
// Note: It doesn't matter how many paths exist from node [A] to node [B] that form a
// word in the language. If at least one path exists, the pair ([A], [B]) will be
// included in the result.
//
// In contrast, the path from node [1] to node [4] forms the word "abb"
// ([1]-a->[2]-b->[3]-b->[4]) and the word "abbb" ([1]-a->[5]-b->[2]-b->[3]-b->[4]).
// The words "aab" and "abbb" are not in the language, so the pair (1, 4) will not be
// included in the result.
//
// With this graph and grammar, we obtain the following results:
// (0, 4) - because there exists a path (0-1-2-3-4) that forms the word "aabb"
// (1, 3) - because there exists a path (1-2-3) that forms "ab"
// (1, 2) - because there exists a path (1-5-2) that forms the word "ab"
// (0, 3) - because there exists a path (0-1-5-2-3) that forms the word "aabb"
GrB_Info LAGraph_CFL_reachability_adv(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to the count
                         // of symbols (symbols_amount).
                         //
                         // Each matrix is square, with size equal to the number of
                         // vertices in the graph. Matrices are allocated by the
                         // caller, not by this method.
                         //
                         // outputs[k]: (i, j) = true if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the symbol 'k' of the specified CFG.
                         //
                         // Note: output[t], where t is the index of a terminal, will
                         // be an exact copy of adj_matrices[t].
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // symbols (symbols_amount).
                                    //
                                    // Each matrix is square, with size equal to the
                                    // number of vertices in the graph.
                                    //
                                    // adj_matrices[i]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the symbol corresponding to index 'i' (where i is
                                    // in the range [0, symbols_amount - 1]).
                                    //
                                    // Note: adj_matrices[N], where N is the index of a
                                    // nonterminal, must be initialized as an empty square
                                    // matrix of size symbols_amount.
    size_t symbols_amount,          // Count of terminal and nonterminals
    const LAGraph_rule_EWCNF *rules, // The rules of the CFG.
                                     // Warning: now not ready for N -> A rules, where
                                     // is N and A are nonterminals.
    size_t rules_count,              // The total number of rules in the CFG.
    char *msg,                       // Message string for error reporting.
    int8_t optimizations             // Optimizations flags
) {
#undef FREE_INNER

#define FREE_INNER()

    // Declare workspace and clear the msg string, if not NULL
    CFL_Matrix **delta_matrices, **matrices, **temp_matrices;
    CFL_Matrix *iden = NULL;
    GrB_Matrix identity_matrix = NULL;

    // for OPT_BLOCK optimization
    size_t new_symbols_amount = 0;
    GrB_Matrix *new_adj_matrices = NULL;
    LAGraph_rule_EWCNF *new_rules = NULL;
    size_t new_rules_count = 0;
    CFL_Symbol *to_new_symbols_map = NULL;

    LG_CLEAR_MSG;
    size_t msg_len = 0; // For error formatting

    TRY(LAGraph_Calloc((void **)&delta_matrices, symbols_amount, sizeof(CFL_Matrix *),
                       msg));
    TRY(LAGraph_Calloc((void **)&matrices, symbols_amount, sizeof(CFL_Matrix *), msg));
    TRY(LAGraph_Calloc((void **)&temp_matrices, symbols_amount, sizeof(CFL_Matrix *),
                       msg));

    LG_ASSERT_MSG(symbols_amount > 0, GrB_INVALID_VALUE,
                  "The number of symbols must be greater than zero.");
    LG_ASSERT_MSG(rules_count > 0, GrB_INVALID_VALUE,
                  "The number of rules must be greater than zero.");
    LG_ASSERT_MSG(outputs != NULL, GrB_NULL_POINTER, "The outputs array cannot be null.");
    LG_ASSERT_MSG(rules != NULL, GrB_NULL_POINTER, "The rules array cannot be null.");
    LG_ASSERT_MSG(adj_matrices != NULL, GrB_NULL_POINTER,
                  "The adjacency matrices array cannot be null.");

    // Find null adjacency matrices
    bool found_null = false;
    for (size_t i = 0; i < symbols_amount; i++) {
        if (adj_matrices[i] != NULL)
            continue;

        if (!found_null) {
            ADD_TO_MSG("Adjacency matrices with these indexes are null: ");
            ADD_TO_MSG("%ld", i);
        } else {
            ADD_TO_MSG(", %ld", i);
        }

        found_null = true;
    }

    if (found_null) {
        LG_FREE_ALL;
        return GrB_NULL_POINTER;
    }

    GrB_Index n;
    TRY(GrB_Matrix_ncols(&n, adj_matrices[0]));

    TRY(get_new_symbols_map(rules, rules_count, symbols_amount, &to_new_symbols_map,
                            &new_symbols_amount, msg, optimizations));
    TRY(get_new_adj_matrices(adj_matrices, to_new_symbols_map, new_symbols_amount,
                             &new_adj_matrices, msg, optimizations));
    TRY(get_new_rules(rules, rules_count, to_new_symbols_map, new_symbols_amount,
                      &new_rules, &new_rules_count, msg, optimizations));

    // Arrays for processing rules
    size_t eps_rules[new_rules_count], eps_rules_count = 0;   // [Variable -> eps]
    size_t term_rules[new_rules_count], term_rules_count = 0; // [Variable -> term]
    size_t bin_rules[new_rules_count], bin_rules_count = 0;   // [Variable -> AB]

    // Process rules
    typedef struct {
        size_t count;
        size_t len_indexes_str;
        char indexes_str[LAGRAPH_MSG_LEN];
    } rule_error_s;
    rule_error_s term_err = {0};
    rule_error_s nonterm_err = {0};
    rule_error_s invalid_err = {0};
    for (size_t i = 0; i < new_rules_count; i++) {
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
            ADD_TO_MSG("Indexes of invalid rules: %s\n", nonterm_err.indexes_str)
        }
        if (term_err.count > 0) {
            ADD_TO_MSG("Terminals must be in range [-1, nonterms_count). ");
            ADD_TO_MSG("Indexes of invalid rules: %s\n", term_err.indexes_str)
        }
        if (invalid_err.count > 0) {
            ADD_TO_MSG("[Variable -> _ B] type of rule is not acceptable. ");
            ADD_TO_MSG("Indexes of invalid rules: %.120s\n", invalid_err.indexes_str)
        }

        LG_FREE_ALL;
        return GrB_INVALID_VALUE;
    }

    // Create symbol matrices
    for (size_t i = 0; i < new_symbols_amount; i++) {
        GrB_Index nrows;
        TRY(GrB_Matrix_nrows(&nrows, new_adj_matrices[i]));
        GrB_Index ncols;
        TRY(GrB_Matrix_ncols(&ncols, new_adj_matrices[i]));

        GrB_Matrix new_adj_matrix;
        TRY(GrB_Matrix_dup(&new_adj_matrix, new_adj_matrices[i]));
        TRY(CFL_matrix_from_base(&delta_matrices[i], new_adj_matrix));

        if (optimizations & OPT_LAZY) {
            TRY(CFL_matrix_create_lazy(&matrices[i], nrows, ncols));
        } else {
            TRY(CFL_matrix_create(&matrices[i], nrows, ncols));
        }

        TRY(CFL_matrix_create(&temp_matrices[i], nrows, ncols));
    }

    // Rule [Variable -> term]
    for (size_t i = 0; i < term_rules_count; i++) {
        LAGraph_rule_EWCNF term_rule = new_rules[term_rules[i]];
        CFL_Matrix *nonterm_matrix = delta_matrices[term_rule.nonterm];
        CFL_Matrix *term_matrix = delta_matrices[term_rule.prod_A];

        TRY(CFL_wise(nonterm_matrix, nonterm_matrix, term_matrix, true, optimizations));
    }

    GrB_Vector v_diag;
    TRY(GrB_Vector_new(&v_diag, GrB_BOOL, n));
    TRY(GrB_Vector_assign_BOOL(v_diag, GrB_NULL, GrB_NULL, true, GrB_ALL, n, NULL));
    TRY(GrB_Matrix_diag(&identity_matrix, v_diag, 0));
    TRY(GrB_Vector_free(&v_diag));
    TRY(CFL_matrix_from_base(&iden, identity_matrix));

    // Rule [Variable -> eps]
    for (size_t i = 0; i < eps_rules_count; i++) {
        LAGraph_rule_EWCNF eps_rule = new_rules[eps_rules[i]];
        CFL_Matrix *nonterm_matrix = delta_matrices[eps_rule.nonterm];

        TRY(CFL_wise(nonterm_matrix, nonterm_matrix, iden, true, optimizations));
    }

    // Rule [Variable -> Variable1 Variable2]
    double start_time, end_time;
    bool changed = true;
    size_t iteration = 0;
    double mxm1 = 0.0;
    double wise1 = 0.0;
    double mxm2 = 0.0;
    double wise2 = 0.0;
    double rsubt = 0.0;
    while (changed) {
        iteration++;
        changed = false;

#if BENCH_CFL_REACHBILITY
        printf("\n--- ITERATARION %ld ---\n", iteration);
#endif

        for (size_t i = 0; i < new_symbols_amount; i++) {
            TRY_I(CFL_matrix_free(&temp_matrices[i]));
            TRY_I(CFL_matrix_create(&temp_matrices[i], matrices[i]->nrows,
                                    matrices[i]->ncols));
        }

        TIMER_START();
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_EWCNF bin_rule = new_rules[bin_rules[i]];
            CFL_Matrix *A = matrices[bin_rule.prod_A];
            CFL_Matrix *B = delta_matrices[bin_rule.prod_B];
            CFL_Matrix *C = temp_matrices[bin_rule.nonterm];

            TRY_I(CFL_mxm(C, A, B, true, false, optimizations));
        }
        TIMER_STOP("MXM 1", &mxm1);

        TIMER_START()
        for (size_t i = 0; i < new_symbols_amount; i++) {
            CFL_Matrix *A = delta_matrices[i];
            CFL_Matrix *C = matrices[i];

            TRY_I(CFL_wise(C, C, A, false, optimizations));
        }
        TIMER_STOP("WISE 1", &wise1);

        TIMER_START()
        for (size_t i = 0; i < bin_rules_count; i++) {
            LAGraph_rule_EWCNF bin_rule = new_rules[bin_rules[i]];
            CFL_Matrix *A = matrices[bin_rule.prod_B];
            CFL_Matrix *B = delta_matrices[bin_rule.prod_A];
            CFL_Matrix *C = temp_matrices[bin_rule.nonterm];

            TRY(CFL_mxm(C, A, B, true, true, optimizations));
        }
        TIMER_STOP("MXM 2", &mxm2);

        // Rule [Variable -> term]
        for (size_t i = 0; i < term_rules_count; i++) {
            LAGraph_rule_EWCNF term_rule = new_rules[term_rules[i]];
            CFL_Matrix *A = temp_matrices[term_rule.nonterm];
            CFL_Matrix *B = delta_matrices[term_rule.prod_A];

            TRY_I(CFL_wise(A, A, B, true, optimizations));
        }

        TIMER_START();
        for (size_t i = 0; i < new_symbols_amount; i++) {
            TRY_I(CFL_dup(delta_matrices[i], temp_matrices[i], optimizations));
        }
        TIMER_STOP("WISE 2 (copy)", &wise2);

        TIMER_START();
        for (size_t i = 0; i < new_symbols_amount; i++) {
            CFL_Matrix *A = matrices[i];
            CFL_Matrix *C = delta_matrices[i];

            TRY_I(CFL_rsub(C, A, optimizations));
        }
        TIMER_STOP("WISE 3 (MASK)", &rsubt);

        size_t new_nnz = 0;
        for (size_t i = 0; i < new_symbols_amount; i++) {
            TRY(CFL_matrix_update(delta_matrices[i]));
            new_nnz += delta_matrices[i]->nvals;
        }

        if (new_nnz != 0) {
            changed = true;
        }
    }

#if BENCH_CFL_REACHBILITY
    printf("MXM1: %.3f, wise1: %.3f, MXM2: %.3f, wise2: %.3f, rsub: %.3f", mxm1, wise1,
           mxm2, wise2, rsubt);
#endif

    // get outputs matrices
    for (size_t i = 0; i < new_symbols_amount; i++) {
        CFL_Symbol sym = to_new_symbols_map[i];
        TRY(split_CFL_matrix(outputs + sym.base_index, matrices[i], optimizations));
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}

// Helper function to free the output matrix of LAGraph_CFL_AllPaths, which contains elements of type AllPathsElem with dynamically allocated arrays of intermediate vertices. 
static void free_AllPaths_matrix(GrB_Matrix* ptr_output) 
{
  GxB_Iterator iterator;
  GxB_Iterator_new(&iterator);
  GrB_Info info = GxB_Matrix_Iterator_attach(iterator, *ptr_output, NULL);
  info = GxB_Matrix_Iterator_seek(iterator, 0);
  AllPathsElem val;

  while (info != GxB_EXHAUSTED)
  {
    GxB_Iterator_get_UDT(iterator, (void*)&val);
    if (val.n > 1 && val.data.middle != NULL) {
      free(val.data.middle);
    }
    info = GxB_Matrix_Iterator_next(iterator);
  }

  GrB_free(&iterator);
  GrB_free(ptr_output);
}

// Free outputs and all_paths_ptr_t after you have finished working with the output matrices from LAGraph_CFL_AllPaths.
// do outputs = NULL, all_paths_ptr_t = NULL after LAGraph_CFL_AllPaths_free_outputs
GrB_Info LAGraph_CFL_AllPaths_free_outputs(GrB_Matrix* outputs, int64_t nonterms_count, GrB_Type* all_paths_ptr_t)
{
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
