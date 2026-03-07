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
        TRY(CFL_matrix_free(&iden));                                                     \
        TRY(LAGraph_Free((void **)&symbols, msg));                                       \
        TRY(LAGraph_Free((void **)&new_rules, msg));                                     \
        for (size_t i = 0; i < new_symbols_amount; i++) {                                \
            TRY(CFL_matrix_free(&temp_matrices[i]));                                     \
            TRY(CFL_matrix_free(&delta_matrices[i]));                                    \
            TRY(CFL_matrix_free(&matrices[i]));                                          \
            if (new_adj_matrices != adj_matrices) {                                      \
                TRY(GrB_free(&new_adj_matrices[i]));                                     \
            }                                                                            \
        }                                                                                \
        if (new_adj_matrices != adj_matrices) {                                          \
            TRY(LAGraph_Free((void **)&new_adj_matrices, msg));                          \
        }                                                                                \
        TRY(LAGraph_Free((void **)&delta_matrices, msg));                                \
        TRY(LAGraph_Free((void **)&matrices, msg));                                      \
        TRY(LAGraph_Free((void **)&temp_matrices, msg));                                 \
    }

#define LG_FREE_ALL                                                                      \
    {                                                                                    \
        for (size_t i = 0; i < symbols_amount; i++) {                                    \
            /* TODO: delete it and free actual matrices */                               \
            /* GrB_free(&T[i]); */                                                       \
        }                                                                                \
                                                                                         \
        LG_FREE_WORK;                                                                    \
    }

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

// If we use OPT_BLOCK optimization we must group indexed symbols together
// After that we get mapping [old_num -> (index, base_index, count_indexed)]
// Where index is new index, base_index is old index and count_indexed count of indexed
// symbols
//
// Example
// (0) S -> (0, 0, 0)
// (1) A_0 -> (1, 1, 3)
// (2) A_1 ->
// (3) A_2 ->
// (4) B_0 -> (2, 4, 2)
// (5) B_1 ->
// (6) C -> (3, 6, 0)
// (7) a -> (4, 7, 0)
//
// This helps us to create new compact array of matrices or explode rules without
// OPT_BLOCK optimization
static GrB_Info get_new_symbols(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                                size_t symbols_amount, CFL_Symbol **symbols, size_t *size,
                                char *msg) {
    bool *checked = NULL;
    int LG_status = 0;

    LG_status = LAGraph_Calloc((void **)&checked, symbols_amount, sizeof(bool), msg);
    if (LG_status < GrB_SUCCESS) {
        free(checked);
        fprintf(stderr, "Calloc error: (%d): file: %s, line: %d\n%s\n", LG_status,
                __FILE__, __LINE__, msg);
        return LG_status;
    };

    for (size_t i = 0; i < symbols_amount; i++) {
        checked[i] = false;
    }

    size_t capacity = 128;
    *size = 0;
    LG_status = LAGraph_Calloc((void **)symbols, capacity, sizeof(CFL_Symbol), msg);
    if (LG_status < GrB_SUCCESS) {
        free(checked);
        free(*symbols);
        fprintf(stderr, "Calloc error: (%d): file: %s, line: %d\n%s\n", LG_status,
                __FILE__, __LINE__, msg);
        return LG_status;
    }

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
                LG_status = LAGraph_Realloc((void **)symbols, capacity, capacity / 2,
                                            sizeof(CFL_Symbol), msg);
                if (LG_status < GrB_SUCCESS) {
                    free(checked);
                    free(*symbols);
                    fprintf(stderr, "Realloc error: (%d): file: %s, line: %d\n%s\n",
                            LG_status, __FILE__, __LINE__, msg);
                    return LG_status;
                }
            }

            (*symbols)[(*size)++] = sym;

            // printf("Inserted (%ld, %ld, %ld)\n", sym.index, sym.base_index, sym.count);
        }
    }

    for (size_t i = 0; i < symbols_amount; i++) {
        if (checked[i])
            continue;

        if (*size == capacity) {
            capacity *= 2;
            LG_status = LAGraph_Realloc((void **)symbols, capacity, capacity / 2,
                                        sizeof(CFL_Symbol), msg);
            if (LG_status < GrB_SUCCESS) {
                free(checked);
                free(*symbols);
                fprintf(stderr, "Realloc error: (%d): file: %s, line: %d\n%s\n",
                        LG_status, __FILE__, __LINE__, msg);
                return LG_status;
            }
        }

        CFL_Symbol sym;
        sym.base_index = i;
        sym.index = *size;
        sym.count = 0;

        (*symbols)[(*size)++] = sym;
        checked[i] = true;
        // printf("Inserted (%ld, %ld, %ld)\n", sym.index, sym.base_index, sym.count);
    }

    free(checked);
    return GrB_SUCCESS;
}

// before: Homka_i
// after: Homka_1
//        Homka_2
//        Homka_3
static GrB_Info explode_rules(const LAGraph_rule_EWCNF *rules, size_t rules_count,
                              LAGraph_rule_EWCNF **new_rules, size_t *new_rules_count,
                              char *msg) {
    size_t new_rules_capacity = 128;
    size_t new_rules_size = 0;
    int LG_status = 0;
    LG_status = LAGraph_Calloc((void **)new_rules, new_rules_capacity,
                               sizeof(LAGraph_rule_EWCNF), msg);
    if (LG_status < GrB_SUCCESS) {
        free(*new_rules);
        fprintf(stderr, "Calloc error: (%d): file: %s, line: %d\n%s\n", LG_status,
                __FILE__, __LINE__, msg);
        return LG_status;
    }

    // explode rules
    for (size_t i_rule = 0; i_rule < rules_count; i_rule++) {
        LAGraph_rule_EWCNF rule = rules[i_rule];

        if (new_rules_size == new_rules_capacity) {
            LG_status =
                LAGraph_Realloc((void **)new_rules, new_rules_capacity * 2,
                                new_rules_capacity, sizeof(LAGraph_rule_EWCNF), msg);
            if (LG_status < GrB_SUCCESS) {
                free(*new_rules);
                fprintf(stderr, "Realloc error: (%d): file: %s, line: %d\n%s\n",
                        LG_status, __FILE__, __LINE__, msg);
                return LG_status;
            }
            new_rules_capacity *= 2;
        }

        if (rule.indexed_count == 0) {
            (*new_rules)[new_rules_size++] = rule;
        } else {
            while (new_rules_size + rule.indexed_count >= new_rules_capacity) {
                LG_status =
                    LAGraph_Realloc((void **)new_rules, new_rules_capacity * 2,
                                    new_rules_capacity, sizeof(LAGraph_rule_EWCNF), msg);
                if (LG_status < GrB_SUCCESS) {
                    free(*new_rules);
                    fprintf(stderr, "Realloc error: (%d): file: %s, line: %d\n%s\n",
                            LG_status, __FILE__, __LINE__, msg);
                    return LG_status;
                }
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

// LAGraph_CFL_reachability_adv: Context-Free Language Reachability Matrix-Based Algorithm
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
                         // outputs[k]: (i, j) = true if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the symbol 'k' of the specified CFG.
                         //
                         // Note: output[t] where t is index of terminal will be just full
                         // copy of adj_matrices[t]
    // Input
    const GrB_Matrix
        *adj_matrices, // Array of adjacency matrices representing the graph.
                       // The length of this array is equal to the count of
                       // symbols (symbols_amount).
                       //
                       // adj_matrices[i]: (i, j) == 1 if and only if there
                       // is an edge between nodes i and j with the label of
                       // the symbol corresponding to index 'i' (where i is
                       // in the range [0, symbols_amount - 1]).
                       //
                       // Note: adj_matrices[N] where N is index of nonterminal doesn't
                       // used for algorithms and may be NULL
    size_t symbols_amount,           // Count of terminal and nonterminals
    const LAGraph_rule_EWCNF *rules, // The rules of the CFG.
                                     // Warning: now not ready for N -> A rules, where is
                                     // N and A are nonterminals.
    size_t rules_count,              // The total number of rules in the CFG.
    char *msg,                       // Message string for error reporting.
    int8_t optimizations             // Optimizations flags
) {
    // Declare workspace and clear the msg string, if not NULL
    CFL_Matrix **delta_matrices, **matrices, **temp_matrices;
    CFL_Matrix *iden = NULL;
    GrB_Matrix identity_matrix = NULL;

    // for OPT_BLOCK optimization
    size_t new_symbols_amount = 0;
    GrB_Matrix *new_adj_matrices = NULL;
    LAGraph_rule_EWCNF *new_rules = NULL;
    size_t new_rules_count = 0;
    CFL_Symbol *symbols = NULL;

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

    if (optimizations & OPT_BLOCK) {
        TRY(get_new_symbols(rules, rules_count, symbols_amount, &symbols,
                            &new_symbols_amount, msg));

        TRY(LAGraph_Calloc((void **)&new_adj_matrices, new_symbols_amount,
                           sizeof(GrB_Matrix), msg));

        for (size_t i = 0; i < new_symbols_amount; i++) {
            CFL_Symbol sym = symbols[i];
            if (sym.count == 0) {
                TRY(GrB_Matrix_dup(&new_adj_matrices[i], adj_matrices[sym.base_index]));
            } else {
                GrB_Matrix new_col_matrix;
                TRY(GrB_Matrix_new(&new_col_matrix, GrB_BOOL, n * sym.count, n));
                GrB_Matrix *Tiles = (GrB_Matrix *)adj_matrices + sym.base_index;
                TRY(GxB_Matrix_concat(new_col_matrix, Tiles, sym.count, 1, GrB_NULL));
                new_adj_matrices[i] = new_col_matrix;
            }
        }

        TRY(LAGraph_Calloc((void **)&new_rules, rules_count, sizeof(LAGraph_rule_EWCNF),
                           msg));
        for (size_t i = 0; i < rules_count; i++) {
            LAGraph_rule_EWCNF rule = rules[i];
            LAGraph_rule_EWCNF new_rule = rule;
            new_rule.indexed = rule.indexed;
            new_rule.indexed_count = rule.indexed_count;

            for (size_t i_sym = 0; i_sym < new_symbols_amount; i_sym++) {
                CFL_Symbol sym = symbols[i_sym];

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

            new_rules[new_rules_count++] = new_rule;
        }
    } else {
        new_symbols_amount = symbols_amount;
        new_adj_matrices = (GrB_Matrix *)adj_matrices;
        TRY(explode_rules(rules, rules_count, &new_rules, &new_rules_count, msg));
    }

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

        GrB_Matrix matrix;
        TRY(GrB_Matrix_new(&matrix, GrB_BOOL, nrows, ncols));
        if (optimizations & OPT_LAZY) {
            TRY(CFL_matrix_from_base_lazy(&matrices[i], matrix));
        } else {
            TRY(CFL_matrix_from_base(&matrices[i], matrix));
        }

        TRY(GrB_Matrix_new(&matrix, GrB_BOOL, nrows, ncols));
        TRY(CFL_matrix_from_base(&temp_matrices[i], matrix));
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
        // matrix_print_lazy(nonterm_matrix, optimizations);
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

    for (size_t i = 0; i < new_symbols_amount; i++) {
        if (optimizations & OPT_BLOCK) {
            CFL_Symbol sym = symbols[i];
            if (sym.count == 0) {
                if (matrices[sym.index]->base_matrices_count == 0) {
                    TRY(GrB_Matrix_free(&outputs[sym.base_index]));
                    TRY(GrB_Matrix_dup(&outputs[sym.base_index],
                                       matrices[sym.index]->base));
                } else {
                    // printf("RESULT");
                    CFL_Matrix *result;
                    TRY(CFL_matrix_to_base(&result, matrices[sym.index], optimizations));
                    TRY(GrB_Matrix_dup(&outputs[sym.base_index], result->base));
                    TRY(CFL_matrix_free(&result));
                }
            } else {
                GrB_Matrix matrix_to_split = NULL;
                if (matrices[sym.index]->base_matrices_count == 0) {
                    TRY(GrB_Matrix_dup(&matrix_to_split, matrices[sym.index]->base));
                } else {
                    CFL_Matrix *result;
                    TRY(CFL_matrix_to_base(&result, matrices[sym.index], optimizations));
                    TRY(GrB_Matrix_dup(&matrix_to_split, result->base));
                    TRY(CFL_matrix_free(&result));
                }

                GrB_Index *nrows;
                TRY(LAGraph_Calloc((void **)&nrows, sym.count, sizeof(GrB_Index), msg));
                for (size_t row_i = 0; row_i < sym.count; row_i++) {
                    nrows[row_i] = n;
                }
                GrB_Index *ncols;
                TRY(LAGraph_Calloc((void **)&ncols, 1, sizeof(GrB_Index), msg));
                ncols[0] = n;

                GrB_Index m = sym.count;
                GrB_Index n = 1;

                GrB_Index matrix_to_split_nrows;
                GrB_Index matrix_to_split_ncols;
                TRY(GrB_Matrix_nrows(&matrix_to_split_nrows, matrix_to_split));
                TRY(GrB_Matrix_ncols(&matrix_to_split_ncols, matrix_to_split));
                // printf("%d %d\n", matrix_to_split_nrows, matrix_to_split_ncols);
                if (matrix_to_split_ncols > matrix_to_split_nrows) {
                    GrB_Index *temp;
                    temp = nrows;
                    nrows = ncols;
                    ncols = temp;

                    GrB_Index temp_n;
                    temp_n = m;
                    m = n;
                    n = temp_n;
                }

                TRY(GxB_Matrix_split(outputs + sym.base_index, m, n, nrows, ncols,
                                     matrix_to_split, GrB_NULL));
                free(nrows);
                free(ncols);
                TRY(GrB_free(&matrix_to_split));
            }
        } else {
            if (matrices[i]->base_matrices_count == 0) {
                TRY(GrB_Matrix_dup(&outputs[i], matrices[i]->base));
            } else {
                CFL_Matrix *result;
                TRY(CFL_matrix_to_base(&result, matrices[i], optimizations));
                TRY(GrB_Matrix_dup(&outputs[i], result->base));
                TRY(CFL_matrix_free(&result));
            }
            // outputs[i] = matrices[i].base;
        }
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}
