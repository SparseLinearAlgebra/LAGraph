#include "LG_internal.h"
#include <LAGraphX.h>

// This function groups rules of a context-free grammar in Weak Chomsky Normal Form
// (represented by LAGraph_rule_WCNF structure) by their type:
//   - epsilon rules (prod_A == -1 && prod_B == -1)
//   - terminal rules (prod_A != -1 && prod_B == -1)
//   - binary rules   (prod_A != -1 && prod_B != -1)
//
// It also counts the number of rules of each type.
//
// Important: This function does not validate the rules.
// It assumes prior validation (e.g., by LG_CFL_CHECK_GRAMMAR_INPUTS macro).
// The memory allocated by the calling function.

void LAGraph_CFL_classify_rules(
    // Output
    size_t *eps_rules,  // Array of eps-rules
    size_t *eps_count,  // Number of terminal rules
    size_t *term_rules, // Array of term-rules
    size_t *term_count, // Number of terminal rules
    size_t *bin_rules,  // Array of bin-rules
    size_t *bin_count,   // Number of terminal rules
    // Input
    const LAGraph_rule_WCNF *rules, // Array of all rules
    int64_t rules_count // Number of terminal rules
)
{
    for (int64_t i = 0; i < rules_count; i++)
    {
        LAGraph_rule_WCNF rule = rules[i];

        bool is_rule_eps = rule.prod_A == -1 && rule.prod_B == -1;
        bool is_rule_term = rule.prod_A != -1 && rule.prod_B == -1;
        bool is_rule_bin = rule.prod_A != -1 && rule.prod_B != -1;

        // [Variable -> eps]
        if (is_rule_eps)
        {
            eps_rules[(*eps_count)++] = i;
            continue;
        }

        // [Variable -> term]
        if (is_rule_term)
        {
            term_rules[(*term_count)++] = i;
            continue;
        }

        // [Variable -> A B]
        if (is_rule_bin)
        {
            bin_rules[(*bin_count)++] = i;
            continue;
        }
    }
}