#include <LAGraphX.h>
#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdio.h>

typedef struct {
    GrB_Index state_count;
    GrB_Index terminal_count;
    GrB_Index nonterminal_count;
    GrB_Index start_nonterminal;
    GrB_Matrix* terminal_matrices;
    GrB_Matrix* nonterminal_matrices;
    GrB_Index* start_states;
    GrB_Vector* final_states;
} RSM;

#define LG_FREE_WORK                                                        \
    {                                                                       \
        GrB_free(&M3);                                                      \
        LAGraph_Free((void**) &i, NULL);                                    \
        LAGraph_Free((void**) &j, NULL);                                    \
        LAGraph_Free((void**) &V, NULL);                                    \
    }

GrB_Info transitive_closure_inplace(GrB_Matrix A) {
    GrB_Index nvals_old = 0, nvals_new = 0;
    GrB_Matrix_nvals(&nvals_old, A);
    while (true) {
        GrB_mxm(A, NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, A, A, NULL);
        GrB_Matrix_nvals(&nvals_new, A);
        if (nvals_new == nvals_old) break;
        nvals_old = nvals_new;
    }
    return GrB_SUCCESS;
}

GrB_Info LAGraph_CFL_AllPaths_Kronecker 
(
    GrB_Matrix *outputs,         
    const GrB_Matrix *adj_matrices,  
    int64_t terms_count,
    RSM *rsm,                    
    char *msg
)
{
    GrB_Matrix M3 = NULL;
    GrB_Index *i = NULL, *j = NULL;
    bool *V = NULL;
    
    GrB_Index g_dim;
    GrB_Matrix_ncols(&g_dim, adj_matrices[0]);
    GrB_Index r_dim = rsm->state_count;
    GrB_Index kronecker_dim = r_dim * g_dim;

    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Matrix_new(&outputs[nt], GrB_BOOL, g_dim, g_dim);
        
        GrB_Index s = rsm->start_states[nt];
        bool is_final = false;
        GrB_Vector_extractElement_BOOL(&is_final, rsm->final_states[nt], s);
        
        if (is_final) {
            for (GrB_Index v = 0; v < g_dim; ++v) {
                GrB_Matrix_setElement_BOOL(outputs[nt], true, v, v);
            }
        }
    }

    GrB_Matrix_new(&M3, GrB_BOOL, kronecker_dim, kronecker_dim);

    GrB_Index nvals_old = 0, nvals_new = 0;
    bool changed = true;

    while (changed) {
        nvals_old = 0;
        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_Index vals;
            GrB_Matrix_nvals(&vals, outputs[nt]);
            nvals_old += vals;
        }

        GrB_Matrix_clear(M3);

        for (int64_t t = 0; t < terms_count; ++t) {
            GrB_kronecker(M3, NULL, GrB_LOR, GrB_LAND,
                          rsm->terminal_matrices[t], adj_matrices[t], NULL);
        }

        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_kronecker(M3, NULL, GrB_LOR, GrB_LAND,
                          rsm->nonterminal_matrices[nt], outputs[nt], NULL);
        }

        transitive_closure_inplace(M3);

        GrB_Index m3_nvals;
        GrB_Matrix_nvals(&m3_nvals, M3);

        if (m3_nvals > 0) {
            i = (GrB_Index*)malloc(m3_nvals * sizeof(GrB_Index));
            j = (GrB_Index*)malloc(m3_nvals * sizeof(GrB_Index));
            V = (bool*)malloc(m3_nvals * sizeof(bool));
            
            GrB_Matrix_extractTuples_BOOL(i, j, V, &m3_nvals, M3);

            for (GrB_Index k = 0; k < m3_nvals; ++k) {
                GrB_Index s = i[k] / g_dim; 
                GrB_Index f = j[k] / g_dim; 
                GrB_Index x = i[k] % g_dim; 
                GrB_Index y = j[k] % g_dim; 

                for (GrB_Index nt = 0; nt < rsm->nonterminal_count; ++nt) {
                    if (rsm->start_states[nt] == s) {
                        bool is_final = false;
                        GrB_Vector_extractElement_BOOL(&is_final, rsm->final_states[nt], f);
                        if (is_final) {
                            GrB_Matrix_setElement_BOOL(outputs[nt], true, x, y);
                        }
                    }
                }
            }
            free(i); free(j); free(V);
            i = NULL; j = NULL; V = NULL;
        }

        nvals_new = 0;
        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_Index vals;
            GrB_Matrix_nvals(&vals, outputs[nt]);
            nvals_new += vals;
        }
        changed = (nvals_new != nvals_old);
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
}

