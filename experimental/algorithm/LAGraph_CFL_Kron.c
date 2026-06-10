#include <LAGraphX.h>
#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdbool.h>
#include <stdlib.h>

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
        GrB_free(&CombinedGraph);                                           \
    }

GrB_Info transitive_closure_inplace(GrB_Matrix A) {
    GrB_Index n;
    GrB_Matrix_ncols(&n, A);
    
    GrB_Matrix Frontier = NULL;
    GrB_Matrix NextFrontier = NULL;
    
    GrB_Matrix_dup(&Frontier, A);
    GrB_Matrix_new(&NextFrontier, GrB_BOOL, n, n);
    
    while (true) {
        GrB_Matrix_clear(NextFrontier);
        
        GrB_mxm(NextFrontier, NULL, NULL, GxB_ANY_PAIR_BOOL, Frontier, A, NULL);
        
        GrB_assign(Frontier, A, NULL, NextFrontier, GrB_ALL, n, GrB_ALL, n, GrB_DESC_RC);
        
        GrB_Index frontier_nvals;
        GrB_Matrix_nvals(&frontier_nvals, Frontier);
        if (frontier_nvals == 0) break; 
        
        GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, Frontier, NULL);
    }
    
    GrB_free(&Frontier);
    GrB_free(&NextFrontier);
    return GrB_SUCCESS;
}

void dfs_find_paths(Node** adj_list, GrB_Index current, GrB_Index target,
                    GrB_Index* path, GrB_Index depth, GrB_Index g_dim,
                    bool* visited)
{
    path[depth] = current;
    visited[current] = true;

    if (current == target) {
        printf("Найден путь: ");
        for (GrB_Index p = 0; p <= depth; p++) {
            printf("%llu ", path[p] % g_dim);
            if (p < depth) printf("-> ");
        }
        printf("\n");
    } else {
        Node* neighbor = adj_list[current];
        while (neighbor != NULL) {
            if (!visited[neighbor->vertex]) {
                dfs_find_paths(adj_list, neighbor->vertex, target, path, depth + 1, g_dim, visited);
            }
            neighbor = neighbor->next;
        }
    }

    visited[current] = false; 
}

void PrintAllPaths(GrB_Matrix M_intersect, RSM* rsm, GrB_Index nt_index, GrB_Index start_v, GrB_Index target_v, GrB_Index g_dim) {
    GrB_Index kronecker_dim;
    GrB_Matrix_ncols(&kronecker_dim, M_intersect);

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, M_intersect);
    
    GrB_Index *i = (GrB_Index*)malloc(nvals * sizeof(GrB_Index));
    GrB_Index *j = (GrB_Index*)malloc(nvals * sizeof(GrB_Index));
    bool *V = (bool*)malloc(nvals * sizeof(bool));
    GrB_Matrix_extractTuples_BOOL(i, j, V, &nvals, M_intersect);

    Node** adj_list = (Node**)calloc(kronecker_dim, sizeof(Node*));
    for (GrB_Index k = 0; k < nvals; k++) {
        Node* new_node = (Node*)malloc(sizeof(Node));
        new_node->vertex = j[k];
        new_node->next = adj_list[i[k]];
        adj_list[i[k]] = new_node;
    }

    GrB_Index s_state = rsm->start_states[nt_index];
    GrB_Index kronecker_start = (s_state * g_dim) + start_v;

    for (GrB_Index f_state = 0; f_state < rsm->state_count; f_state++) {
        bool is_final = false;
        GrB_Vector_extractElement_BOOL(&is_final, rsm->final_states[nt_index], f_state);
        
        if (is_final) {
            GrB_Index kronecker_target = (f_state * g_dim) + target_v;
            
            GrB_Index* path = (GrB_Index*)malloc(kronecker_dim * sizeof(GrB_Index));
            bool* visited = (bool*)calloc(kronecker_dim, sizeof(bool));
            
            printf("Ищем пути из вершины %llu в %llu для нетерминала %llu (состояния %llu -> %llu)\n",
                   start_v, target_v, nt_index, s_state, f_state);
            
            dfs_find_paths(adj_list, kronecker_start, kronecker_target, path, 0, g_dim, visited);
            
            free(path);
            free(visited);
        }
    }

    for (GrB_Index idx = 0; idx < kronecker_dim; idx++) {
        Node* curr = adj_list[idx];
        while (curr != NULL) {
            Node* temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
    free(adj_list); free(i); free(j); free(V);
}

GrB_Info LAGraph_CFL_AllPaths_Kronecker_orig
(
    GrB_Matrix *M_intersect_out,
    GrB_Matrix *outputs,
    const GrB_Matrix *adj_matrices,
    int64_t terms_count,
    RSM *rsm,
    char *msg
)
{
    GrB_Index g_dim;
    GrB_Matrix_ncols(&g_dim, adj_matrices[0]);
    GrB_Index r_dim = rsm->state_count;
    GrB_Index kronecker_dim = r_dim * g_dim;

    GrB_Matrix *delta_outputs = (GrB_Matrix *)malloc(rsm->nonterminal_count * sizeof(GrB_Matrix));
    GrB_Matrix temp_output;
    GrB_Matrix_new(&temp_output, GrB_BOOL, g_dim, g_dim);

    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Matrix_new(&outputs[nt], GrB_BOOL, g_dim, g_dim);
        GrB_Matrix_new(&delta_outputs[nt], GrB_BOOL, g_dim, g_dim);

        GrB_Index s = rsm->start_states[nt];
        bool is_final = false;
        GrB_Vector_extractElement_BOOL(&is_final, rsm->final_states[nt], s);

        if (is_final) {
            for (GrB_Index v = 0; v < g_dim; ++v) {
                GrB_Matrix_setElement_BOOL(outputs[nt], true, v, v);
                GrB_Matrix_setElement_BOOL(delta_outputs[nt], true, v, v);
            }
        }
    }

    GrB_Matrix M_graph = NULL;
    GrB_Matrix CombinedGraph = NULL; 
    GrB_Matrix DeltaM = NULL;

    GrB_Matrix_new(&M_graph, GrB_BOOL, kronecker_dim, kronecker_dim);

    for (int64_t t = 0; t < terms_count; ++t) {
        GrB_kronecker(M_graph, NULL, GrB_LOR, GxB_PAIR_BOOL,
                      rsm->terminal_matrices[t], adj_matrices[t], NULL);
    }

    GrB_Matrix_new(&DeltaM, GrB_BOOL, kronecker_dim, kronecker_dim);

    GrB_Index max_finals = 0;
    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Index f_nvals = 0;
        GrB_Vector_nvals(&f_nvals, rsm->final_states[nt]);
        if (f_nvals > max_finals) max_finals = f_nvals;
    }

    GrB_Index *row_indices = (GrB_Index *)malloc(g_dim * sizeof(GrB_Index));
    GrB_Index *col_indices = (GrB_Index *)malloc(g_dim * sizeof(GrB_Index));
    GrB_Index *f_indices   = (GrB_Index *)malloc((max_finals > 0 ? max_finals : 1) * sizeof(GrB_Index));
    bool *f_values         = (bool *)malloc((max_finals > 0 ? max_finals : 1) * sizeof(bool));

    bool changed = true;
    bool is_first_iteration = true;

    while (changed) {
        changed = false;
        GrB_Matrix_clear(DeltaM);
        
        bool has_new_edges = false;
        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_Index delta_vals;
            GrB_Matrix_nvals(&delta_vals, delta_outputs[nt]);
            if (delta_vals > 0) {
                has_new_edges = true;
                GrB_kronecker(DeltaM, NULL, GrB_LOR, GxB_PAIR_BOOL,
                              rsm->nonterminal_matrices[nt], delta_outputs[nt], NULL);
            }
        }

        if (!has_new_edges && !is_first_iteration) break;
      
        if (has_new_edges) {
            GrB_eWiseAdd(M_graph, NULL, NULL, GrB_LOR, M_graph, DeltaM, NULL);
        }

        if (CombinedGraph != NULL) {
            GrB_free(&CombinedGraph);
        }
        GrB_Matrix_dup(&CombinedGraph, M_graph);
        
        transitive_closure_inplace(CombinedGraph);

        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_Matrix_clear(temp_output);
            GrB_Index s = rsm->start_states[nt];

            for (GrB_Index i = 0; i < g_dim; ++i) {
                row_indices[i] = s * g_dim + i;
            }

            GrB_Index f_nvals = max_finals;
            GrB_Vector_extractTuples_BOOL(f_indices, f_values, &f_nvals, rsm->final_states[nt]);
            for (GrB_Index k = 0; k < f_nvals; ++k) {
                if (f_values[k]) {
                    GrB_Index f = f_indices[k];
                    for (GrB_Index j = 0; j < g_dim; ++j) {
                        col_indices[j] = f * g_dim + j;
                    }
                    GrB_Matrix_extract(temp_output, NULL, GrB_LOR, CombinedGraph,
                                       row_indices, g_dim, col_indices, g_dim, NULL);
                }
            }

            GrB_Matrix_clear(delta_outputs[nt]);
            GrB_assign(delta_outputs[nt], outputs[nt], NULL, temp_output, 
                       GrB_ALL, g_dim, GrB_ALL, g_dim, GrB_DESC_RC);

            GrB_Index new_vals;
            GrB_Matrix_nvals(&new_vals, delta_outputs[nt]);
            
            if (new_vals > 0) {
                changed = true;
                GrB_eWiseAdd(outputs[nt], NULL, NULL, GrB_LOR, outputs[nt], delta_outputs[nt], NULL);
            }
        }
        is_first_iteration = false; 
    }

    free(row_indices); free(col_indices);
    free(f_indices); free(f_values);
    
    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_free(&delta_outputs[nt]);
    }
    free(delta_outputs);
    
    GrB_free(&temp_output);

    GrB_Matrix_dup(M_intersect_out, M_graph);
    if (M_graph != NULL) GrB_free(&M_graph);
    if (CombinedGraph != NULL) GrB_free(&CombinedGraph); 
    if (DeltaM != NULL) GrB_free(&DeltaM);

    return GrB_SUCCESS;
}

GrB_Info LAGraph_CFL_AllPaths_Kronecker
(
    GrB_Matrix *M_intersect_out,
    GrB_Matrix *outputs,
    const GrB_Matrix *adj_matrices,
    int64_t terms_count,
    RSM *rsm,
    char *msg
)
{
    GrB_Index g_dim;
    GrB_Matrix_ncols(&g_dim, adj_matrices[0]);
    GrB_Index r_dim = rsm->state_count;
    GrB_Index kronecker_dim = r_dim * g_dim;

    GrB_Matrix M_work = NULL;
    GrB_Matrix_new(&M_work, GrB_BOOL, kronecker_dim, kronecker_dim);

    GrB_Matrix *delta_outputs = (GrB_Matrix *)malloc(rsm->nonterminal_count * sizeof(GrB_Matrix));
    GrB_Matrix temp_output;
    GrB_Matrix_new(&temp_output, GrB_BOOL, g_dim, g_dim);

    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Matrix_new(&outputs[nt], GrB_BOOL, g_dim, g_dim);
        GrB_Matrix_new(&delta_outputs[nt], GrB_BOOL, g_dim, g_dim);

        GrB_Index s = rsm->start_states[nt];
        bool is_final = false;
        GrB_Vector_extractElement_BOOL(&is_final, rsm->final_states[nt], s);

        if (is_final) {
            for (GrB_Index v = 0; v < g_dim; ++v) {
                GrB_Matrix_setElement_BOOL(outputs[nt], true, v, v);
            }
            GrB_kronecker(M_work, NULL, GrB_LOR, GxB_PAIR_BOOL,
                          rsm->nonterminal_matrices[nt], outputs[nt], NULL);
        }
    }

    for (int64_t t = 0; t < terms_count; ++t) {
        GrB_kronecker(M_work, NULL, GrB_LOR, GxB_PAIR_BOOL,
                      rsm->terminal_matrices[t], adj_matrices[t], NULL);
    }

    GrB_Matrix S_diag;
    GrB_Matrix_new(&S_diag, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Index num_start_nodes = rsm->nonterminal_count * g_dim;
    GrB_Index *s_row = (GrB_Index *)malloc(num_start_nodes * sizeof(GrB_Index));
    GrB_Index *s_col = (GrB_Index *)malloc(num_start_nodes * sizeof(GrB_Index));
    bool *s_val = (bool *)malloc(num_start_nodes * sizeof(bool));

    GrB_Index idx = 0;
    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Index s = rsm->start_states[nt];
        for (GrB_Index i = 0; i < g_dim; ++i) {
            s_row[idx] = s * g_dim + i;
            s_col[idx] = s * g_dim + i;
            s_val[idx] = true;
            idx++;
        }
    }
    GrB_Matrix_build_BOOL(S_diag, s_row, s_col, s_val, num_start_nodes, GrB_LOR);
    free(s_row); free(s_col); free(s_val);

    GrB_Matrix Reach, Frontier, NextFrontier, DeltaM, Reach_x_Delta, TempKron;
    GrB_Matrix_new(&Reach, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Matrix_new(&Frontier, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Matrix_new(&NextFrontier, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Matrix_new(&DeltaM, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Matrix_new(&Reach_x_Delta, GrB_BOOL, kronecker_dim, kronecker_dim);
    GrB_Matrix_new(&TempKron, GrB_BOOL, kronecker_dim, kronecker_dim);

    GrB_Matrix_dup(&Frontier, S_diag);
    GrB_Matrix_dup(&Reach, S_diag);

    GrB_Index *row_indices = (GrB_Index *)malloc(g_dim * sizeof(GrB_Index));
    GrB_Index *col_indices = (GrB_Index *)malloc(g_dim * sizeof(GrB_Index));
    
    GrB_Index max_finals = 0;
    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_Index f_nvals = 0;
        GrB_Vector_nvals(&f_nvals, rsm->final_states[nt]);
        if (f_nvals > max_finals) max_finals = f_nvals;
    }
    GrB_Index *f_indices = (GrB_Index *)malloc((max_finals > 0 ? max_finals : 1) * sizeof(GrB_Index));
    bool *f_values = (bool *)malloc((max_finals > 0 ? max_finals : 1) * sizeof(bool));

    while (true) {
        GrB_Index f_nvals;
        GrB_Matrix_nvals(&f_nvals, Frontier);
        if (f_nvals == 0) break;

        GrB_Matrix_clear(DeltaM);
        bool has_new_edges = false;

        for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
            GrB_Matrix_clear(temp_output);
            GrB_Index s = rsm->start_states[nt];
            for (GrB_Index v = 0; v < g_dim; v++) row_indices[v] = s * g_dim + v;

            GrB_Index curr_finals = max_finals;
            GrB_Vector_extractTuples_BOOL(f_indices, f_values, &curr_finals, rsm->final_states[nt]);

            for (GrB_Index k = 0; k < curr_finals; ++k) {
                if (f_values[k]) {
                    GrB_Index f = f_indices[k];
                    for (GrB_Index j = 0; j < g_dim; ++j) col_indices[j] = f * g_dim + j;
                    
                    GrB_Matrix_extract(temp_output, NULL, GrB_LOR, Frontier,
                                       row_indices, g_dim, col_indices, g_dim, NULL);
                }
            }

            GrB_Matrix_clear(delta_outputs[nt]);
            GrB_assign(delta_outputs[nt], outputs[nt], NULL, temp_output,
                       GrB_ALL, g_dim, GrB_ALL, g_dim, GrB_DESC_RC);

            GrB_Index new_vals;
            GrB_Matrix_nvals(&new_vals, delta_outputs[nt]);
            
            if (new_vals > 0) {
                has_new_edges = true;
                GrB_eWiseAdd(outputs[nt], NULL, NULL, GrB_LOR, outputs[nt], delta_outputs[nt], NULL);
                
                GrB_Matrix_clear(TempKron);
                GrB_kronecker(TempKron, NULL, NULL, GxB_PAIR_BOOL,
                              rsm->nonterminal_matrices[nt], delta_outputs[nt], NULL);
                GrB_eWiseAdd(DeltaM, NULL, NULL, GrB_LOR, DeltaM, TempKron, NULL);
            }
        }

        GrB_Matrix_clear(Reach_x_Delta);

        if (has_new_edges) {
            GrB_eWiseAdd(M_work, NULL, NULL, GrB_LOR, M_work, DeltaM, NULL);

            GrB_mxm(Reach_x_Delta, Reach, NULL, GxB_ANY_PAIR_BOOL, Reach, DeltaM, GrB_DESC_RC);
            
            GrB_eWiseAdd(Reach, NULL, NULL, GrB_LOR, Reach, Reach_x_Delta, NULL);
        }

        GrB_Matrix_clear(NextFrontier);
        GrB_mxm(NextFrontier, Reach, NULL, GxB_ANY_PAIR_BOOL, Frontier, M_work, GrB_DESC_RC);
        GrB_eWiseAdd(Reach, NULL, NULL, GrB_LOR, Reach, NextFrontier, NULL);

        if (has_new_edges) {
            GrB_eWiseAdd(NextFrontier, NULL, NULL, GrB_LOR, NextFrontier, Reach_x_Delta, NULL);
        }

        GrB_Matrix temp = Frontier;
        Frontier = NextFrontier;
        NextFrontier = temp;
    }

    GrB_Matrix_dup(M_intersect_out, M_work);

    free(row_indices); free(col_indices);
    free(f_indices); free(f_values);
    for (int64_t nt = 0; nt < rsm->nonterminal_count; ++nt) {
        GrB_free(&delta_outputs[nt]);
    }
    free(delta_outputs);
    GrB_free(&temp_output); GrB_free(&M_work);
    GrB_free(&Reach); GrB_free(&Frontier);
    GrB_free(&NextFrontier); GrB_free(&DeltaM);
    GrB_free(&S_diag); GrB_free(&Reach_x_Delta); GrB_free(&TempKron);

    return GrB_SUCCESS;
}
