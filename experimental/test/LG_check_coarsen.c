//------------------------------------------------------------------------------
// LG_check_coarsen: construct the kcore of a graph (simple method)
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

// A naive coarsening implementation that individually traverses edges in the original graph 
// and updates the corresponding edge in the coarsened graph

#include "LG_internal.h"
#include "LG_test.h"
#include "LG_test.h"
#include "LG_Xtest.h"

// #undef LG_FREE_ALL
// #undef LG_FREE_WORK

int LG_check_coarsen
(
    // outputs:
    GrB_Matrix *coarsened,    // coarsened adjacency
    // inputs:
    GrB_Matrix A,          // input adjacency (for the purposes of testing, is FP64)
    GrB_Vector parent,     // parent mapping
    GrB_Vector mapping,    // new labels of nodes, used to populate resulting adjacency matrix, NULL if preserve_mapping = 1
    int preserve_mapping,  // whether to preserve the original namespace of nodes
    int combine_weights,   // whether to combine the weights of edges that collapse together
    char *msg
)
{
    GrB_Matrix result = NULL ;

    GrB_Index n ;
    GrB_Index nvals ;

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;

    GrB_Index n_new = n ;
    
    // check that parent mapping is compressed (i.e. p[p[i]] = p[i] for all i)
    // also calculate n_new
    for (GrB_Index i = 0 ; i < n ; i++) {
        uint64_t par ;
        GRB_TRY (GrB_Vector_extractElement (&par, parent, i)) ;
        // TODO: What error code to use here?
        LG_ASSERT (par >= 0 && par < n, -1) ;
        uint64_t grandpa ;
        GRB_TRY (GrB_Vector_extractElement (&grandpa, parent, par)) ;
        LG_ASSERT (grandpa == par, -1) ;

        if (par != i) {
            n_new-- ;
        }
    }

    // check that mapping vector is valid
    for (GrB_Index i = 0 ; i < n ; i++) {
        uint64_t new_label ;
        GRB_TRY (GrB_Vector_extractElement (&new_label, mapping, i)) ;
        // TODO: What error code to use here?
        LG_ASSERT (new_label >= 0 && new_label < n_new, -1) ;
    }


    GRB_TRY (GrB_Matrix_new (&result, GrB_FP64, n_new, n_new)) ;

    GrB_Index *rows = NULL, *cols = NULL ;
    double *vals = NULL ;

    LG_TRY (LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(double), msg)) ;

    GRB_TRY (GrB_Matrix_extractTuples (rows, cols, vals, &nvals, A)) ;

    for (GrB_Index i = 0 ; i < nvals ; i++) {
        GrB_Index u = rows [i] ;
        GrB_Index v = cols [i] ;

        if (u > v) {
            // only consider upper-triangular entries (don't duplicate edges)
            continue ;
        }

        uint64_t u_par, v_par ;

        GRB_TRY (GrB_Vector_extractElement (&u_par, parent, u)) ;
        GRB_TRY (GrB_Vector_extractElement (&v_par, parent, v)) ;

        if (u_par == v_par) {
            // both nodes part of same coarsened component
            continue ;
        }
        uint64_t u_par_newlabel, v_par_newlabel ;
        if (preserve_mapping) {
            // new labels are the same as old labels
            u_par_newlabel = u_par ;
            v_par_newlabel = v_par ;
        } else {
            // find new labels
            GRB_TRY (GrB_Vector_extractElement (&u_par_newlabel, mapping, u_par)) ;
            GRB_TRY (GrB_Vector_extractElement (&v_par_newlabel, mapping, v_par)) ;
        }
        double res_weight = 1 ;

        if (combine_weights) {
            // current weight of edge between u_par and v_par
            double curr_weight ;
            GRB_TRY (GrB_Matrix_extractElement (&curr_weight, result, u_par_newlabel, v_par_newlabel)) ;

            // weight of the current edge being added
            double my_weight = vals [i] ;
            res_weight = curr_weight + my_weight ;
        }

        GRB_TRY (GrB_Matrix_setElement (result, res_weight, u_par_newlabel, v_par_newlabel)) ;
        GRB_TRY (GrB_Matrix_setElement (result, res_weight, v_par_newlabel, u_par_newlabel)) ;
    }
    (*coarsened) = result ;

    LG_TRY (LAGraph_Free ((void**)(&rows), msg)) ;
    LG_TRY (LAGraph_Free ((void**)(&cols), msg)) ;
    LG_TRY (LAGraph_Free ((void**)(&vals), msg)) ;

    return (GrB_SUCCESS) ;
}